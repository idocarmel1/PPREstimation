"""Export one Excel workbook per Ecopath model, holding every SPPR method side by side.

For each model JSON in a directory this builds a workbook with:

    groups_df       -- the per-group parameters the SPPR / PPR / NPP calculations consume
    sppr_PP         -- SPPR per group x method, summed over primary-producer sources only
    sppr_inner      -- the same, over within-system sources (primary producers + detritus)
    sppr_all        -- the same, over every basal source (adds Import)
    model_health    -- diagnose_sppr output, one row per TE_option, under a fixed detritus config
    footprint       -- PPR / NPP footprint per method
    mc_diagnostics  -- Monte-Carlo accept/reject breakdown per MC method
    run_notes       -- the conventions and caveats needed to read the numbers correctly

and a plain-text run report next to them listing every warning, skip and failure.

`read_pprs_excel` reads a workbook back into the same tables, index and column structure.

Conventions used throughout, and repeated in the run_notes sheet:
  * NaN means "not available", never zero. A source a method does not resolve is NaN; a method
    that raised leaves its whole block NaN.
  * The three SPPR sheets follow get_PPR's own source semantics: all sources, all minus Import,
    and all minus Import and Detritus respectively. Each is a flat groups x methods table with
    no MultiIndex on either axis.
  * Methods that return a single un-attributed SPPR column (SPPR_1986, SPPR_1995) are NaN in
    sppr_PP and in ppr_pp_only, because their value cannot honestly be attributed to primary
    production alone.
  * Every method runs in a worker process under a wall-clock budget (3 minutes by default).
    A method that overruns it is killed, left NaN and reported as a timeout rather than as a
    failure -- the run_notes 'status' column and the run report both say which.

Requires `openpyxl` (the .xlsx engine pandas uses here); everything else is already a project
dependency.

Entry point:
    python create_PPRS_excel.py [json_dir] [out_dir] [--resume] [--timeout SECONDS]

Defaults to real_models/EwE_jsons -> output/, with a 180 s per-method budget
(--timeout none runs every method in-process with no budget).
"""
from __future__ import annotations

import multiprocessing as mp
import os
import sys
import traceback
import warnings
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Callable, Optional, Sequence

import numpy as np
import pandas as pd
from tqdm import tqdm

from ModelData import ModelData
from PPRCalculator import PPRCalculator

# --------------------------------------------------------------------------- constants

# The three SPPR sheets and the source sum each one holds. Every sheet is a flat
# groups x methods table: no MultiIndex anywhere.
SPPR_SHEET_SUMS = {"pp": "SUM_PP", "inner": "SUM_INNER", "all": "SUM_ALL"}

# diagnose_sppr is always run under this exact detritus configuration, so the health verdicts are
# comparable across models and TE options. det_theta is inert while det_open_mode='none' (it does
# not move the recycling gain b at all) -- it is pinned here only to make the run reproducible.
HEALTH_TE_OPTIONS = ("GE", "TE", "With Egestion")
HEALTH_DET_CONFIG = {
    "det_open_mode": "none",
    "det_theta": 1.0,
    "det_external_sppr": 0.0,
    "det_collapse_mode": "never",
}

FOOTPRINT_COLUMNS = ("ppr_all", "ppr_inner", "ppr_pp_only", "npp", "ppr2npp", "ppr2npp_pp_only")

SHEET_GROUPS = "groups_df"
SHEET_SPPR_PP = "sppr_PP"
SHEET_SPPR_INNER = "sppr_inner"
SHEET_SPPR_ALL = "sppr_all"
SHEET_HEALTH = "model_health"
SHEET_FOOTPRINT = "footprint"
SHEET_MC = "mc_diagnostics"
SHEET_NOTES = "run_notes"

LOG_FILENAME = "sppr_export_report.txt"
DEFAULT_JSON_DIR = os.path.join("real_models", "EwE_jsons")
DEFAULT_OUT_DIR = "output/EwE_jsons/models"
DEFAULT_MC_SAMPLES = 100

# Wall-clock budget for one SPPR method (and for one diagnose_sppr row). SPPR_EwE enumerates
# every path in the food web and the Monte-Carlo drivers solve the whole system once per draw;
# on a large web either can run for hours, and neither can be interrupted from inside the
# interpreter once it is down in numpy/sympy. So each one runs in a worker process that is killed
# outright when the budget expires, and the export moves on to the next method.
DEFAULT_METHOD_TIMEOUT = 180.0

# A restarted worker has to be handed the model again before it can run anything. That handshake
# pays for a process spawn plus the numpy/pandas/sympy imports, so it gets its own, far more
# generous budget -- a slow spawn must never be mistaken for a hung method.
MODEL_HANDSHAKE_TIMEOUT = 300.0

RUN_NOTES = [
    ("NaN convention",
     "NaN means not available, never zero. A basal source a method does not resolve is NaN, and "
     "a method that raised leaves its entire block NaN."),
    ("The three SPPR sheets",
     "sppr_all, sppr_inner and sppr_PP are the same table under three source scopes, mirroring "
     "get_PPR: sppr_all sums every basal source, sppr_inner drops Import, and sppr_PP drops "
     "Import and Detritus. Each is a flat groups x methods table -- rows are group seq with the "
     "group name in the first column, one column per method, no MultiIndex. Sums are taken over "
     "each method's own basal-source columns before aggregation; the per-source breakdown itself "
     "is not written to the workbook."),
    ("Un-attributed methods are blank in sppr_PP",
     "SPPR_1986 and SPPR_1995 return one un-attributed SPPR column, so their value cannot be "
     "split into a PP-only part: they are NaN in sppr_PP (and in ppr_pp_only) while sppr_inner "
     "and sppr_all carry the total. A method that resolves no Detritus source (SPPR_2015) has "
     "sppr_inner equal to sppr_PP by construction."),
    ("SPPR_2015 takes no arguments",
     "SPPR_2015() has no only_pp_det parameter (only_pp belongs to get_PPR / get_PPR2NPP_ratio). "
     "It is called bare here and PP-only filtering is applied downstream at the get_PPR layer, "
     "which is what the ppr_pp_only column and SUM_PP already do."),
    ("det_theta is inert here",
     "Under det_open_mode='none' det_theta does not affect the recycling gain b at all (verified: "
     "identical b at theta = 1, 3 and 10). Do not read the health sheet as a theta sensitivity."),
    ("Monte-Carlo rejections",
     "MC runs use exclude_diverged=True, so draws graded FAIL on divergence are dropped as well "
     "as draws with a negative SPPR. The two are often the same draws rather than additive -- see "
     "the mc_diagnostics sheet for the split by cause."),
    ("Solver progress output is suppressed",
     "SPPR_EwE and the Monte-Carlo methods are called with silent=True rather than silent=False: "
     "PPRCalculator imports tqdm.notebook, whose bar needs ipywidgets and raises "
     "ImportError('IProgress not found') outside Jupyter, which makes those methods fail outright "
     "in a terminal run. Nothing is lost -- this exporter has its own progress bar, and the MC "
     "rejection counts appear in mc_diagnostics and in the run report."),
    ("model_health is per configuration",
     "diagnose_sppr grades a CONFIGURATION, not a model: a model can be healthy under one "
     "TE_option and divergent under another. Read each row with its config_ columns."),
    ("Per-method timeout",
     "Every SPPR method, and every diagnose_sppr row, runs in a worker process with a "
     "wall-clock budget (180 s by default; --timeout SECONDS on the command line, "
     "method_timeout= from Python, and None or 0 to disable it). A method still running when "
     "the budget expires is killed and left NaN, exactly like a method that raised -- the "
     "status column of this sheet and the run report are the only places that say which of the "
     "two happened. A timeout is a statement about this machine and this model, not about the "
     "method: rerun it with a larger budget before reading anything into it."),
]


# --------------------------------------------------------------------------- method registry

@dataclass(frozen=True)
class MethodSpec:
    """One SPPR method as it is run for the report.

    Attributes:
        key: the short indicative name used as the method label everywhere in the workbook.
        description: what the method is, for the notes / log.
        run: callable(model, mc_samples) -> (sppr DataFrame, extras dict).
        source_resolved: False for methods that return a single un-attributed SPPR column.
    """
    key: str
    description: str
    run: Callable[..., tuple]
    source_resolved: bool = True


def _first(result):
    """SPPR getter for methods that return (SPPR, ...) tuples."""
    return result[0]


def _spec_1986() -> MethodSpec:
    return MethodSpec(
        "SPPR_1986", "Pauly & Christensen (1986): one catch-weighted TL, fixed TE=0.1.",
        lambda m, mc_samples: (m.SPPR_1986(), {}), source_resolved=False)


def _spec_1995(key: str, global_TE) -> MethodSpec:
    return MethodSpec(
        key, f"Christensen & Pauly (1995): per-group TE^(1-TL) with global_TE={global_TE!r}.",
        lambda m, mc_samples, te=global_TE: (m.SPPR_1995(global_TE=te), {}),
        source_resolved=False)


def _spec_ulanowicz(key: str, TE_option, global_TE, description) -> MethodSpec:
    return MethodSpec(
        key, description,
        lambda m, mc_samples, t=TE_option, g=global_TE: (
            _first(m.SPPR_EwE_Ulanowicz(TE_option=t, global_TE=g, use_EE=False)), {}))


def _spec_new(key: str, TE_option, description, **extra) -> MethodSpec:
    return MethodSpec(
        key, description,
        lambda m, mc_samples, t=TE_option, e=extra: (
            _first(m.SPPR_new(TE_option=t, **e)), {}))


def _spec_symbolic(TE_option: str, diet_import_option: str) -> MethodSpec:
    label = {"With Egestion": "WithEgestion"}.get(TE_option, TE_option)
    diet_label = {"as_PP": "asPP", "as_DC": "asDC"}[diet_import_option]

    def run(m, mc_samples, t=TE_option, d=diet_import_option):
        # SPPR_symbolic returns (sppr_symbolic, sppr_mat, equations, variables); the second
        # element is the numeric matrix.
        _, sppr, eqs, vars_ = m.SPPR_symbolic(TE_option=t, diet_import_option=d)
        return sppr, {"equations": eqs, "variables": vars_}

    return MethodSpec(
        f"sym_{label}_{diet_label}",
        f"Symbolic solver, TE_option={TE_option!r}, diet_import_option={diet_import_option!r}.",
        run)


def _spec_monte_carlo(key: str, TE_option: str, method_kwargs: Optional[dict],
                      description: str) -> MethodSpec:
    def run(m, mc_samples, t=TE_option, mk=method_kwargs):
        # silent=True deliberately: PPRCalculator imports tqdm.notebook, whose progress bar needs
        # ipywidgets and raises outside Jupyter, which would take the whole method down in a plain
        # terminal run. This exporter has its own progress bar and reports every rejection itself.
        mean, samples, reject_frac, _, _, diag = m.monte_carlo_SPPR(
            n_samples=mc_samples, TE_error_percent=10, TE_error_cut_percent=20,
            kind="new", TE_option=t, silent=True,
            method_kwargs=dict(mk) if mk else None,
            exclude_diverged=True, return_diagnostics=True)
        extras = {
            "mc": {
                "n_samples": mc_samples,
                "n_accepted": diag["n_accepted"],
                "n_rejected_negative": diag["n_rejected_negative"],
                "n_rejected_diverged": diag["n_rejected_diverged"],
                "reject_frac": reject_frac,
                "reject_frac_negative": diag["reject_frac_negative"],
                "reject_frac_diverged": diag["reject_frac_diverged"],
            }
        }
        return mean, extras

    return MethodSpec(key, description, run)


METHOD_SPECS: tuple[MethodSpec, ...] = (
    _spec_1986(),
    _spec_1995("SPPR_1995_TE0.1", 0.1),
    _spec_1995("SPPR_1995_TEmean", "mean"),
    _spec_ulanowicz("Ulanowicz_globalTEmean", "global", "mean",
                    "Nullspace form of the EwE path sum, global TE at the true mean "
                    "(Jensen-able on TL)."),
    _spec_ulanowicz("Ulanowicz_TE", "TE", None,
                    "Nullspace form of the EwE path sum, per-group TE (Jensen-able on TE)."),
    MethodSpec("EwE_TE_noEE",
               "EwE path enumeration without EE weighting -- heavy method.",
               # silent=True for the same reason as the Monte-Carlo methods: silent=False routes
               # through tqdm.notebook, which raises ImportError('IProgress not found') outside
               # Jupyter and takes the whole method down.
               lambda m, mc_samples: (
                   _first(m.SPPR_EwE(TE_option="TE", use_EE=False, return_paths=False,
                                     silent=True, max_paths=100_000)), {})),
    MethodSpec("EwE_TE_EE",
                   "EwE path enumeration with EE weighting -- heavy method.",
                   # silent=True for the same reason as the Monte-Carlo methods: silent=False routes
                   # through tqdm.notebook, which raises ImportError('IProgress not found') outside
                   # Jupyter and takes the whole method down.
                   lambda m, mc_samples: (
                       _first(m.SPPR_EwE(TE_option="TE", use_EE=True, return_paths=True,
                                         silent=True, max_paths=100_000)), {})),
    MethodSpec("SPPR_2015",
               "2015 Leontief method, includes cycles. Takes no arguments: PP-only filtering is "
               "applied downstream via get_PPR.",
               lambda m, mc_samples: (_first(m.SPPR_2015()), {})),
    _spec_new("new_TE_noEEfix", "TE",
              "Primary solver, TE_option='TE', EE=0 groups left unrepaired (reproduces 2015).",
              fix_EE_0_cases=False),
    _spec_new("new_TE_EEfix", "TE",
              "Primary solver, TE_option='TE', EE=0 groups repaired.",
              fix_EE_0_cases=True),
    _spec_new("new_GE", "GE", "Primary solver, TE_option='GE' (should sit below TE)."),
    _spec_new("new_WithEgestion", "With Egestion",
              "Primary solver, TE_option='With Egestion' (should sit below GE)."),
    _spec_symbolic("TE", "as_PP"),
    _spec_symbolic("TE", "as_DC"),
    _spec_symbolic("GE", "as_PP"),
    _spec_symbolic("GE", "as_DC"),
    _spec_symbolic("With Egestion", "as_PP"),
    _spec_symbolic("With Egestion", "as_DC"),
    _spec_monte_carlo("MC_new_GE", "GE", None,
                      "Monte-Carlo over TE (10%/20%), kind='new', TE_option='GE', diverged draws "
                      "excluded."),
    _spec_monte_carlo("MC_new_TE_EEfix", "TE", {"fix_EE_0_cases": True},
                      "Monte-Carlo over TE (10%/20%), kind='new', TE_option='TE' with "
                      "fix_EE_0_cases=True, diverged draws excluded."),
)

SPEC_BY_KEY = {s.key: s for s in METHOD_SPECS}
ALL_METHOD_KEYS: tuple[str, ...] = tuple(s.key for s in METHOD_SPECS)

# SPPR_EwE enumerates every path in the food web; it is orders of magnitude slower than the rest.
HEAVY_METHOD_KEYS: tuple[str, ...] = ("EwE_TE_EE",)
FAST_METHOD_KEYS: tuple[str, ...] = tuple(k for k in ALL_METHOD_KEYS if k not in HEAVY_METHOD_KEYS)


# --------------------------------------------------------------------------- timed execution

# A unit of work the runner can time out. Both kinds take the model and return something
# picklable, so either can be handed to a worker process.
TASK_METHOD = "method"
TASK_HEALTH = "health"

# Every task comes back as (status, payload, warning_notes):
#   'ok'          payload is the task's return value
#   'error'       payload is (exception_name, message, traceback_text)
#   'timeout'     payload is None -- the worker was killed mid-task
#   'crashed'     payload is a string -- the worker died without answering
#   'unavailable' payload is a string -- no worker could be used at all (caller falls back)
# warning_notes is a list of (category_name, message) captured while the task ran.


def _dispatch_task(model: PPRCalculator, kind: str, payload, mc_samples: int):
    """Run one timed unit of work against `model`.

    Kept at module level and keyed by name rather than by callable: the worker process imports
    this module and looks the spec up itself, because MethodSpec.run is a closure and closures
    do not survive the pickling a spawned process needs.
    """
    if kind == TASK_METHOD:
        return SPEC_BY_KEY[payload].run(model, mc_samples)
    if kind == TASK_HEALTH:
        return model.diagnose_sppr(TE_option=payload, **HEALTH_DET_CONFIG)
    raise ValueError(f"unknown task kind {kind!r}")


def _run_task_inline(model: PPRCalculator, kind: str, payload, mc_samples: int) -> tuple:
    """Run one task in this process, with no time budget, in the standard envelope."""
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = _dispatch_task(model, kind, payload, mc_samples)
        return "ok", result, [(w.category.__name__, str(w.message)) for w in caught]
    except Exception as exc:
        return "error", (type(exc).__name__, str(exc), traceback.format_exc()), []


def _drop_unpicklable_extras(result, notes: list, exc: Exception) -> tuple:
    """Salvage a method result whose extras dict cannot cross the process boundary.

    Only SPPR_symbolic is really at risk: it hands back sympy equations and symbols that this
    exporter never writes to the workbook anyway. Dropping them keeps the SPPR itself, which is
    the part every sheet is built from.
    """
    note = ("PicklingError",
            f"extras dropped, they could not be sent back from the worker "
            f"({type(exc).__name__}: {exc})")
    if isinstance(result, tuple) and len(result) == 2:
        return "ok", (result[0], {}), list(notes) + [note]
    return "error", ("PicklingError", f"result could not be sent back: {exc}", ""), list(notes)


def _worker_main(conn) -> None:
    """Worker-process entry point: hold one model and answer one task at a time over `conn`.

    The loop is deliberately dumb -- it decides nothing, it only runs what it is told and reports
    back -- because the parent kills it without warning whenever a task overruns its budget.
    """
    model = None
    try:
        while True:
            try:
                message = conn.recv()
            except (EOFError, KeyboardInterrupt):
                return

            command = message[0]
            if command == "stop":
                return
            if command == "model":
                model = message[1]
                conn.send(("ok", None, []))
                continue
            if command != "task":
                conn.send(("error", ("ValueError", f"unknown command {command!r}", ""), []))
                continue

            _, kind, payload, mc_samples = message
            status, result, notes = _run_task_inline(model, kind, payload, mc_samples)
            try:
                conn.send((status, result, notes))
            except Exception as exc:  # the result itself did not survive pickling
                conn.send(_drop_unpicklable_extras(result, notes, exc))
    except Exception:
        # Never die silently: the parent reads an unexplained exit as a crash, which is correct,
        # but the traceback belongs on stderr where the run can be debugged from.
        traceback.print_exc()
    finally:
        try:
            conn.close()
        except Exception:
            pass


class _TaskRunner:
    """Runs SPPR tasks in a worker process, killing the worker when one overruns its budget.

    One worker serves as many models and methods as it can: the model is sent once and every
    task afterwards runs against it, so a directory run pays for a single process spawn rather
    than one per method. When a task overruns `timeout` the worker is killed -- the only way to
    stop a CPU-bound call -- and the next task transparently starts a fresh worker and re-sends
    the model to it.

    If a worker cannot be used at all (the model will not pickle, the process will not spawn),
    the runner disables itself and every later task reports 'unavailable', so the caller can fall
    back to running in-process. That costs the time budget, never the results.

    Attributes:
        timeout: wall-clock seconds a single task may run for.
        disabled_reason: why the worker was abandoned, or None while it is usable.
    """

    def __init__(self, timeout: float = DEFAULT_METHOD_TIMEOUT) -> None:
        self.timeout = float(timeout)
        self.disabled_reason: Optional[str] = None
        # 'spawn' explicitly: it is the only start method on Windows, and a forked worker would
        # inherit the parent's numpy/BLAS thread state, which is not safe to run in.
        self._ctx = mp.get_context("spawn")
        self._proc = None
        self._conn = None
        self._model = None
        self._model_sent = False

    # -- lifecycle ---------------------------------------------------------

    def set_model(self, model: PPRCalculator) -> None:
        """Point the runner at a model. Re-setting the same object is a no-op.

        A new model also clears `disabled_reason`: whatever went wrong was most likely a
        property of the previous model (an attribute that would not pickle), and one awkward
        model must not cost every later model in the directory its time budget.
        """
        if model is self._model:
            return
        self._model = model
        self._model_sent = False
        self.disabled_reason = None

    def close(self) -> None:
        """Shut the worker down, politely first and then by force."""
        if self._proc is not None and self._conn is not None and self._proc.is_alive():
            try:
                self._conn.send(("stop",))
                self._proc.join(5)
            except Exception:
                pass
        self._kill()
        self._model = None

    def __enter__(self) -> "_TaskRunner":
        return self

    def __exit__(self, *exc_info) -> None:
        self.close()

    # -- the one call that matters -----------------------------------------

    def run(self, kind: str, payload, mc_samples: int) -> tuple:
        """Run one task under the time budget and return the standard envelope."""
        if self.disabled_reason is not None:
            return "unavailable", self.disabled_reason, []
        try:
            self._ensure_ready()
        except Exception as exc:
            self._disable(f"{type(exc).__name__}: {exc}")
            return "unavailable", self.disabled_reason, []

        try:
            self._conn.send(("task", kind, payload, mc_samples))
            answered = self._conn.poll(self.timeout)
        except Exception as exc:
            self._kill()
            return "crashed", f"could not reach the worker ({type(exc).__name__}: {exc})", []

        if not answered:
            self._kill()
            return "timeout", None, []

        try:
            status, result, notes = self._conn.recv()
        except Exception as exc:
            code = self._proc.exitcode if self._proc is not None else None
            self._kill()
            return "crashed", (f"the worker exited (exitcode {code}) while running the task "
                               f"[{type(exc).__name__}]"), []
        return status, result, notes

    # -- internals ---------------------------------------------------------

    def _ensure_ready(self) -> None:
        if self._proc is None or not self._proc.is_alive():
            self._start()
        if not self._model_sent:
            self._send_model()

    def _start(self) -> None:
        self._kill()
        parent_conn, child_conn = self._ctx.Pipe(duplex=True)
        proc = self._ctx.Process(target=_worker_main, args=(child_conn,), daemon=True)
        proc.start()
        # The parent must drop its own copy of the child end, or the pipe never reports EOF when
        # the worker dies and a recv on a dead worker would block instead of raising.
        child_conn.close()
        self._proc, self._conn = proc, parent_conn
        self._model_sent = False

    def _send_model(self) -> None:
        if self._model is None:
            raise RuntimeError("no model has been set on the task runner")
        self._conn.send(("model", self._model))
        if not self._conn.poll(MODEL_HANDSHAKE_TIMEOUT):
            self._kill()
            raise TimeoutError(
                f"the worker did not take the model within {MODEL_HANDSHAKE_TIMEOUT:.0f}s")
        self._conn.recv()
        self._model_sent = True

    def _kill(self) -> None:
        proc, conn = self._proc, self._conn
        self._proc, self._conn = None, None
        self._model_sent = False
        if conn is not None:
            try:
                conn.close()
            except Exception:
                pass
        if proc is None:
            return
        try:
            if proc.is_alive():
                proc.terminate()
                proc.join(5)
            if proc.is_alive():
                proc.kill()
                proc.join(5)
        except Exception:
            pass
        try:
            proc.close()
        except Exception:
            pass

    def _disable(self, reason: str) -> None:
        self._kill()
        self.disabled_reason = reason


def _run_task(model: PPRCalculator, kind: str, payload, mc_samples: int,
              runner: Optional[_TaskRunner]) -> tuple:
    """Run one task under `runner`'s budget, falling back to this process if there is none."""
    if runner is None:
        return _run_task_inline(model, kind, payload, mc_samples)
    status, result, notes = runner.run(kind, payload, mc_samples)
    if status == "unavailable":
        return _run_task_inline(model, kind, payload, mc_samples)
    return status, result, notes


def _budget_text(runner: Optional[_TaskRunner]) -> str:
    return f"{runner.timeout:.0f}s" if runner is not None else "no time budget"


# --------------------------------------------------------------------------- result container

@dataclass
class ModelTables:
    """Every table for one model, plus everything worth telling the user about the run."""
    model_label: str
    source_file: str
    groups: pd.DataFrame
    sppr_pp: pd.DataFrame
    sppr_inner: pd.DataFrame
    sppr_all: pd.DataFrame
    health: pd.DataFrame
    footprint: pd.DataFrame
    mc_diagnostics: pd.DataFrame
    notes: pd.DataFrame
    issues: list = field(default_factory=list)
    # method key -> 'ok' | 'failed' | 'timeout' | 'crashed'. The SPPR sheets cannot tell those
    # apart -- every one of them is a column of NaN -- so the outcome is carried separately and
    # written into the run_notes sheet.
    method_status: dict = field(default_factory=dict)

    @property
    def sheets(self) -> dict:
        return {
            SHEET_GROUPS: self.groups,
            SHEET_SPPR_PP: self.sppr_pp,
            SHEET_SPPR_INNER: self.sppr_inner,
            SHEET_SPPR_ALL: self.sppr_all,
            SHEET_HEALTH: self.health,
            SHEET_FOOTPRINT: self.footprint,
            SHEET_MC: self.mc_diagnostics,
            SHEET_NOTES: self.notes,
        }


def _issue(kind: str, detail: str, method: Optional[str] = None) -> dict:
    return {"kind": kind, "method": method, "detail": detail}


# --------------------------------------------------------------------------- helpers

def _group_type_map(model: PPRCalculator) -> dict:
    types = {}
    for seq in model.get_PP_seq():
        types[int(seq)] = "PP"
    for seq in model.get_DET_seq():
        types[int(seq)] = "DET"
    for seq in model.get_Import_seq():
        types[int(seq)] = "Import"
    for seq in model.get_Regular_seq():
        types[int(seq)] = "Regular"
    return types


def _to_seq_frame(model: PPRCalculator, sppr) -> pd.DataFrame:
    """Normalise any SPPR result onto seq-indexed rows in groups_df order."""
    if isinstance(sppr, pd.Series):
        sppr = sppr.to_frame()
    sppr = PPRCalculator.rename_results(sppr.copy(), model.name2seq)
    return sppr.reindex(model.get_groups_df().index)


def _source_sums(model: PPRCalculator, raw: pd.DataFrame, source_resolved: bool) -> dict:
    """SUM_ALL / SUM_INNER / SUM_PP from a method's own columns, per get_PPR semantics."""
    det = [c for c in model.get_DET_seq() if c in raw.columns]
    imp = [c for c in model.get_Import_seq() if c in raw.columns]

    sums = {
        "SUM_ALL": raw.sum(axis=1, min_count=1),
        "SUM_INNER": raw.drop(columns=imp).sum(axis=1, min_count=1),
    }
    if source_resolved:
        sums["SUM_PP"] = raw.drop(columns=det + imp).sum(axis=1, min_count=1)
    else:
        # An un-attributed total cannot be split into a PP-only part.
        sums["SUM_PP"] = pd.Series(np.nan, index=raw.index)
    return sums


def _footprint_row(model: PPRCalculator, raw: pd.DataFrame, source_resolved: bool) -> dict:
    """The same six quantities diagnose_sppr reports under 'footprint', for any method."""
    npp = float(model.get_NPP(only_inner=True))
    ppr_all = float(model.get_PPR(raw, only_inner=False).sum(axis=1).sum())
    ppr_inner = float(model.get_PPR(raw, only_inner=True).sum(axis=1).sum())
    ppr_pp_only = (float(model.get_PPR(raw, only_pp=True).sum(axis=1).sum())
                   if source_resolved else np.nan)
    return {
        "ppr_all": ppr_all,
        "ppr_inner": ppr_inner,
        "ppr_pp_only": ppr_pp_only,
        "npp": npp,
        "ppr2npp": (ppr_inner / npp) if npp else np.nan,
        "ppr2npp_pp_only": (ppr_pp_only / npp) if (npp and source_resolved) else np.nan,
    }


# --------------------------------------------------------------------------- table builders

def build_groups_table(model: PPRCalculator) -> pd.DataFrame:
    """The per-group parameter table, with the group-type flag the calculators key off."""
    groups = model.get_groups_df()
    types = _group_type_map(model)
    groups.insert(1, "group_type", [types.get(int(s), "Regular") for s in groups.index])
    groups.index.name = "seq"
    return groups


def build_sppr_tables(model: PPRCalculator, results: dict,
                      method_keys: Sequence[str]) -> dict:
    """The three flat SPPR tables, one per source scope.

    Each is a groups x methods table: the index is group seq, the first column is the group name,
    and there is one column per method holding that method's SPPR summed over the sources in
    scope ('pp' drops Import and Detritus, 'inner' drops Import, 'all' keeps everything). No
    MultiIndex is used on either axis. A method that failed, or whose result cannot be attributed
    to the scope, is NaN.

    Returns:
        dict: {'pp': DataFrame, 'inner': DataFrame, 'all': DataFrame}.
    """
    groups = model.get_groups_df()
    index = pd.Index(list(groups.index), name="seq")
    names = [model.seq2name.get(int(s), str(s)) for s in groups.index]

    tables = {}
    for scope in SPPR_SHEET_SUMS:
        table = pd.DataFrame(index=index)
        table["group_name"] = names
        tables[scope] = table

    for key in method_keys:
        raw = results.get(key)
        if raw is None:
            for table in tables.values():
                table[key] = np.nan
            continue
        sums = _source_sums(model, raw, SPEC_BY_KEY[key].source_resolved)
        for scope, sum_key in SPPR_SHEET_SUMS.items():
            tables[scope][key] = sums[sum_key].to_numpy(dtype=float)

    return tables


def build_footprint_table(model: PPRCalculator, results: dict,
                          method_keys: Sequence[str]) -> pd.DataFrame:
    rows = {}
    for key in method_keys:
        raw = results.get(key)
        if raw is None:
            rows[key] = {c: np.nan for c in FOOTPRINT_COLUMNS}
        else:
            rows[key] = _footprint_row(model, raw, SPEC_BY_KEY[key].source_resolved)
    table = pd.DataFrame.from_dict(rows, orient="index")
    table = table.reindex(index=list(method_keys), columns=list(FOOTPRINT_COLUMNS))
    table.index.name = "method"
    return table


def build_health_table(model: PPRCalculator, issues: list,
                       runner: Optional[_TaskRunner] = None) -> pd.DataFrame:
    """diagnose_sppr under the fixed detritus config, one row per TE_option.

    Each row runs under the same time budget as an SPPR method -- diagnose_sppr solves the whole
    system once, so on a pathological model it can hang for exactly the same reasons. A row that
    times out is marked TIMEOUT rather than ERROR, because nothing was actually diagnosed.
    """
    rows = {}
    for te in HEALTH_TE_OPTIONS:
        status, payload, notes = _run_task(model, TASK_HEALTH, te, 0, runner)
        for category, message in notes:
            issues.append(_issue("warning", f"[diagnose_sppr {te}] {category}: {message}"))

        if status == "ok":
            report = payload
            for text in report.get("warnings", []):
                issues.append(_issue("health_warning", f"[{te}] {text}"))
            rows[te] = PPRCalculator._flatten_diagnostics(report)
        elif status == "timeout":
            detail = f"still running after {_budget_text(runner)}; the worker was killed"
            issues.append(_issue("health_timeout", f"[{te}] {detail}"))
            rows[te] = {"status": "TIMEOUT", "error": detail}
        elif status == "crashed":
            issues.append(_issue("health_failed", f"[{te}] worker process died: {payload}"))
            rows[te] = {"status": "ERROR", "error": f"worker process died: {payload}"}
        else:
            name, message, _ = payload
            issues.append(_issue("health_failed", f"[{te}] {name}: {message}"))
            rows[te] = {"status": "ERROR", "error": f"{name}: {message}"}

    table = pd.DataFrame.from_dict(rows, orient="index").reindex(list(HEALTH_TE_OPTIONS))
    table.index.name = "TE_option"

    # status first, then the echoed configuration, then everything else -- readable left to right.
    lead = [c for c in ("status",) if c in table.columns]
    cfg = sorted(c for c in table.columns if c.startswith("config_"))
    rest = [c for c in table.columns if c not in lead + cfg]
    return table[lead + cfg + rest]


def build_mc_table(results_extras: dict, method_keys: Sequence[str]) -> pd.DataFrame:
    columns = ("n_samples", "n_accepted", "n_rejected_negative", "n_rejected_diverged",
               "reject_frac", "reject_frac_negative", "reject_frac_diverged")
    rows = {k: results_extras[k]["mc"] for k in method_keys
            if k in results_extras and "mc" in results_extras[k]}
    table = pd.DataFrame.from_dict(rows, orient="index")
    if table.empty:
        table = pd.DataFrame(columns=list(columns))
    table = table.reindex(columns=list(columns))
    table.index.name = "method"
    return table


def build_notes_table(method_keys: Sequence[str],
                      method_status: Optional[dict] = None) -> pd.DataFrame:
    """The conventions sheet, plus one row per method with its description and run outcome.

    The status column is the only place inside the workbook that separates a method which timed
    out from one that raised or was never run: all three leave nothing but NaN in the SPPR
    sheets. The run report carries the same information with the full reason attached.
    """
    notes = [{"topic": topic, "note": text, "status": ""} for topic, text in RUN_NOTES]
    notes += [{"topic": f"method: {k}", "note": SPEC_BY_KEY[k].description,
               "status": "" if method_status is None else method_status.get(k, "not run")}
              for k in method_keys]
    return pd.DataFrame(notes).set_index("topic")


# --------------------------------------------------------------------------- orchestration

def build_model_tables(model: PPRCalculator, *, method_keys: Optional[Sequence[str]] = None,
                       mc_samples: int = DEFAULT_MC_SAMPLES,
                       model_label: str = "", source_file: str = "",
                       method_timeout: Optional[float] = DEFAULT_METHOD_TIMEOUT,
                       runner: Optional[_TaskRunner] = None) -> ModelTables:
    """Run the selected SPPR methods on one model and assemble every table.

    A method that raises is recorded in `issues` and leaves its whole block NaN; it never aborts
    the run and its cells are never zero-filled. A method that overruns `method_timeout` is
    treated exactly the same way, except that its worker is killed and the outcome is recorded as
    a timeout rather than a failure -- the two mean very different things when the numbers are
    read later, and neither is visible in the sheets themselves.

    Args:
        model: the calculator to run.
        method_keys: which methods to run, defaulting to all of them (including the heavy
            SPPR_EwE). Pass FAST_METHOD_KEYS to skip the heavy one.
        mc_samples: draws per Monte-Carlo method.
        model_label: human-readable model name for the report.
        source_file: the JSON this model came from.
        method_timeout: wall-clock seconds one method (or one diagnose_sppr row) may run for
            before its worker is killed. None or 0 runs everything in this process, unbudgeted.
        runner: an existing _TaskRunner to reuse, so that one worker process can serve many
            models. When omitted, one is created for this call and closed at the end of it;
            `method_timeout` is then what sets its budget.

    Returns:
        ModelTables: the eight sheets, the issue list and the per-method outcome map.
    """
    keys = list(ALL_METHOD_KEYS if method_keys is None else method_keys)
    unknown = [k for k in keys if k not in SPEC_BY_KEY]
    if unknown:
        raise KeyError(f"unknown method key(s) {unknown}; known keys: {list(ALL_METHOD_KEYS)}")

    issues: list = []
    results: dict = {}
    extras: dict = {}
    method_status: dict = {}

    owns_runner = False
    if runner is None and method_timeout:
        runner = _TaskRunner(method_timeout)
        owns_runner = True
    if runner is not None:
        runner.set_model(model)

    try:
        for key in keys:
            status, payload, notes = _run_task(model, TASK_METHOD, key, mc_samples, runner)
            for category, message in notes:
                issues.append(_issue("warning", f"{category}: {message}", key))

            if status == "timeout":
                method_status[key] = "timeout"
                issues.append(_issue(
                    "method_timeout",
                    f"still running after {_budget_text(runner)}; the worker was killed and the "
                    f"method left empty", key))
                continue
            if status == "crashed":
                method_status[key] = "crashed"
                issues.append(_issue("method_failed", f"worker process died: {payload}", key))
                continue
            if status == "error":
                name, message, tb = payload
                method_status[key] = "failed"
                issues.append(_issue("method_failed", f"{name}: {message}", key))
                issues.append(_issue("traceback", tb, key))
                continue

            try:
                sppr, extra = payload
                frame = _to_seq_frame(model, sppr)
                values = frame.to_numpy(dtype=float)
            except Exception as exc:
                method_status[key] = "failed"
                issues.append(_issue(
                    "method_failed",
                    f"result could not be read back: {type(exc).__name__}: {exc}", key))
                issues.append(_issue("traceback", traceback.format_exc(), key))
                continue

            results[key] = frame
            extras[key] = extra
            method_status[key] = "ok"

            if np.any(values < -1e-10):
                issues.append(_issue(
                    "negative_sppr",
                    f"{int(np.sum(values < -1e-10))} negative SPPR cell(s), min "
                    f"{np.nanmin(values):.6g}", key))
            if not np.isfinite(values[~np.isnan(values)]).all():
                issues.append(_issue("non_finite_sppr", "contains infinite values", key))

        health = build_health_table(model, issues, runner)
    finally:
        if owns_runner and runner is not None:
            runner.close()

    if runner is not None and runner.disabled_reason:
        issues.append(_issue(
            "timeout_unavailable",
            f"the time budget could not be enforced, every method ran in-process: "
            f"{runner.disabled_reason}"))

    try:
        balanced, _, _ = model.is_model_balanced()
        if not balanced:
            issues.append(_issue("mass_balance", "is_model_balanced() is False for this model"))
    except Exception as exc:
        issues.append(_issue("mass_balance", f"is_model_balanced() raised {type(exc).__name__}: {exc}"))

    for key in keys:
        if key not in results:
            continue
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ok, inflow, outflow = model.is_sppr_balanced(results[key])
            if not ok:
                issues.append(_issue(
                    "sppr_balance", f"inflow {inflow:.6g} != outflow {outflow:.6g}", key))
        except Exception as exc:
            issues.append(_issue("sppr_balance",
                                 f"is_sppr_balanced raised {type(exc).__name__}: {exc}", key))

    sppr_tables = build_sppr_tables(model, results, keys)

    return ModelTables(
        model_label=model_label,
        source_file=source_file,
        groups=build_groups_table(model),
        sppr_pp=sppr_tables["pp"],
        sppr_inner=sppr_tables["inner"],
        sppr_all=sppr_tables["all"],
        health=health,
        footprint=build_footprint_table(model, results, keys),
        mc_diagnostics=build_mc_table(extras, keys),
        notes=build_notes_table(keys, method_status),
        issues=issues,
        method_status=method_status,
    )


def load_model(model_path: str) -> tuple[PPRCalculator, str]:
    """Load one Ecopath JSON into a calculator, returning it with a human-readable label."""
    model_data = ModelData(model_path)
    model = PPRCalculator.from_modeldata(
        model_data, underdetermined=True, zero_biomass_accum=False)
    label = f"{model_data.model_name} ({model_data.model_year})"
    return model, label


# --------------------------------------------------------------------------- excel writing

_COLUMN_WIDTH = 10
_INDEX_WIDTH = 34


def _autoformat(writer, sheet_name: str, table: pd.DataFrame, n_index_cols: int,
                n_header_rows: int) -> None:
    """Freeze the labels, widen the columns and give the numbers a readable format."""
    try:
        worksheet = writer.sheets[sheet_name]
    except (AttributeError, KeyError):  # pragma: no cover - engine without sheet access
        return
    try:
        from openpyxl.utils import get_column_letter
    except ImportError:  # pragma: no cover
        return

    worksheet.freeze_panes = worksheet.cell(row=n_header_rows + 1, column=n_index_cols + 1)
    for i in range(n_index_cols):
        worksheet.column_dimensions[get_column_letter(1 + i)].width = _INDEX_WIDTH
    for i in range(table.shape[1]):
        worksheet.column_dimensions[get_column_letter(1 + n_index_cols + i)].width = _COLUMN_WIDTH

    for row in worksheet.iter_rows(min_row=n_header_rows + 1,
                                   min_col=n_index_cols + 1):
        for cell in row:
            if isinstance(cell.value, float):
                cell.number_format = "0.000000"


def write_model_excel(model_path: str, out_dir: str = DEFAULT_OUT_DIR, *,
                      method_keys: Optional[Sequence[str]] = None,
                      mc_samples: int = DEFAULT_MC_SAMPLES,
                      method_timeout: Optional[float] = DEFAULT_METHOD_TIMEOUT) -> str:
    """Build every table for one model JSON and write it to a single workbook.

    Args:
        model_path: the Ecopath model JSON to run.
        out_dir: where the workbook is written.
        method_keys: which methods to run; defaults to all of them.
        mc_samples: draws per Monte-Carlo method.
        method_timeout: wall-clock seconds per method; None or 0 disables the budget.

    Returns:
        str: the path of the workbook written.
    """
    model, label = load_model(model_path)
    tables = build_model_tables(model, method_keys=method_keys, mc_samples=mc_samples,
                                model_label=label, source_file=os.path.basename(model_path),
                                method_timeout=method_timeout)
    return write_tables_excel(tables, model_path, out_dir)


def write_tables_excel(tables: ModelTables, model_path: str,
                       out_dir: str = DEFAULT_OUT_DIR) -> str:
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    out_path = os.path.join(out_dir, Path(model_path).stem + ".xlsx")

    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        for sheet, table in tables.sheets.items():
            table.to_excel(writer, sheet_name=sheet)
            # Every sheet is flat: one index column, one header row.
            _autoformat(writer, sheet, table, 1, 1)
    return out_path


SHEET_READ_LAYOUT = {
    SHEET_GROUPS: dict(index_col=0, header=0),
    SHEET_SPPR_PP: dict(index_col=0, header=0),
    SHEET_SPPR_INNER: dict(index_col=0, header=0),
    SHEET_SPPR_ALL: dict(index_col=0, header=0),
    SHEET_HEALTH: dict(index_col=0, header=0),
    SHEET_FOOTPRINT: dict(index_col=0, header=0),
    SHEET_MC: dict(index_col=0, header=0),
    SHEET_NOTES: dict(index_col=0, header=0),
}


def read_pprs_excel(path: str) -> dict:
    """Read a workbook written by `write_model_excel` back into its tables.

    Returns:
        dict: sheet name -> DataFrame, with the same index and column structure that was written
        (the sppr_table keeps its (method, basal_source) MultiIndex columns and (seq, group_name)
        MultiIndex index).
    """
    available = set(pd.ExcelFile(path).sheet_names)
    out = {}
    for sheet, kwargs in SHEET_READ_LAYOUT.items():
        if sheet not in available:
            continue
        out[sheet] = pd.read_excel(path, sheet_name=sheet, **kwargs)
    return out


# --------------------------------------------------------------------------- run report

_ISSUE_HEADINGS = [
    ("method_timeout", "METHODS THAT TIMED OUT (no result -- their cells are empty, not zero)"),
    ("method_failed", "METHODS THAT FAILED (their cells are empty, not zero)"),
    ("health_timeout", "diagnose_sppr TIMEOUTS"),
    ("health_failed", "diagnose_sppr FAILURES"),
    ("timeout_unavailable", "TIME BUDGET COULD NOT BE ENFORCED"),
    ("negative_sppr", "NEGATIVE SPPR VALUES"),
    ("non_finite_sppr", "NON-FINITE SPPR VALUES"),
    ("mass_balance", "MODEL MASS-BALANCE PROBLEMS"),
    ("sppr_balance", "SPPR BALANCE PROBLEMS (inflow != outflow)"),
    ("health_warning", "diagnose_sppr WARNINGS"),
    ("warning", "PYTHON WARNINGS RAISED DURING THE RUN"),
]


def _format_model_section(entry: dict) -> list:
    lines = [f"MODEL: {entry['source_file']}"]
    if entry.get("model_label"):
        lines.append(f"  name           : {entry['model_label']}")
    if entry.get("output"):
        lines.append(f"  workbook       : {entry['output']}")
    if entry.get("load_error"):
        lines.append(f"  !! COULD NOT LOAD: {entry['load_error']}")
        return lines + [""]
    if entry.get("skipped"):
        lines.append("  skipped        : workbook already present (resume)")
        return lines + [""]

    lines.append(f"  health         : {entry.get('health_summary', 'n/a')}")
    ok_line = f"  methods ok     : {entry.get('n_ok', 0)}/{entry.get('n_methods', 0)}"
    if entry.get("n_timeout"):
        ok_line += f"   ({entry['n_timeout']} timed out)"
    lines.append(ok_line)
    if entry.get("mc_summary"):
        lines.append(f"  monte-carlo    : {entry['mc_summary']}")

    issues = entry.get("issues", [])
    for kind, heading in _ISSUE_HEADINGS:
        matching = [i for i in issues if i["kind"] == kind]
        if not matching:
            continue
        lines.append(f"  {heading}:")
        seen = set()
        for item in matching:
            text = f"{item['method'] + ': ' if item['method'] else ''}{item['detail']}"
            text = " ".join(str(text).split())
            if text in seen:
                continue
            seen.add(text)
            lines.append(f"    - {text}")
    lines.append("")
    return lines


def _health_summary(health: pd.DataFrame) -> str:
    if "status" not in health.columns:
        return "n/a"
    return ", ".join(f"{te}={health.loc[te, 'status']}" for te in health.index)


def _mc_summary(mc: pd.DataFrame) -> str:
    if mc.empty:
        return ""
    parts = []
    for key, row in mc.iterrows():
        parts.append(
            f"{key}: {int(row['n_accepted'])}/{int(row['n_samples'])} accepted "
            f"({int(row['n_rejected_negative'])} negative, "
            f"{int(row['n_rejected_diverged'])} diverged)")
    return "; ".join(parts)


def _timeout_label(method_timeout: Optional[float]) -> str:
    return f"{method_timeout:.0f}s per method" if method_timeout else "none (methods run in-process)"


def write_run_report(entries: list, out_dir: str, method_keys: Sequence[str],
                     mc_samples: int,
                     method_timeout: Optional[float] = DEFAULT_METHOD_TIMEOUT) -> str:
    """Write the human-readable run report and return its path."""
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    path = os.path.join(out_dir, LOG_FILENAME)

    n_skipped = sum(1 for e in entries if e.get("skipped"))
    n_written = sum(1 for e in entries if e.get("output") and not e.get("skipped"))
    n_failed = len(entries) - n_written - n_skipped

    header = [
        "=" * 88,
        "SPPR EXPORT RUN REPORT",
        "=" * 88,
        f"generated      : {datetime.now().isoformat(timespec='seconds')}",
        f"models found   : {len(entries)}",
        f"workbooks ok   : {n_written}",
        f"models skipped : {n_skipped} (resume: workbook already present)",
        f"models failed  : {n_failed}",
        f"methods run    : {len(method_keys)} -> {', '.join(method_keys)}",
        f"MC draws       : {mc_samples}",
        f"time budget    : {_timeout_label(method_timeout)}",
        "",
        "READ THIS BEFORE USING THE NUMBERS",
        "-" * 88,
    ]
    for topic, text in RUN_NOTES:
        header.append(f"* {topic}: {' '.join(text.split())}")
    header += ["", "=" * 88, "PER-MODEL DETAIL", "=" * 88, ""]

    body = []
    for entry in entries:
        body += _format_model_section(entry)

    text = "\n".join(header + body)
    Path(path).write_text(text, encoding="utf-8")
    return path


# --------------------------------------------------------------------------- directory run

def run_directory(json_dir: str = DEFAULT_JSON_DIR, out_dir: str = DEFAULT_OUT_DIR, *,
                  method_keys: Optional[Sequence[str]] = None,
                  mc_samples: int = DEFAULT_MC_SAMPLES,
                  resume: bool = False,
                  method_timeout: Optional[float] = DEFAULT_METHOD_TIMEOUT,
                  silent: bool = False) -> dict:
    """Write one workbook per model JSON in `json_dir`, plus a run report.

    A model that cannot be loaded, a method that raises, and a method that runs past its time
    budget are all recorded and the run continues.

    One worker process is shared by the whole directory: it is spawned on the first method of
    the first model and only replaced when a method has to be killed for overrunning, so the
    time budget costs a single process spawn rather than one per model.

    Args:
        json_dir: directory of Ecopath model JSONs.
        out_dir: where the workbooks and the report are written.
        method_keys: which methods to run; defaults to all of them.
        mc_samples: draws per Monte-Carlo method.
        resume: continue a previous run -- skip any model whose workbook is already in
            `out_dir`. The report still lists them, marked as skipped, but their per-model
            detail (health, issues, Monte-Carlo counts) is not re-derived.
        method_timeout: wall-clock seconds one method (or one diagnose_sppr row) may run for
            before it is killed and recorded as a timeout. None or 0 runs every method in this
            process with no budget, which is the old behaviour.
        silent: suppress the progress bar.

    Returns:
        dict: {'n_models', 'n_written', 'n_skipped', 'n_failed', 'n_timeout', 'report',
        'entries'}.
    """
    keys = list(ALL_METHOD_KEYS if method_keys is None else method_keys)
    paths = sorted(str(p) for p in Path(json_dir).glob("*.json"))
    entries = []

    runner = _TaskRunner(method_timeout) if method_timeout else None
    try:
        for model_path in tqdm(paths, disable=silent, desc="models"):
            entry = {"source_file": os.path.basename(model_path)}
            existing = os.path.join(out_dir, Path(model_path).stem + ".xlsx")
            if resume and os.path.exists(existing):
                entry["skipped"] = True
                entry["output"] = existing
                entries.append(entry)
                continue
            try:
                model, label = load_model(model_path)
                entry["model_label"] = label
                tables = build_model_tables(
                    model, method_keys=keys, mc_samples=mc_samples,
                    model_label=label, source_file=entry["source_file"],
                    method_timeout=method_timeout, runner=runner)
                entry["output"] = write_tables_excel(tables, model_path, out_dir)
                entry["issues"] = tables.issues
                entry["n_methods"] = len(keys)
                entry["n_ok"] = sum(1 for k in keys if tables.method_status.get(k) == "ok")
                entry["n_timeout"] = sum(
                    1 for k in keys if tables.method_status.get(k) == "timeout")
                entry["health_summary"] = _health_summary(tables.health)
                entry["mc_summary"] = _mc_summary(tables.mc_diagnostics)
            except Exception as exc:
                entry["load_error"] = f"{type(exc).__name__}: {exc}"
            entries.append(entry)
    finally:
        if runner is not None:
            runner.close()

    report = write_run_report(entries, out_dir, keys, mc_samples, method_timeout)
    n_skipped = sum(1 for e in entries if e.get("skipped"))
    n_written = sum(1 for e in entries if e.get("output") and not e.get("skipped"))
    return {
        "n_models": len(entries),
        "n_written": n_written,
        "n_skipped": n_skipped,
        "n_failed": len(entries) - n_written - n_skipped,
        "n_timeout": sum(e.get("n_timeout", 0) for e in entries),
        "report": report,
        "entries": entries,
    }


def _pop_timeout_flag(argv: list) -> Optional[float]:
    """Pull --timeout SECONDS / --timeout=SECONDS out of argv, in place.

    'none', 'off' and 0 all mean the same thing: run every method in-process with no budget.
    """
    value = None
    for i, arg in enumerate(list(argv)):
        if arg in ("--timeout", "-t"):
            if i + 1 >= len(argv):
                raise SystemExit("--timeout needs a value in seconds ('none' or 0 to disable it)")
            value = argv[i + 1]
            del argv[i:i + 2]
            break
        if arg.startswith("--timeout="):
            value = arg.split("=", 1)[1]
            del argv[i]
            break

    if value is None:
        return DEFAULT_METHOD_TIMEOUT
    if str(value).strip().lower() in ("none", "off"):
        return None
    try:
        seconds = float(value)
    except ValueError:
        raise SystemExit(f"--timeout expects seconds, got {value!r}")
    return seconds if seconds > 0 else None


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    resume = False
    for flag in ("--resume", "-r"):
        if flag in argv:
            argv.remove(flag)
            resume = True
    method_timeout = _pop_timeout_flag(argv)
    json_dir = argv[0] if argv else DEFAULT_JSON_DIR
    out_dir = argv[1] if len(argv) > 1 else DEFAULT_OUT_DIR

    summary = run_directory(json_dir, out_dir, resume=resume, method_timeout=method_timeout)
    print(f"wrote {summary['n_written']}/{summary['n_models']} workbooks to {out_dir}")
    if summary["n_skipped"]:
        print(f"skipped {summary['n_skipped']} model(s) already present (--resume)")
    print(f"report: {summary['report']}")
    if summary["n_timeout"]:
        print(f"{summary['n_timeout']} method run(s) hit the "
              f"{_timeout_label(method_timeout)} budget and were left empty -- see the report")
    if summary["n_failed"]:
        print(f"{summary['n_failed']} model(s) could not be processed -- see the report")
    return 0


if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    json_dir = os.path.join("real_models", "EwE_jsons")
    out_dir = "output/Ecobase_models"

    run_directory(json_dir=json_dir, out_dir=out_dir, resume=True, silent=False,
                  method_timeout=DEFAULT_METHOD_TIMEOUT)

    # raise SystemExit(main())
