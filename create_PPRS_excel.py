"""Export one Excel workbook per Ecopath model, holding every SPPR method side by side.

For each model JSON in a directory this builds a workbook with:

    groups_df       -- the per-group parameters the SPPR / PPR / NPP calculations consume
    sppr_table      -- SPPR per group x basal source x method, plus per-method source sums
    model_health    -- diagnose_sppr output, one row per TE_option, under a fixed detritus config
    footprint       -- PPR / NPP footprint per method
    mc_diagnostics  -- Monte-Carlo accept/reject breakdown per MC method
    run_notes       -- the conventions and caveats needed to read the numbers correctly

and a plain-text run report next to them listing every warning, skip and failure.

`read_pprs_excel` reads a workbook back into the same tables, index and column structure.

Conventions used throughout, and repeated in the run_notes sheet:
  * NaN means "not available", never zero. A source a method does not resolve is NaN; a method
    that raised leaves its whole block NaN.
  * SUM_ALL / SUM_INNER / SUM_PP follow get_PPR's own source semantics: all sources, all minus
    Import, and all minus Import and Detritus respectively.
  * Methods that return a single un-attributed SPPR column (SPPR_1986, SPPR_1995) are stored
    under the single basal source 'AGGREGATE' and get NaN for SUM_PP and ppr_pp_only, because
    their value cannot honestly be attributed to primary production alone.

Requires `openpyxl` (the .xlsx engine pandas uses here); everything else is already a project
dependency.

Entry point:
    python create_PPRS_excel.py [json_dir] [out_dir]

Defaults to real_models/EwE_jsons -> output/.
"""
from __future__ import annotations

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

AGGREGATE_SOURCE = "AGGREGATE"
SUM_COLUMNS = ("SUM_PP", "SUM_INNER", "SUM_ALL")

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
SHEET_SPPR = "sppr_table"
SHEET_HEALTH = "model_health"
SHEET_FOOTPRINT = "footprint"
SHEET_MC = "mc_diagnostics"
SHEET_NOTES = "run_notes"

LOG_FILENAME = "sppr_export_report.txt"
DEFAULT_JSON_DIR = os.path.join("real_models", "EwE_jsons")
DEFAULT_OUT_DIR = "output"
DEFAULT_MC_SAMPLES = 100

RUN_NOTES = [
    ("NaN convention",
     "NaN means not available, never zero. A basal source a method does not resolve is NaN, and "
     "a method that raised leaves its entire block NaN."),
    ("SUM_ALL / SUM_INNER / SUM_PP",
     "Source sums mirror get_PPR: SUM_ALL is every source, SUM_INNER drops Import, SUM_PP drops "
     "Import and Detritus. They are computed from each method's own columns, not from the "
     "NaN-padded union."),
    ("AGGREGATE source",
     "SPPR_1986 and SPPR_1995 return one un-attributed SPPR column. It is stored under the basal "
     "source AGGREGATE, and SUM_PP / ppr_pp_only are NaN because the value cannot be attributed "
     "to primary production alone."),
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
    MethodSpec("EwE_TE_EE",
               "EwE path enumeration with EE weighting -- the heavy one.",
               # silent=True for the same reason as the Monte-Carlo methods: silent=False routes
               # through tqdm.notebook, which raises ImportError('IProgress not found') outside
               # Jupyter and takes the whole method down.
               lambda m, mc_samples: (
                   _first(m.SPPR_EwE(TE_option="TE", use_EE=True, return_paths=True,
                                     silent=True)), {})),
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


# --------------------------------------------------------------------------- result container

@dataclass
class ModelTables:
    """Every table for one model, plus everything worth telling the user about the run."""
    model_label: str
    source_file: str
    groups: pd.DataFrame
    sppr: pd.DataFrame
    health: pd.DataFrame
    footprint: pd.DataFrame
    mc_diagnostics: pd.DataFrame
    notes: pd.DataFrame
    issues: list = field(default_factory=list)

    @property
    def sheets(self) -> dict:
        return {
            SHEET_GROUPS: self.groups,
            SHEET_SPPR: self.sppr,
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


def _ordered_sources(model: PPRCalculator) -> list:
    """Basal sources in reading order: PP, then detritus, then import."""
    ordered = (sorted(int(s) for s in model.get_PP_seq())
               + sorted(int(s) for s in model.get_DET_seq())
               + sorted(int(s) for s in model.get_Import_seq()))
    return ordered


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


def build_sppr_table(model: PPRCalculator, results: dict, method_keys: Sequence[str]) -> pd.DataFrame:
    """SPPR per group x basal source x method, with the three source sums per method.

    Columns are a (method, basal_source) MultiIndex over the union of basal sources, so methods
    are directly comparable; a source a given method does not resolve stays NaN.
    """
    groups = model.get_groups_df()
    index = pd.MultiIndex.from_arrays(
        [groups.index, [model.seq2name.get(int(s), str(s)) for s in groups.index]],
        names=["seq", "group_name"])

    sources = _ordered_sources(model)
    source_names = [model.seq2name.get(s, str(s)) for s in sources]
    needs_aggregate = any(not SPEC_BY_KEY[k].source_resolved for k in method_keys)
    column_sources = source_names + ([AGGREGATE_SOURCE] if needs_aggregate else [])

    blocks = {}
    for key in method_keys:
        spec = SPEC_BY_KEY[key]
        frame = pd.DataFrame(np.nan, index=index,
                             columns=list(column_sources) + list(SUM_COLUMNS), dtype=float)
        raw = results.get(key)
        if raw is not None:
            if spec.source_resolved:
                for seq, name in zip(sources, source_names):
                    if seq in raw.columns:
                        frame[name] = raw[seq].to_numpy(dtype=float)
            else:
                # single un-attributed column, whatever it is called
                frame[AGGREGATE_SOURCE] = raw.iloc[:, 0].to_numpy(dtype=float)
            for name, values in _source_sums(model, raw, spec.source_resolved).items():
                frame[name] = values.to_numpy(dtype=float)
        blocks[key] = frame

    table = pd.concat([blocks[k] for k in method_keys], axis=1,
                      keys=list(method_keys), names=["method", "basal_source"])
    table.index = index
    return table


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


def build_health_table(model: PPRCalculator, issues: list) -> pd.DataFrame:
    """diagnose_sppr under the fixed detritus config, one row per TE_option."""
    rows = {}
    for te in HEALTH_TE_OPTIONS:
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                report = model.diagnose_sppr(TE_option=te, **HEALTH_DET_CONFIG)
            for w in caught:
                issues.append(_issue("warning", f"[diagnose_sppr {te}] {w.message}"))
            for text in report.get("warnings", []):
                issues.append(_issue("health_warning", f"[{te}] {text}"))
            rows[te] = PPRCalculator._flatten_diagnostics(report)
        except Exception as exc:
            issues.append(_issue("health_failed", f"[{te}] {type(exc).__name__}: {exc}"))
            rows[te] = {"status": "ERROR", "error": f"{type(exc).__name__}: {exc}"}

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


def build_notes_table(method_keys: Sequence[str]) -> pd.DataFrame:
    notes = [{"topic": topic, "note": text} for topic, text in RUN_NOTES]
    notes += [{"topic": f"method: {k}", "note": SPEC_BY_KEY[k].description} for k in method_keys]
    return pd.DataFrame(notes).set_index("topic")


# --------------------------------------------------------------------------- orchestration

def build_model_tables(model: PPRCalculator, *, method_keys: Optional[Sequence[str]] = None,
                       mc_samples: int = DEFAULT_MC_SAMPLES,
                       model_label: str = "", source_file: str = "") -> ModelTables:
    """Run the selected SPPR methods on one model and assemble every table.

    A method that raises is recorded in `issues` and leaves its whole block NaN; it never aborts
    the run and its cells are never zero-filled.

    Args:
        model: the calculator to run.
        method_keys: which methods to run, defaulting to all of them (including the heavy
            SPPR_EwE). Pass FAST_METHOD_KEYS to skip the heavy one.
        mc_samples: draws per Monte-Carlo method.
        model_label: human-readable model name for the report.
        source_file: the JSON this model came from.

    Returns:
        ModelTables: the six sheets plus the issue list.
    """
    keys = list(ALL_METHOD_KEYS if method_keys is None else method_keys)
    unknown = [k for k in keys if k not in SPEC_BY_KEY]
    if unknown:
        raise KeyError(f"unknown method key(s) {unknown}; known keys: {list(ALL_METHOD_KEYS)}")

    issues: list = []
    results: dict = {}
    extras: dict = {}

    for key in keys:
        spec = SPEC_BY_KEY[key]
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                sppr, extra = spec.run(model, mc_samples)
            for w in caught:
                issues.append(_issue("warning", f"{w.category.__name__}: {w.message}", key))
            frame = _to_seq_frame(model, sppr)
            results[key] = frame
            extras[key] = extra

            values = frame.to_numpy(dtype=float)
            if np.any(values < -1e-10):
                issues.append(_issue(
                    "negative_sppr",
                    f"{int(np.sum(values < -1e-10))} negative SPPR cell(s), min "
                    f"{np.nanmin(values):.6g}", key))
            if not np.isfinite(values[~np.isnan(values)]).all():
                issues.append(_issue("non_finite_sppr", "contains infinite values", key))
        except Exception as exc:
            issues.append(_issue(
                "method_failed", f"{type(exc).__name__}: {exc}", key))
            issues.append(_issue("traceback", traceback.format_exc(), key))

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

    return ModelTables(
        model_label=model_label,
        source_file=source_file,
        groups=build_groups_table(model),
        sppr=build_sppr_table(model, results, keys),
        health=build_health_table(model, issues),
        footprint=build_footprint_table(model, results, keys),
        mc_diagnostics=build_mc_table(extras, keys),
        notes=build_notes_table(keys),
        issues=issues,
    )


def load_model(model_path: str) -> tuple[PPRCalculator, str]:
    """Load one Ecopath JSON into a calculator, returning it with a human-readable label."""
    model_data = ModelData(model_path)
    model = PPRCalculator.from_modeldata(
        model_data, underdetermined=True, zero_biomass_accum=False)
    label = f"{model_data.model_name} ({model_data.model_year})"
    return model, label


# --------------------------------------------------------------------------- excel writing

_COLUMN_WIDTH = 18
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
                      mc_samples: int = DEFAULT_MC_SAMPLES) -> str:
    """Build every table for one model JSON and write it to a single workbook.

    Returns:
        str: the path of the workbook written.
    """
    model, label = load_model(model_path)
    tables = build_model_tables(model, method_keys=method_keys, mc_samples=mc_samples,
                               model_label=label, source_file=os.path.basename(model_path))
    return write_tables_excel(tables, model_path, out_dir)


def write_tables_excel(tables: ModelTables, model_path: str,
                       out_dir: str = DEFAULT_OUT_DIR) -> str:
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    out_path = os.path.join(out_dir, Path(model_path).stem + ".xlsx")

    # (index columns, header rows) per sheet, so the reader and the freeze panes agree.
    layout = {
        SHEET_GROUPS: (1, 1),
        SHEET_SPPR: (2, 2),
        SHEET_HEALTH: (1, 1),
        SHEET_FOOTPRINT: (1, 1),
        SHEET_MC: (1, 1),
        SHEET_NOTES: (1, 1),
    }

    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        for sheet, table in tables.sheets.items():
            table.to_excel(writer, sheet_name=sheet)
            n_index, n_header = layout[sheet]
            _autoformat(writer, sheet, table, n_index, n_header)
    return out_path


SHEET_READ_LAYOUT = {
    SHEET_GROUPS: dict(index_col=0, header=0),
    SHEET_SPPR: dict(index_col=[0, 1], header=[0, 1]),
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
        table = pd.read_excel(path, sheet_name=sheet, **kwargs)
        if sheet == SHEET_SPPR:
            table.index.names = ["seq", "group_name"]
            table.columns.names = ["method", "basal_source"]
        out[sheet] = table
    return out


# --------------------------------------------------------------------------- run report

_ISSUE_HEADINGS = [
    ("method_failed", "METHODS THAT FAILED (their cells are empty, not zero)"),
    ("health_failed", "diagnose_sppr FAILURES"),
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

    lines.append(f"  health         : {entry.get('health_summary', 'n/a')}")
    lines.append(f"  methods ok     : {entry.get('n_ok', 0)}/{entry.get('n_methods', 0)}")
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


def write_run_report(entries: list, out_dir: str, method_keys: Sequence[str],
                     mc_samples: int) -> str:
    """Write the human-readable run report and return its path."""
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    path = os.path.join(out_dir, LOG_FILENAME)

    n_written = sum(1 for e in entries if e.get("output"))
    n_failed = len(entries) - n_written

    header = [
        "=" * 88,
        "SPPR EXPORT RUN REPORT",
        "=" * 88,
        f"generated      : {datetime.now().isoformat(timespec='seconds')}",
        f"models found   : {len(entries)}",
        f"workbooks ok   : {n_written}",
        f"models failed  : {n_failed}",
        f"methods run    : {len(method_keys)} -> {', '.join(method_keys)}",
        f"MC draws       : {mc_samples}",
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
                  silent: bool = False) -> dict:
    """Write one workbook per model JSON in `json_dir`, plus a run report.

    A model that cannot be loaded, or a method that raises, is recorded and the run continues.

    Args:
        json_dir: directory of Ecopath model JSONs.
        out_dir: where the workbooks and the report are written.
        method_keys: which methods to run; defaults to all of them.
        mc_samples: draws per Monte-Carlo method.
        silent: suppress the progress bar.

    Returns:
        dict: {'n_models', 'n_written', 'n_failed', 'report', 'entries'}.
    """
    keys = list(ALL_METHOD_KEYS if method_keys is None else method_keys)
    paths = sorted(str(p) for p in Path(json_dir).glob("*.json"))
    entries = []

    for model_path in tqdm(paths, disable=silent, desc="models"):
        entry = {"source_file": os.path.basename(model_path)}
        try:
            model, label = load_model(model_path)
            entry["model_label"] = label
            tables = build_model_tables(
                model, method_keys=keys, mc_samples=mc_samples,
                model_label=label, source_file=entry["source_file"])
            entry["output"] = write_tables_excel(tables, model_path, out_dir)
            entry["issues"] = tables.issues
            entry["n_methods"] = len(keys)
            entry["n_ok"] = len(keys) - len(
                [i for i in tables.issues if i["kind"] == "method_failed"])
            entry["health_summary"] = _health_summary(tables.health)
            entry["mc_summary"] = _mc_summary(tables.mc_diagnostics)
        except Exception as exc:
            entry["load_error"] = f"{type(exc).__name__}: {exc}"
        entries.append(entry)

    report = write_run_report(entries, out_dir, keys, mc_samples)
    n_written = sum(1 for e in entries if e.get("output"))
    return {
        "n_models": len(entries),
        "n_written": n_written,
        "n_failed": len(entries) - n_written,
        "report": report,
        "entries": entries,
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    json_dir = argv[0] if argv else DEFAULT_JSON_DIR
    out_dir = argv[1] if len(argv) > 1 else DEFAULT_OUT_DIR

    summary = run_directory(json_dir, out_dir)
    print(f"wrote {summary['n_written']}/{summary['n_models']} workbooks to {out_dir}")
    print(f"report: {summary['report']}")
    if summary["n_failed"]:
        print(f"{summary['n_failed']} model(s) could not be processed -- see the report")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
