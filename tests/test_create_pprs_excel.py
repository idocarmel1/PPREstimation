"""Tests for create_PPRS_excel.py -- the per-model SPPR/PPR Excel exporter and its reader.

Run with either:
    pytest tests/test_create_pprs_excel.py
    python tests/test_create_pprs_excel.py

Everything here runs on the multi-DET toy (6 groups, 2 detritus pools) with a deliberately small
method subset. SPPR_EwE is never invoked: it is the one genuinely heavy method, so it is excluded
from ALL the fast paths and only the full-sweep run touches it.

Toy layout, relied on below:
    PP=[1] 'Phytoplankton'   DET=[4,5] 'Fast/Slow detritus'   Import=[6] 'diet_import'
    Regular=[2,3] 'Herbivore', 'Detritivore'
"""
from __future__ import annotations

import os
import sys
import warnings

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import create_PPRS_excel as cpe  # noqa: E402
from ModelData import ModelData  # noqa: E402
from PPRCalculator import PPRCalculator  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOY_DIR = os.path.join(ROOT, "real_models", "ToyModels")
TOY = os.path.join(TOY_DIR, "900_900_Multi_DET_Toy_(2026).json")

# Cheap, source-resolved + non-source-resolved + no-DET-column methods, so the shape logic is
# fully exercised without paying for the symbolic or Monte-Carlo paths.
CHEAP = ("SPPR_1986", "new_GE", "SPPR_2015")


@pytest.fixture(scope="module")
def toy_model() -> PPRCalculator:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return PPRCalculator.from_modeldata(
            ModelData(TOY), underdetermined=True, zero_biomass_accum=False)


@pytest.fixture(scope="module")
def cheap_tables(toy_model):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return cpe.build_model_tables(toy_model, method_keys=CHEAP)


# --------------------------------------------------------------------- method registry

def test_registry_holds_the_nineteen_specified_methods():
    assert len(cpe.ALL_METHOD_KEYS) == 19, cpe.ALL_METHOD_KEYS
    assert len(set(cpe.ALL_METHOD_KEYS)) == 19, "method keys must be unique"
    # the six symbolic combinations must all be present
    symbolic = [k for k in cpe.ALL_METHOD_KEYS if k.startswith("sym_")]
    assert len(symbolic) == 6, symbolic
    assert {"MC_new_GE", "MC_new_TE_EEfix"} <= set(cpe.ALL_METHOD_KEYS)


def test_heavy_methods_are_excluded_from_the_fast_key_set():
    assert "EwE_TE_EE" in cpe.HEAVY_METHOD_KEYS
    assert "EwE_TE_EE" not in cpe.FAST_METHOD_KEYS
    assert set(cpe.FAST_METHOD_KEYS) | set(cpe.HEAVY_METHOD_KEYS) == set(cpe.ALL_METHOD_KEYS)


# --------------------------------------------------------------------- the three sppr sheets

def test_there_are_three_flat_sppr_sheets(cheap_tables):
    sheets = cheap_tables.sheets
    for name in (cpe.SHEET_SPPR_PP, cpe.SHEET_SPPR_INNER, cpe.SHEET_SPPR_ALL):
        assert name in sheets, name


@pytest.mark.parametrize("attr", ["sppr_pp", "sppr_inner", "sppr_all"])
def test_sppr_sheets_are_flat_groups_by_methods(cheap_tables, toy_model, attr):
    table = getattr(cheap_tables, attr)
    assert not isinstance(table.columns, pd.MultiIndex), "columns must be flat"
    assert not isinstance(table.index, pd.MultiIndex), "index must be flat"
    assert table.index.name == "seq"
    assert list(table.columns) == ["group_name"] + list(CHEAP)
    assert list(table.index) == list(toy_model.get_groups_df().index)
    assert table.loc[1, "group_name"] == "Phytoplankton"


def test_sheets_hold_the_matching_source_sums(cheap_tables, toy_model):
    """all = every source, inner = drop Import, PP = drop Import and DET."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        raw, _, _ = toy_model.SPPR_new(TE_option="GE")

    raw = raw.reindex(toy_model.get_groups_df().index)
    det = list(toy_model.get_DET_seq())
    imp = list(toy_model.get_Import_seq())

    expected_all = raw.sum(axis=1)
    expected_inner = raw.drop(columns=[c for c in imp if c in raw.columns]).sum(axis=1)
    expected_pp = raw.drop(columns=[c for c in det + imp if c in raw.columns]).sum(axis=1)

    assert np.allclose(cheap_tables.sppr_all["new_GE"].values, expected_all.values)
    assert np.allclose(cheap_tables.sppr_inner["new_GE"].values, expected_inner.values)
    assert np.allclose(cheap_tables.sppr_pp["new_GE"].values, expected_pp.values)


def test_method_without_detritus_columns_has_equal_inner_and_pp(cheap_tables):
    """SPPR_2015 resolves only PP and Import, so dropping Import already leaves PP alone."""
    inner = cheap_tables.sppr_inner["SPPR_2015"]
    pp = cheap_tables.sppr_pp["SPPR_2015"]
    assert np.allclose(inner.values, pp.values, equal_nan=True)


def test_non_source_resolved_method_is_nan_in_the_pp_sheet_only(cheap_tables):
    """SPPR_1986 returns one un-attributed column, so PP cannot be separated out of it."""
    assert cheap_tables.sppr_pp["SPPR_1986"].isna().all(), (
        "an un-attributed SPPR must not claim a PP-only value")
    assert cheap_tables.sppr_all["SPPR_1986"].notna().any()
    assert np.allclose(cheap_tables.sppr_all["SPPR_1986"].values,
                       cheap_tables.sppr_inner["SPPR_1986"].values)


# --------------------------------------------------------------------- footprint

def test_footprint_columns_match_the_diagnose_sppr_footprint_keys(cheap_tables, toy_model):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        keys = list(toy_model.diagnose_sppr(TE_option="GE")["footprint"])

    assert list(cheap_tables.footprint.columns) == keys


def test_footprint_index_is_the_selected_method_names(cheap_tables):
    assert list(cheap_tables.footprint.index) == list(CHEAP)


def test_footprint_ratio_matches_get_PPR2NPP_ratio(cheap_tables, toy_model):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        raw, _, _ = toy_model.SPPR_new(TE_option="GE")
        expected = toy_model.get_PPR2NPP_ratio(raw)

    assert cheap_tables.footprint.loc["new_GE", "ppr2npp"] == pytest.approx(expected)


# --------------------------------------------------------------------- model_health

def test_health_table_has_one_row_per_TE_option(cheap_tables):
    health = cheap_tables.health
    assert set(cpe.HEALTH_TE_OPTIONS) <= set(health.index)


def test_health_table_echoes_the_locked_detritus_config(cheap_tables):
    health = cheap_tables.health
    cfg = {c: health[c] for c in health.columns if c.startswith("config_")}
    assert cfg, f"expected echoed config columns, got {list(health.columns)}"
    assert (health["config_det_open_mode"] == "none").all()
    assert (health["config_det_theta"] == 1.0).all()
    assert (health["config_det_external_sppr"] == 0.0).all()
    assert (health["config_det_collapse_mode"] == "never").all()


# --------------------------------------------------------------------- groups_df

def test_groups_table_keeps_the_parameter_columns_and_adds_group_types(cheap_tables):
    groups = cheap_tables.groups
    for col in ("group_name", "tl", "ge", "ee", "catch", "biomass", "pb", "qb", "p", "q"):
        assert col in groups.columns, col
    assert "group_type" in groups.columns
    assert set(groups["group_type"]) <= {"PP", "DET", "Import", "Regular"}
    assert (groups.loc[1, "group_type"] == "PP")


# --------------------------------------------------------------------- failures

def test_a_failing_method_is_recorded_and_leaves_its_cells_empty(toy_model, monkeypatch):
    def boom(*args, **kwargs):
        raise RuntimeError("synthetic solver failure")

    monkeypatch.setattr(toy_model, "SPPR_2015", boom)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tables = cpe.build_model_tables(toy_model, method_keys=CHEAP)

    failed = [i for i in tables.issues if i["kind"] == "method_failed"]
    assert [i["method"] for i in failed] == ["SPPR_2015"]
    assert "synthetic solver failure" in failed[0]["detail"]

    for table in (tables.sppr_pp, tables.sppr_inner, tables.sppr_all):
        assert table["SPPR_2015"].isna().all(), (
            "a failed method must leave empty cells, not zeros")
    assert tables.footprint.loc["SPPR_2015"].isna().all()


# --------------------------------------------------------------------- monte carlo

def test_monte_carlo_diagnostics_are_collected(toy_model):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tables = cpe.build_model_tables(
            toy_model, method_keys=("MC_new_GE",), mc_samples=2)

    mc = tables.mc_diagnostics
    assert "MC_new_GE" in mc.index
    for col in ("n_accepted", "n_rejected_negative", "n_rejected_diverged", "reject_frac"):
        assert col in mc.columns, col
    assert mc.loc["MC_new_GE", "n_accepted"] <= 2


# --------------------------------------------------------------------- excel round trip

@pytest.mark.parametrize("sheet,attr", [
    ("sppr_PP", "sppr_pp"),
    ("sppr_inner", "sppr_inner"),
    ("sppr_all", "sppr_all"),
])
def test_round_trip_preserves_each_sppr_sheet_flat(toy_model, tmp_path, sheet, attr):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        path = cpe.write_model_excel(TOY, out_dir=str(tmp_path), method_keys=CHEAP)
        loaded = cpe.read_pprs_excel(path)
        original = cpe.build_model_tables(toy_model, method_keys=CHEAP)

    assert "groups_df" in loaded
    back = loaded[sheet]
    expected = getattr(original, attr)

    assert not isinstance(back.columns, pd.MultiIndex)
    assert not isinstance(back.index, pd.MultiIndex)
    assert back.index.name == "seq"
    assert list(back.columns) == list(expected.columns)
    assert list(back["group_name"]) == list(expected["group_name"])
    np.testing.assert_allclose(
        back[list(CHEAP)].values.astype(float),
        expected[list(CHEAP)].values.astype(float), equal_nan=True)


def test_round_trip_preserves_the_footprint_table(toy_model, tmp_path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        path = cpe.write_model_excel(TOY, out_dir=str(tmp_path), method_keys=CHEAP)
        loaded = cpe.read_pprs_excel(path)
        original = cpe.build_model_tables(toy_model, method_keys=CHEAP)

    back = loaded["footprint"]
    assert list(back.index) == list(original.footprint.index)
    assert list(back.columns) == list(original.footprint.columns)
    np.testing.assert_allclose(
        back.values.astype(float), original.footprint.values.astype(float), equal_nan=True)


# --------------------------------------------------------------------- directory run

def test_run_directory_writes_one_workbook_per_model_and_a_log(tmp_path):
    out = tmp_path / "out"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        summary = cpe.run_directory(TOY_DIR, out_dir=str(out), method_keys=CHEAP)

    written = sorted(p.name for p in out.glob("*.xlsx"))
    assert len(written) == 1, written
    assert summary["n_models"] == 1
    assert summary["n_written"] == 1

    log = out / cpe.LOG_FILENAME
    assert log.exists(), "run_directory must leave a human-readable report"
    text = log.read_text(encoding="utf-8")
    assert "SPPR EXPORT RUN REPORT" in text
    assert "900_900_Multi_DET_Toy" in text


def test_log_reports_method_failures_per_model(tmp_path, monkeypatch):
    real = cpe.build_model_tables

    def with_one_failure(model, **kwargs):
        monkeypatch.setattr(model, "SPPR_2015",
                            lambda *a, **k: (_ for _ in ()).throw(RuntimeError("synthetic boom")))
        return real(model, **kwargs)

    monkeypatch.setattr(cpe, "build_model_tables", with_one_failure)
    out = tmp_path / "out"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cpe.run_directory(TOY_DIR, out_dir=str(out), method_keys=CHEAP)

    text = (out / cpe.LOG_FILENAME).read_text(encoding="utf-8")
    assert "synthetic boom" in text
    assert "SPPR_2015" in text


def test_run_directory_records_a_model_that_cannot_be_loaded(tmp_path):
    bad_dir = tmp_path / "models"
    bad_dir.mkdir()
    (bad_dir / "broken.json").write_text("{not valid json", encoding="utf-8")

    out = tmp_path / "out"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        summary = cpe.run_directory(str(bad_dir), out_dir=str(out), method_keys=CHEAP)

    assert summary["n_written"] == 0
    assert summary["n_failed"] == 1
    text = (out / cpe.LOG_FILENAME).read_text(encoding="utf-8")
    assert "broken.json" in text


# --------------------------------------------------------------------- per-method timeout

# A budget no round trip can meet: the worker still has to receive the task, solve the whole
# system and answer over a pipe, which is milliseconds at the very best. Nothing about the toy
# model has to be slow for this to be deterministic -- the IPC alone outlasts the budget.
IMPOSSIBLE_BUDGET = 0.001


@pytest.fixture
def one_health_row(monkeypatch):
    """Trim model_health to a single TE option: each row costs a worker restart under a timeout."""
    monkeypatch.setattr(cpe, "HEALTH_TE_OPTIONS", ("GE",))


def test_a_method_over_its_budget_is_a_timeout_and_not_a_failure(toy_model, one_health_row):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tables = cpe.build_model_tables(toy_model, method_keys=("new_GE",),
                                        method_timeout=IMPOSSIBLE_BUDGET)

    assert tables.method_status["new_GE"] == "timeout"
    timed_out = [i for i in tables.issues if i["kind"] == "method_timeout"]
    assert [i["method"] for i in timed_out] == ["new_GE"]
    assert "still running after" in timed_out[0]["detail"]
    # a timeout is reported as its own thing, never folded into the failures
    assert not [i for i in tables.issues if i["kind"] == "method_failed"]


def test_a_timed_out_method_leaves_empty_cells_not_zeros(toy_model, one_health_row):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tables = cpe.build_model_tables(toy_model, method_keys=("new_GE",),
                                        method_timeout=IMPOSSIBLE_BUDGET)

    for table in (tables.sppr_pp, tables.sppr_inner, tables.sppr_all):
        assert table["new_GE"].isna().all()
    assert tables.footprint.loc["new_GE"].isna().all()
    # the workbook itself has to say which of NaN's several meanings this one is
    assert tables.notes.loc["method: new_GE", "status"] == "timeout"


def test_a_diagnose_sppr_row_over_its_budget_is_marked_TIMEOUT(toy_model, one_health_row):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tables = cpe.build_model_tables(toy_model, method_keys=("new_GE",),
                                        method_timeout=IMPOSSIBLE_BUDGET)

    assert tables.health.loc["GE", "status"] == "TIMEOUT"
    assert [i["kind"] for i in tables.issues if i["kind"] == "health_timeout"]


def test_the_worker_gives_the_same_numbers_as_running_in_process(cheap_tables, toy_model):
    """cheap_tables ran under the default budget, i.e. in a worker. Inline must match it."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        inline = cpe.build_model_tables(toy_model, method_keys=CHEAP, method_timeout=None)

    assert inline.method_status == cheap_tables.method_status
    for attr in ("sppr_pp", "sppr_inner", "sppr_all", "footprint"):
        pd.testing.assert_frame_equal(getattr(inline, attr), getattr(cheap_tables, attr))


def test_a_method_that_finishes_in_time_is_untouched_by_the_budget(cheap_tables):
    assert set(cheap_tables.method_status) == set(CHEAP)
    assert all(v == "ok" for v in cheap_tables.method_status.values()), cheap_tables.method_status
    assert not [i for i in cheap_tables.issues if i["kind"] == "method_timeout"]


def test_run_directory_counts_timeouts_and_names_them_in_the_report(tmp_path, one_health_row):
    out = tmp_path / "out"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        summary = cpe.run_directory(TOY_DIR, out_dir=str(out), method_keys=("new_GE",),
                                    method_timeout=IMPOSSIBLE_BUDGET)

    # a timeout must not cost the workbook: every other sheet is still written
    assert summary["n_written"] == 1
    assert summary["n_failed"] == 0
    assert summary["n_timeout"] == 1

    text = (out / cpe.LOG_FILENAME).read_text(encoding="utf-8")
    assert "METHODS THAT TIMED OUT" in text
    assert "new_GE" in text
    assert "time budget" in text


@pytest.mark.parametrize("argv,expected", [
    (["--timeout", "42"], 42.0),
    (["--timeout=90"], 90.0),
    (["-t", "7.5"], 7.5),
    (["--timeout", "none"], None),
    (["--timeout=off"], None),
    (["-t", "0"], None),
    ([], cpe.DEFAULT_METHOD_TIMEOUT),
])
def test_the_timeout_flag_is_parsed_and_removed_from_argv(argv, expected):
    rest = list(argv)
    assert cpe._pop_timeout_flag(rest) == expected
    assert rest == [], "the flag and its value must be consumed"


def test_the_timeout_flag_leaves_the_positional_arguments_alone():
    argv = ["real_models/ToyModels", "out", "--timeout", "12"]
    assert cpe._pop_timeout_flag(argv) == 12.0
    assert argv == ["real_models/ToyModels", "out"]


def test_an_unusable_worker_falls_back_to_running_in_process(toy_model, monkeypatch,
                                                             one_health_row):
    """Losing the worker costs the time budget, never the results."""
    class _Doomed(cpe._TaskRunner):
        def _start(self):
            raise OSError("synthetic spawn failure")

    monkeypatch.setattr(cpe, "_TaskRunner", _Doomed)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tables = cpe.build_model_tables(toy_model, method_keys=("new_GE",), method_timeout=60)

    assert tables.method_status["new_GE"] == "ok"
    assert tables.sppr_all["new_GE"].notna().any()
    unavailable = [i for i in tables.issues if i["kind"] == "timeout_unavailable"]
    assert unavailable and "synthetic spawn failure" in unavailable[0]["detail"]


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
