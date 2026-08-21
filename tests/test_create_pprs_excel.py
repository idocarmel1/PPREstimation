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


# --------------------------------------------------------------------- sppr_table

def test_sppr_table_columns_are_a_method_source_multiindex(cheap_tables):
    sppr = cheap_tables.sppr
    assert isinstance(sppr.columns, pd.MultiIndex)
    assert sppr.columns.nlevels == 2
    assert sppr.columns.names == ["method", "basal_source"]
    assert set(sppr.columns.get_level_values("method")) == set(CHEAP)


def test_sppr_table_rows_are_every_group_as_seq_and_name(cheap_tables, toy_model):
    sppr = cheap_tables.sppr
    assert sppr.index.names == ["seq", "group_name"]
    assert list(sppr.index.get_level_values("seq")) == list(toy_model.get_groups_df().index)
    assert "Phytoplankton" in sppr.index.get_level_values("group_name")


def test_each_method_carries_the_three_sum_columns(cheap_tables):
    for key in CHEAP:
        cols = set(cheap_tables.sppr[key].columns)
        assert {"SUM_PP", "SUM_INNER", "SUM_ALL"} <= cols, (key, cols)


def test_sums_follow_get_PPR_source_semantics(cheap_tables, toy_model):
    """SUM_ALL = every source, SUM_INNER = drop Import, SUM_PP = drop Import and DET."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        raw, _, _ = toy_model.SPPR_new(TE_option="GE")

    raw = raw.reindex(toy_model.get_groups_df().index)
    det = list(toy_model.get_DET_seq())
    imp = list(toy_model.get_Import_seq())

    block = cheap_tables.sppr["new_GE"]
    expected_all = raw.sum(axis=1)
    expected_inner = raw.drop(columns=[c for c in imp if c in raw.columns]).sum(axis=1)
    expected_pp = raw.drop(columns=[c for c in det + imp if c in raw.columns]).sum(axis=1)

    assert np.allclose(block["SUM_ALL"].values, expected_all.values)
    assert np.allclose(block["SUM_INNER"].values, expected_inner.values)
    assert np.allclose(block["SUM_PP"].values, expected_pp.values)


def test_a_source_a_method_does_not_resolve_is_nan_not_zero(cheap_tables):
    """SPPR_2015 returns no detritus columns; those cells must be NaN, never a silent 0."""
    block = cheap_tables.sppr["SPPR_2015"]
    det_cols = [c for c in block.columns if "detritus" in str(c).lower()]
    assert det_cols, f"expected detritus source columns in the union, got {list(block.columns)}"
    assert block[det_cols].isna().all().all(), (
        "a source the method does not resolve must be NaN, not 0")


def test_non_source_resolved_method_is_marked_aggregate_with_nan_sum_pp(cheap_tables):
    """SPPR_1986 returns one un-attributed 'sppr' column, so PP cannot be separated out."""
    block = cheap_tables.sppr["SPPR_1986"]
    assert cpe.AGGREGATE_SOURCE in block.columns
    assert block["SUM_PP"].isna().all(), "an un-attributed SPPR must not claim a PP-only sum"
    assert np.allclose(block["SUM_ALL"].values, block[cpe.AGGREGATE_SOURCE].values)


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

    block = tables.sppr["SPPR_2015"]
    assert block.isna().all().all(), "a failed method must leave empty cells, not zeros"
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

def test_round_trip_preserves_values_and_multiindex(toy_model, tmp_path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        path = cpe.write_model_excel(TOY, out_dir=str(tmp_path), method_keys=CHEAP)
        loaded = cpe.read_pprs_excel(path)
        original = cpe.build_model_tables(toy_model, method_keys=CHEAP)

    assert "sppr_table" in loaded and "groups_df" in loaded
    back = loaded["sppr_table"]

    assert isinstance(back.columns, pd.MultiIndex)
    assert back.columns.names == ["method", "basal_source"]
    assert back.index.names == ["seq", "group_name"]
    assert list(back.columns) == list(original.sppr.columns)
    np.testing.assert_allclose(
        back.values.astype(float), original.sppr.values.astype(float), equal_nan=True)


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


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
