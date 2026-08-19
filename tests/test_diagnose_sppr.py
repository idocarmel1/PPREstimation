"""Tests for PPRCalculator.diagnose_sppr.

Run with either:
    pytest tests/test_diagnose_sppr.py
    python tests/test_diagnose_sppr.py

Models used and why:
  * the multi-DET toy    -- the only mass-balanced model available, so the one case that can
                            legitimately return status 'OK'; also exercises k=2 detritus pools.
  * Black Sea 1960-1969  -- b = 0.83, i.e. converged but close to divergence, with a 305% PP
                            balance gap: the converged-but-untrustworthy case.
  * the toy with a flat, very low explicit TE -- forces b > 1 so the divergence FAIL path and
                            the negative-source count are actually exercised.
"""
from __future__ import annotations

import os
import sys
import warnings

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from PPRCalculator import DEFAULT_DIAGNOSTIC_THRESHOLDS, PPRCalculator  # noqa: E402

TOY = "real_models/ToyModels/900_900_Multi_DET_Toy_(2026).json"
BLACK_SEA = "real_models/new_EwE_jsons/62_10062_Black_Sea_(1960-1969).json"

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _calc(rel_path: str) -> PPRCalculator:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return PPRCalculator(os.path.join(ROOT, rel_path))


def _flat_te(pc: PPRCalculator, value: float) -> pd.DataFrame:
    """A constant TE matrix at `value`, with detritus rows left at 1 (as SPPR_new expects)."""
    DC = pc.get_DC(DET_as_PP=True, normalize=False)
    TE = pd.DataFrame(np.full(DC.shape, value), index=DC.index, columns=DC.columns)
    for d in pc.get_DET_seq():
        TE.loc[d, :] = 1.0
    return TE


@pytest.fixture(scope="module")
def toy() -> PPRCalculator:
    return _calc(TOY)


@pytest.fixture(scope="module")
def black_sea() -> PPRCalculator:
    return _calc(BLACK_SEA)


@pytest.fixture(scope="module")
def toy_report(toy) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return toy.diagnose_sppr()


@pytest.fixture(scope="module")
def black_sea_report(black_sea) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return black_sea.diagnose_sppr()


@pytest.fixture(scope="module")
def divergent_report(toy) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return toy.diagnose_sppr(TE=_flat_te(toy, 0.005))


# --------------------------------------------------------------------------- structure

def test_full_report_has_all_sections(toy_report):
    assert set(toy_report) == {
        "status", "model_input", "divergence", "balance", "footprint", "config", "warnings"
    }
    for section in ("model_input", "divergence", "balance"):
        assert toy_report[section]["status"] in {"OK", "WARN", "FAIL"}
    assert "status" not in toy_report["footprint"], "footprint must never be graded"


def test_short_report_drops_ungraded_sections(toy):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        report = toy.diagnose_sppr(short=True)
    assert set(report) == {"status", "model_input", "divergence", "balance"}


def test_top_status_is_worst_section(black_sea_report):
    rank = {"OK": 0, "WARN": 1, "FAIL": 2}
    section_worst = max(
        (black_sea_report[s]["status"] for s in ("model_input", "divergence", "balance")),
        key=lambda s: rank[s],
    )
    assert black_sea_report["status"] == section_worst


# --------------------------------------------------------------------------- healthy case

def test_toy_model_is_healthy(toy_report):
    """The toy is mass-balanced, acyclic among living groups and far from divergence."""
    assert toy_report["status"] == "OK"
    assert toy_report["model_input"]["is_model_balanced"] is True
    assert toy_report["divergence"]["b"] < DEFAULT_DIAGNOSTIC_THRESHOLDS["b_warn"]
    assert toy_report["divergence"]["b_converges"] is True
    assert toy_report["divergence"]["rho_living"] == pytest.approx(0.0)
    assert toy_report["divergence"]["n_negative_sources"] == 0
    assert toy_report["divergence"]["expect_negatives"] is False
    assert toy_report["warnings"] == []


def test_toy_reports_one_sppr_det_per_pool(toy, toy_report):
    det_seq = [int(d) for d in toy.get_DET_seq()]
    assert sorted(toy_report["divergence"]["sppr_det"]) == sorted(det_seq)
    assert len(det_seq) == 2, "the toy is the multi-DET fixture; k should be 2"


def test_sppr_det_matches_detritus_resolution_info(toy):
    """sppr_det is read off SPPR.loc[d, d]; it must equal _solve_det_scaling's x_vec."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        report = toy.diagnose_sppr()
    x_vec = toy.detritus_resolution_info["x_vec"]
    det_seq = [int(d) for d in toy.get_DET_seq()]
    for seq, x in zip(det_seq, x_vec):
        assert report["divergence"]["sppr_det"][seq] == pytest.approx(float(x))


def test_group_landmarks_have_the_documented_shape(toy, toy_report):
    for key in ("max_sppr_group", "max_tl_group"):
        rec = toy_report["divergence"][key]
        assert set(rec) == {"seq", "tl", "sppr"}, key
        assert isinstance(rec["seq"], int)
        assert rec["seq"] in [int(g) for g in toy.get_DC(DET_as_PP=True).index]


def test_max_sppr_group_really_is_the_maximum(toy):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        report, SPPR, _, _ = toy.diagnose_sppr(return_sppr=True)
    totals = SPPR.sum(axis=1)
    rec = report["divergence"]["max_sppr_group"]
    assert rec["sppr"] == pytest.approx(totals.max())
    assert rec["seq"] == int(totals.idxmax())


def test_max_tl_group_really_is_the_top_of_the_trophic_ordering(toy):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        report, SPPR, _, _ = toy.diagnose_sppr(return_sppr=True)
        tl = toy.get_TL(break_cycles=True, DET_as_PP=True)
    rec = report["divergence"]["max_tl_group"]
    assert rec["tl"] == pytest.approx(tl.max())
    assert rec["sppr"] == pytest.approx(SPPR.sum(axis=1)[rec["seq"]])


def test_trophic_levels_are_not_degenerate(black_sea_report):
    """Guards the self.TL trap: that attribute reads 1.0 for every group on real models, so
    the landmarks must come from get_TL, not from it."""
    assert black_sea_report["divergence"]["max_tl_group"]["tl"] > 1.5


def test_ee0_warning_names_the_groups(black_sea, black_sea_report):
    mi = black_sea_report["model_input"]
    assert mi["n_ee0"] > 0, "Black Sea is the EE=0 fixture"
    ee0_warnings = [w for w in black_sea_report["warnings"] if "EE=0" in w]
    assert len(ee0_warnings) == 1
    for seq in mi["ee0_groups"]:
        assert f"{seq} ({black_sea.seq2name.get(seq)})" in ee0_warnings[0]


def test_b_matches_recorded_rho_B(toy):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        report = toy.diagnose_sppr()
    assert report["divergence"]["b"] == pytest.approx(toy.detritus_resolution_info["rho_B"])


# ------------------------------------------------------- converged but untrustworthy

def test_black_sea_converges_but_fails(black_sea_report):
    div = black_sea_report["divergence"]
    assert div["b_converges"] is True, "b < 1: this model does converge"
    assert div["n_negative_sources"] == 0, "and produces no negatives"
    assert div["status"] == "WARN", "yet b > 0.7 must still be flagged"
    assert black_sea_report["balance"]["status"] == "FAIL"
    assert black_sea_report["balance"]["rel_gap"] > 1.0
    assert black_sea_report["status"] == "FAIL"


def test_black_sea_footprint_reports_all_three_ppr_views(black_sea_report):
    fp = black_sea_report["footprint"]
    assert fp["ppr_all"] >= fp["ppr_inner"] >= fp["ppr_pp_only"]
    assert fp["npp"] > 0
    assert fp["ppr2npp"] == pytest.approx(fp["ppr_inner"] / fp["npp"])
    assert fp["ppr2npp_pp_only"] == pytest.approx(fp["ppr_pp_only"] / fp["npp"])


# --------------------------------------------------------------------------- divergence

def test_forced_divergence_is_caught(divergent_report):
    div = divergent_report["divergence"]
    assert div["b"] > 1.0
    assert div["b_converges"] is False
    assert div["expect_negatives"] is True
    assert div["n_negative_sources"] >= 1
    assert div["status"] == "FAIL"
    assert divergent_report["status"] == "FAIL"
    assert any("diverges" in w for w in divergent_report["warnings"])


def test_divergence_is_monotone_in_te(toy):
    """Lower TE amplifies every path, so b must rise as TE falls."""
    bs = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for te in (0.05, 0.02, 0.01):
            bs.append(toy.diagnose_sppr(TE=_flat_te(toy, te), short=True)["divergence"]["b"])
    assert bs == sorted(bs), f"b should increase as TE decreases, got {bs}"
    assert bs[0] < 1.0 < bs[-1], "the sweep should straddle the divergence threshold"


def test_theta_damping_rescues_a_divergent_case(toy):
    """b is measured on diag(theta) @ B, so retention loss must lower it."""
    TE = _flat_te(toy, 0.01)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        closed = toy.diagnose_sppr(TE=TE, short=True)
        opened = toy.diagnose_sppr(TE=TE, short=True,
                                  det_open_mode="recycling_loss", det_theta=0.2)
    assert closed["divergence"]["b"] > 1.0
    assert opened["divergence"]["b"] < 1.0
    assert opened["divergence"]["b"] == pytest.approx(0.2 * closed["divergence"]["b"])


def test_collapse_mode_does_not_soften_b(toy):
    """det_collapse_mode is a remedy, not a diagnosis: b is measured pre-decision.

    On this toy the pooled fallback is *worse* than the coupled system (q_combined is
    dominated by the fast pool, so pooling raises b from 1.04 to 1.13), which makes 'auto'
    raise. The recovery re-solve must still report the same b as 'never'.
    """
    TE = _flat_te(toy, 0.01)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        never = toy.diagnose_sppr(TE=TE, det_collapse_mode="never")
        auto = toy.diagnose_sppr(TE=TE, det_collapse_mode="auto")
    assert never["divergence"]["b"] > 1.0
    assert never["divergence"]["solve_error"] is None
    assert never["config"]["would_pool"] is False
    # 'auto' pools, the pooled system also diverges, SPPR_new raises -- and the report still
    # measures the identical b instead of losing it.
    assert auto["divergence"]["solve_error"] is not None
    assert auto["divergence"]["b"] == pytest.approx(never["divergence"]["b"])
    assert auto["divergence"]["status"] == "FAIL"
    assert auto["balance"]["rel_gap"] is None, "no SPPR was produced, so nothing to balance"
    assert auto["config"]["would_pool"] is True, "auto must pool once b >= 1"


# --------------------------------------------------------------------------- input checks

# Editing catch also breaks Ecopath's production identity
# (production = catch + predation + growth + net_migration + M0), which would grade FAIL and
# mask the catch check under test. These thresholds isolate the catch branch.
NO_BALANCE_GRADING = {"model_balance_warn": np.inf, "model_balance_fail": np.inf}


def test_zero_catch_is_warned_without_invalidating_divergence():
    pc = _calc(TOY)
    pc.catch = pc.catch * 0.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        report = pc.diagnose_sppr(thresholds=NO_BALANCE_GRADING)
    mi = report["model_input"]
    assert mi["has_catch"] is False
    assert mi["total_catch"] == 0.0
    assert mi["status"] == "WARN"
    assert report["divergence"]["status"] == "OK", "a zero catch says nothing about convergence"
    assert report["footprint"]["ppr_all"] == pytest.approx(0.0)
    assert any("total catch is 0" in w for w in report["warnings"])


def test_negative_catch_is_warned():
    pc = _calc(TOY)
    pc.catch = pc.catch.copy()
    pc.catch.iloc[0] = -1.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        report = pc.diagnose_sppr(thresholds=NO_BALANCE_GRADING)
    assert report["model_input"]["n_negative_catch"] == 1
    assert report["model_input"]["status"] == "WARN"
    assert any("negative catch" in w for w in report["warnings"])


def test_dc_rows_sum_to_one_at_load_tolerance(toy_report):
    """ModelData.validate_DC guards 1e-3 at load, so any loaded model passes at 1e-3."""
    assert toy_report["model_input"]["dc_max_deviation"] <= 1e-3


def test_ee_flags_are_reported(black_sea_report):
    mi = black_sea_report["model_input"]
    assert mi["has_ee_issues"] is (mi["n_ee0"] > 0 or mi["n_ee_marginal"] > 0)
    assert mi["n_ee0"] == len(mi["ee0_groups"])
    assert mi["n_ee_marginal"] == len(mi["ee_marginal_groups"])


def test_model_balance_residuals_are_reported(black_sea_report):
    mi = black_sea_report["model_input"]
    assert mi["is_model_balanced"] is False
    assert max(mi["p_max_rel_residual"], mi["q_max_rel_residual"]) > \
        DEFAULT_DIAGNOSTIC_THRESHOLDS["model_balance_fail"]
    assert mi["status"] == "FAIL"


# --------------------------------------------------------------------------- options

def test_thresholds_override(black_sea):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        strict = black_sea.diagnose_sppr(short=True)
        loose = black_sea.diagnose_sppr(short=True, thresholds={"b_warn": 0.9})
    assert strict["divergence"]["status"] == "WARN"
    assert loose["divergence"]["status"] == "OK", "b=0.83 is below an overridden 0.9 warn"


def test_flat_report_is_single_level(toy):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        report = toy.diagnose_sppr(flat=True)
    assert not any(isinstance(v, (dict, list, tuple)) for v in report.values())
    assert "divergence_b" in report
    assert "balance_inflow" in report and "balance_outflow" in report
    assert "n_warnings" in report
    det_seq = [int(d) for d in toy.get_DET_seq()]
    for d in det_seq:
        assert f"divergence_sppr_det_{d}" in report


def test_flat_reports_concatenate_into_a_dataframe(toy, black_sea):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rows = [toy.diagnose_sppr(flat=True, short=True),
                black_sea.diagnose_sppr(flat=True, short=True)]
    df = pd.DataFrame(rows)
    assert len(df) == 2
    assert df["divergence_b"].notna().all()


def test_return_sppr_gives_back_the_solve(toy):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        report, SPPR, A, L = toy.diagnose_sppr(return_sppr=True)
    assert isinstance(SPPR, pd.DataFrame) and isinstance(A, pd.DataFrame)
    assert L.shape == A.shape
    det_seq = [int(d) for d in toy.get_DET_seq()]
    for d in det_seq:
        assert float(SPPR.loc[d, d]) == pytest.approx(report["divergence"]["sppr_det"][d])


def test_te_option_te_has_no_recycling_matrix(black_sea):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        report = black_sea.diagnose_sppr(TE_option="TE")
    assert report["divergence"]["b"] == 0.0
    assert report["config"]["TE_option"] == "TE"
    assert any("no detritus recycling matrix" in w for w in report["warnings"])


def test_config_echoes_what_was_evaluated(toy):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        report = toy.diagnose_sppr(det_open_mode="recycling_loss", det_theta=0.5)
    cfg = report["config"]
    assert cfg["det_open_mode"] == "recycling_loss"
    assert cfg["det_theta"] == 0.5
    assert cfg["explicit_TE"] is False
    assert cfg["model"] is not None


# --------------------------------------------------------------------------- error paths

@pytest.mark.parametrize("kwargs", [
    {"TE_option": "nonsense"},
    {"det_open_mode": "nonsense"},
    {"det_collapse_mode": "nonsense"},
])
def test_caller_mistakes_raise(toy, kwargs):
    with pytest.raises(ValueError):
        toy.diagnose_sppr(**kwargs)


def test_does_not_raise_when_the_pooled_fallback_diverges(toy):
    """det_collapse_mode='always' on a divergent model makes _collapse_det_scaling raise.

    diagnose_sppr must report FAIL instead of propagating, and the recovery re-solve must
    still measure b so the report explains the failure rather than reporting nothing.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        report = toy.diagnose_sppr(TE=_flat_te(toy, 0.001), det_collapse_mode="always")
    assert report["status"] == "FAIL"
    assert report["divergence"]["status"] == "FAIL"
    assert report["divergence"]["solve_error"] is not None
    assert any("SPPR_new failed" in w for w in report["warnings"])
    # recovered diagnostics
    assert report["divergence"]["b"] > 1.0
    assert report["divergence"]["rho_living"] is not None
    assert any("diagnostic re-solve" in w for w in report["warnings"])
    # but nothing the caller's own configuration could have produced
    assert report["balance"]["status"] == "FAIL"
    assert report["balance"]["inflow"] is None
    assert report["footprint"]["ppr_all"] is None
    assert report["config"]["method"] is None
    assert report["config"]["would_pool"] is True


def test_b_is_none_only_when_even_the_recovery_solve_fails(toy):
    """b == 1 exactly is singular for the coupled solve too, so nothing can be measured."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        report = toy.diagnose_sppr(short=True)
    # sanity: the healthy path always has a b
    assert report["divergence"]["b"] is not None


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
