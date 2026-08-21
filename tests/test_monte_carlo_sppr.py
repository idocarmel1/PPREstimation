"""Tests for the method_kwargs / exclude_diverged / return_diagnostics additions to
PPRCalculator.monte_carlo_SPPR.

Run with either:
    pytest tests/test_monte_carlo_sppr.py
    python tests/test_monte_carlo_sppr.py

Models used and why:
  * the multi-DET toy    -- healthy (b = 0.069), so nothing is ever rejected for divergence;
                            used for the cheap structural tests (return arity, diagnostics keys).
  * Black Sea 1960-1969  -- b = 0.83 at TE_option='GE', i.e. close enough to divergence that
                            resampling TE pushes some draws over b = 1 and others not, which is
                            the only way to exercise a genuine accept/reject mix. It is also the
                            one model on hand where fix_EE_0_cases visibly changes SPPR_new's
                            result (under TE_option='TE'), which is what makes the
                            method_kwargs forwarding observable rather than a no-op.

Every test seeds numpy first: monte_carlo_SPPR draws its TE samples through scipy's gamma,
which reads the numpy global random state, so seeding makes the accept/reject split reproducible.
Sample counts are kept tiny (2-6) -- these tests must never become a heavy run.
"""
from __future__ import annotations

import os
import sys
import warnings

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from PPRCalculator import PPRCalculator  # noqa: E402

TOY = "real_models/ToyModels/900_900_Multi_DET_Toy_(2026).json"
BLACK_SEA = "real_models/new_EwE_jsons/62_10062_Black_Sea_(1960-1969).json"

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Resampling band that produces a genuine mix of diverged / accepted draws on Black Sea.
MIXED_ERR, MIXED_CUT = 30, 40

DIAGNOSTIC_KEYS = {
    "n_accepted",
    "n_rejected_negative",
    "n_rejected_diverged",
    "reject_frac_negative",
    "reject_frac_diverged",
    "per_sample_reason",
}


def _calc(rel_path: str) -> PPRCalculator:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return PPRCalculator(os.path.join(ROOT, rel_path))


@pytest.fixture(scope="module")
def toy() -> PPRCalculator:
    return _calc(TOY)


@pytest.fixture(scope="module")
def black_sea() -> PPRCalculator:
    return _calc(BLACK_SEA)


# ------------------------------------------------------------------ method_kwargs forwarding

def test_method_kwargs_forwards_fix_EE_0_cases_to_the_solver(black_sea):
    """method_kwargs must actually reach SPPR_new, not be silently dropped.

    fix_EE_0_cases changes Black Sea's SPPR_new result under TE_option='TE' (the EE=0 groups
    get their severed TE rows repaired), so forwarding it has to move the Monte-Carlo mean.
    """
    kw = dict(n_samples=2, TE_error_percent=5, TE_error_cut_percent=10,
              TE_option="TE", silent=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        np.random.seed(0)
        fixed, _, _, _, _ = black_sea.monte_carlo_SPPR(
            method_kwargs={"fix_EE_0_cases": True}, **kw)
        np.random.seed(0)
        unfixed, _, _, _, _ = black_sea.monte_carlo_SPPR(
            method_kwargs={"fix_EE_0_cases": False}, **kw)

    assert not np.allclose(fixed.values, unfixed.values, equal_nan=True), (
        "method_kwargs={'fix_EE_0_cases': ...} did not change the result -- it was not forwarded"
    )


def test_unknown_method_kwarg_raises_naming_the_offending_key(toy):
    """A typo must fail loudly with the key named, not be swallowed or hit the solver."""
    with pytest.raises(TypeError, match="fix_EE_cases"):
        toy.monte_carlo_SPPR(n_samples=2, method_kwargs={"fix_EE_cases": True}, silent=True)


@pytest.mark.parametrize("key,value", [
    ("TE_option", "GE"),
    ("TE", None),
    ("DET_TE_vals", 1),
    ("det_theta", 1.0),
])
def test_method_kwargs_collision_with_wrapper_owned_key_raises(toy, key, value):
    """Keys the wrapper already controls must not be settable twice from two places."""
    with pytest.raises(ValueError, match=key):
        toy.monte_carlo_SPPR(n_samples=2, method_kwargs={key: value}, silent=True)


def test_method_kwargs_defaults_to_no_change(toy):
    """Omitting method_kwargs must behave exactly like passing an empty dict."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        np.random.seed(0)
        without, _, _, _, _ = toy.monte_carlo_SPPR(n_samples=2, silent=True)
        np.random.seed(0)
        empty, _, _, _, _ = toy.monte_carlo_SPPR(n_samples=2, method_kwargs={}, silent=True)

    assert np.allclose(without.values, empty.values, equal_nan=True)


# ------------------------------------------------------------------ exclude_diverged

def test_exclude_diverged_rejects_draws_that_fail_the_divergence_grade(black_sea):
    """With exclude_diverged=True, FAIL-divergence draws must be dropped from the average."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        np.random.seed(0)
        _, _, _, _, _, diag = black_sea.monte_carlo_SPPR(
            n_samples=6, TE_error_percent=MIXED_ERR, TE_error_cut_percent=MIXED_CUT,
            TE_option="GE", exclude_diverged=True, return_diagnostics=True, silent=True)

    assert diag["n_rejected_diverged"] > 0, "no draw was rejected for divergence"
    assert diag["n_accepted"] > 0, "every draw was rejected -- config is not a usable mix"
    assert (diag["n_accepted"] + diag["n_rejected_negative"]
            + diag["n_rejected_diverged"]) == 6


def test_exclude_diverged_never_accepts_a_draw_the_ungated_run_rejected(black_sea):
    """The gate may only ever remove draws, and its rejections must be labelled 'diverged'.

    Note it is NOT strictly additional in general: on Black Sea every FAIL-divergence draw also
    happens to produce a negative SPPR, so the gated and ungated accepted sets coincide here and
    only the attribution changes ('diverged' rather than 'negative'). The invariant worth pinning
    is therefore the subset relation, not a smaller count -- a draw the ungated run rejected must
    never come back once the stricter gate is on.
    """
    kw = dict(n_samples=6, TE_error_percent=MIXED_ERR, TE_error_cut_percent=MIXED_CUT,
              TE_option="GE", return_diagnostics=True, silent=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        np.random.seed(0)
        _, _, _, _, _, gated = black_sea.monte_carlo_SPPR(exclude_diverged=True, **kw)
        np.random.seed(0)
        _, _, _, _, _, ungated = black_sea.monte_carlo_SPPR(exclude_diverged=False, **kw)

    assert ungated["n_rejected_diverged"] == 0, (
        "exclude_diverged=False must not reject anything for divergence")

    gated_ok = {i for i, r in enumerate(gated["per_sample_reason"]) if r == "accepted"}
    ungated_ok = {i for i, r in enumerate(ungated["per_sample_reason"]) if r == "accepted"}
    assert gated_ok <= ungated_ok, "the divergence gate resurrected a previously rejected draw"
    assert gated["n_accepted"] <= ungated["n_accepted"]
    assert gated["n_rejected_diverged"] > 0, "config no longer exercises the divergence gate"


def test_exclude_diverged_defaults_to_off(toy):
    """Default behaviour must be unchanged from before this feature existed."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        np.random.seed(0)
        _, _, _, _, _, diag = toy.monte_carlo_SPPR(
            n_samples=2, return_diagnostics=True, silent=True)

    assert diag["n_rejected_diverged"] == 0


def test_exclude_diverged_rejects_symbolic_kind(toy):
    """diagnose_sppr only grades SPPR_new, so the gate is meaningless for kind='symbolic'."""
    with pytest.raises(ValueError, match="symbolic"):
        toy.monte_carlo_SPPR(n_samples=2, kind="symbolic", exclude_diverged=True, silent=True)


def test_all_draws_rejected_warns_and_returns_nan_instead_of_crashing(black_sea, monkeypatch):
    """Zero accepted draws must degrade to an explicit warning + NaN, not a bare RuntimeWarning.

    No real TE band rejects every draw (clipping keeps at least one healthy), so the FAIL verdict
    is stubbed here. The assertions are on monte_carlo_SPPR's own return value, not on the stub.
    """
    real = black_sea.diagnose_sppr

    def always_diverged(*args, **kwargs):
        report, sppr, A, L = real(*args, **{**kwargs, "return_sppr": True})
        report["divergence"] = {**report["divergence"], "status": "FAIL"}
        return (report, sppr, A, L) if kwargs.get("return_sppr") else report

    monkeypatch.setattr(black_sea, "diagnose_sppr", always_diverged)

    np.random.seed(0)
    with pytest.warns(UserWarning, match="accepted"):
        mean, samples, reject_frac, _, _, diag = black_sea.monte_carlo_SPPR(
            n_samples=2, TE_option="GE", exclude_diverged=True,
            return_diagnostics=True, silent=True)

    assert diag["n_accepted"] == 0
    assert reject_frac == 1.0
    assert len(samples) == 0
    assert mean.isna().all().all(), "mean over zero accepted draws must be all-NaN"


# ------------------------------------------------------------------ symbolic path

def test_symbolic_kind_still_returns_its_equation_system(toy):
    """Regression guard: the symbolic branch was refactored into solve_draw alongside the rest."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        np.random.seed(0)
        mean, samples, reject_frac, eqs, vars_ = toy.monte_carlo_SPPR(
            n_samples=2, kind="symbolic", silent=True)

    assert eqs is not None and vars_ is not None, "symbolic must return its equations/variables"
    assert np.isfinite(mean.values).all()
    assert len(samples) == 2 and reject_frac == 0.0


def test_method_kwargs_reach_the_symbolic_solver_too(toy):
    """method_kwargs must be validated against SPPR_symbolic, not always against SPPR_new."""
    with pytest.raises(TypeError, match="nonexistent_param"):
        toy.monte_carlo_SPPR(n_samples=2, kind="symbolic",
                             method_kwargs={"nonexistent_param": 1}, silent=True)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        np.random.seed(0)
        mean, _, _, _, _ = toy.monte_carlo_SPPR(
            n_samples=2, kind="symbolic",
            method_kwargs={"sppr_det_value": None, "fix_EE_0_cases": True}, silent=True)

    assert np.isfinite(mean.values).all()


# ------------------------------------------------------------------ return shape

def test_default_call_still_returns_five_elements(toy):
    """Existing notebook call sites unpack exactly five values -- that must not break."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        np.random.seed(0)
        result = toy.monte_carlo_SPPR(n_samples=2, silent=True)

    assert isinstance(result, tuple)
    assert len(result) == 5


def test_return_diagnostics_appends_a_sixth_element(toy):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        np.random.seed(0)
        result = toy.monte_carlo_SPPR(n_samples=2, return_diagnostics=True, silent=True)

    assert len(result) == 6
    diag = result[-1]
    assert set(diag) == DIAGNOSTIC_KEYS


def test_diagnostics_counts_agree_with_the_returned_samples_and_fraction(toy):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        np.random.seed(0)
        _, samples, reject_frac, _, _, diag = toy.monte_carlo_SPPR(
            n_samples=4, return_diagnostics=True, silent=True)

    assert diag["n_accepted"] == len(samples)
    assert diag["reject_frac_negative"] + diag["reject_frac_diverged"] == pytest.approx(reject_frac)
    assert len(diag["per_sample_reason"]) == 4
    assert set(diag["per_sample_reason"]) <= {"accepted", "negative", "diverged"}


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
