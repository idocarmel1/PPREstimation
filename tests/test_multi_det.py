import pytest
import numpy as np
from PPRCalculator import PPRCalculator

ICELAND = 'real_models/EwE_jsons/227_Iceland_(1950).json'
HUMBOLDT = 'real_models/EwE_jsons/10013_Humboldt_Current_(1995-2004).json'


class TestIcelandRegression:
    """Iceland (single-DET) must produce identical results throughout all code changes.

    Baseline recorded 2026-06-13:
      n_groups      : 25
      DET seq       : [np.int64(24)]
      SPPR_1986 mean: 224.5091368959694
      SPPR_2015 sum : 407203.9410639084
    """

    def setup_method(self):
        self.pc = PPRCalculator(ICELAND)

    def test_loads(self):
        assert self.pc.n_groups == 25
        assert len(self.pc.get_DET_seq()) == 1

    def test_sppr_1986_positive(self):
        s = self.pc.SPPR_1986()
        assert (s['sppr'] > 0).all()
        assert np.isfinite(s['sppr']).all()

    def test_sppr_1986_mean_regression(self):
        s = self.pc.SPPR_1986()
        assert abs(s['sppr'].mean() - 224.5091368959694) < 1e-6

    def test_sppr_2015_runs(self):
        sppr, _, _ = self.pc.SPPR_2015()
        assert not sppr.empty

    def test_sppr_2015_sum_regression(self):
        sppr, _, _ = self.pc.SPPR_2015()
        assert abs(sppr.sum().sum() - 407203.9410639084) < 1e-3

    @pytest.mark.xfail(
        reason="SPPR_new has a pre-existing read-only array bug (ValueError: assignment "
               "destination is read-only at PPRCalculator.py:1027). "
               "Will be fixed in Task 8.",
        raises=ValueError,
        strict=True,
    )
    def test_sppr_new_runs(self):
        sppr, _, _ = self.pc.SPPR_new(TE_option='GE')
        assert not sppr.empty


class TestHumboldtLoads:
    """Humboldt (multi-DET) must load without exception once guard is removed.
    These tests FAIL until Task 1 removes the guard in PPRCalculator."""

    def test_loads_after_guard_removed(self):
        pc = PPRCalculator(HUMBOLDT)
        assert len(pc.get_DET_seq()) > 1

    def test_sppr_1986_runs(self):
        pc = PPRCalculator(HUMBOLDT)
        s = pc.SPPR_1986()
        assert (s['sppr'] >= 0).all()
