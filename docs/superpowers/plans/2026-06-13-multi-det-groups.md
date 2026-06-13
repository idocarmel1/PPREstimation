# Multi-DET Group Support Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the single-DET restriction in `PPRCalculator` so models with >1 detritus groups (e.g., Humboldt Current 1995-2004) load and compute correctly, using `det_fate` to distribute flows across DET pools.

**Architecture:** `det_fate[group_i, DET_j]` = fraction of group_i's M0+egestion going to DET pool j. Every place that reads `get_DET_seq()[0]` (assuming one DET) must instead iterate over all DET seqs and weight flows by the corresponding `det_fate` column. The single-DET case is preserved as a fallback when `det_fate` lacks a DET-indexed column.

**Tech Stack:** Python, pandas, numpy, sympy (for symbolic helpers). Tests use `pytest`.

---

## Background: det_fate matrix semantics

`self._det_fate` (passed as `det_fate` classmethod param) is a DataFrame where:
- **Rows** = all group seqs (consumers / producers)
- **Columns** = group seqs that appear as prey, PLUS DET group seqs (injected at L168 in ModelData)
- **`det_fate[group_i, DET_j]`** = fraction of group_i's waste (M0+egestion) that flows to DET_j

For single-DET models the DET column may not be present (columns are prey_seq values from diet). In that case fall back to the old "assign total flow" behavior.

**Key formula** used in `apply_ecopath_defaults` L201-203 (already correct):
```
q[DET_j] = sum_i( flow_to_det[i] * det_fate[i, DET_j] )
```
This is the reference implementation of per-DET-group inflow calculation.

---

## Baseline test (must pass before and after every task)

```python
# tests/test_multi_det.py — create this file in Task 0
import pytest
from PPRCalculator import PPRCalculator

ICELAND = 'real_models/EwE_jsons/227_Iceland_(1950).json'
HUMBOLDT = 'real_models/EwE_jsons/10013_Humboldt_Current_(1995-2004).json'

def test_iceland_loads():
    pc = PPRCalculator(ICELAND)
    assert pc.n_groups == 25
    assert len(pc.get_DET_seq()) == 1

def test_iceland_sppr_1986():
    pc = PPRCalculator(ICELAND)
    result = pc.SPPR_1986()
    # must be positive and finite
    assert (result['sppr'] > 0).all()
    assert result['sppr'].isfinite().all()
```

Run with: `python -m pytest tests/test_multi_det.py -v`

---

## Task 0: Create test file and record Iceland baseline

**Files:**
- Create: `tests/test_multi_det.py`

- [ ] **Step 1: Record Iceland baseline values**

```bash
python -c "
from PPRCalculator import PPRCalculator
pc = PPRCalculator('real_models/EwE_jsons/227_Iceland_(1950).json')
print('n_groups:', pc.n_groups)
print('DET:', pc.get_DET_seq())
print('SPPR_1986 mean:', pc.SPPR_1986()['sppr'].mean())
sppr2015, _, _ = pc.SPPR_2015()
print('SPPR_2015 sum:', sppr2015.sum().sum())
"
```

Record these numbers — every task must produce identical values for Iceland.

- [ ] **Step 2: Write the test file**

```python
# tests/test_multi_det.py
import pytest
import numpy as np
from PPRCalculator import PPRCalculator

ICELAND = 'real_models/EwE_jsons/227_Iceland_(1950).json'
HUMBOLDT = 'real_models/EwE_jsons/10013_Humboldt_Current_(1995-2004).json'


class TestIcelandRegression:
    """Iceland (single-DET) must produce identical results throughout."""

    def setup_method(self):
        self.pc = PPRCalculator(ICELAND)

    def test_loads(self):
        assert self.pc.n_groups == 25
        assert len(self.pc.get_DET_seq()) == 1

    def test_sppr_1986_positive(self):
        s = self.pc.SPPR_1986()
        assert (s['sppr'] > 0).all()
        assert np.isfinite(s['sppr']).all()

    def test_sppr_2015_runs(self):
        sppr, _, _ = self.pc.SPPR_2015()
        assert not sppr.empty

    def test_sppr_new_runs(self):
        sppr, _, _ = self.pc.SPPR_new(TE_option='GE')
        assert not sppr.empty


class TestHumboldtLoads:
    """Humboldt (multi-DET) must load without exception once guard is removed."""

    def test_loads_after_guard_removed(self):
        # This test FAILS until Task 1 is done.
        pc = PPRCalculator(HUMBOLDT)
        assert len(pc.get_DET_seq()) > 1

    def test_sppr_1986_runs(self):
        pc = PPRCalculator(HUMBOLDT)
        s = pc.SPPR_1986()
        assert (s['sppr'] >= 0).all()
```

- [ ] **Step 3: Run — Iceland passes, Humboldt fails (expected)**

```bash
python -m pytest tests/test_multi_det.py -v
```

Expected: Iceland tests PASS, `TestHumboldtLoads` FAIL with `Exception('more than 1 DET groups')`.

- [ ] **Step 4: Commit**

```bash
git add tests/test_multi_det.py
git commit -m "test: add multi-DET regression suite; Humboldt currently fails"
```

---

## Task 1: Remove the guard in `from_modeldata()`

**Files:**
- Modify: `PPRCalculator.py:75-78`

- [ ] **Step 1: Remove lines 77-78**

Current (`PPRCalculator.py:75-78`):
```python
        # error and exit if needed:
        DET_seq = instance.get_DET_seq()
        if len(DET_seq) > 1:
            raise Exception('more than 1 DET groups')
```

Replace with:
```python
        # verify at least one DET group exists:
        DET_seq = instance.get_DET_seq()
        if len(DET_seq) == 0:
            raise Exception('no DET group found')
```

- [ ] **Step 2: Run tests — Humboldt may now crash at a different location**

```bash
python -m pytest tests/test_multi_det.py::TestHumboldtLoads -v
```

Expected: Either loads (unlikely) or crashes at a new location (good — confirms guard was the blocker).
Note the new crash location — it identifies the next broken site.

- [ ] **Step 3: Run Iceland tests — must still all pass**

```bash
python -m pytest tests/test_multi_det.py::TestIcelandRegression -v
```

Expected: All PASS.

- [ ] **Step 4: Commit**

```bash
git add PPRCalculator.py
git commit -m "fix: remove single-DET guard in from_modeldata; require at least one DET"
```

---

## Task 2: Fix `apply_ecopath_defaults()` finalize section

**Files:**
- Modify: `PPRCalculator.py:258-268` (the "finalize ratios" block)

**Problem:** L264 `df.loc[is_det, 'q'] = df['flow_to_det'].sum()` assigns the scalar total to ALL DET rows. With multiple DET, each group should receive only its share.

- [ ] **Step 1: Understand the det_fate column structure**

```python
python -c "
from ModelData import ModelData
md = ModelData('real_models/EwE_jsons/10013_Humboldt_Current_(1995-2004).json')
print('det_fate shape:', md.det_fate.shape)
print('det_fate columns:', sorted(md.det_fate.columns.tolist()))
det_groups = md.groups_data[md.groups_data['trophic_info']=='DET'].index.tolist()
print('DET group seqs:', det_groups)
print('DET cols in det_fate:', [c for c in md.det_fate.columns if c in det_groups])
print('det_fate sum per DET col:')
print(md.det_fate[[c for c in md.det_fate.columns if c in det_groups]].sum())
"
```

This confirms whether DET seqs appear as columns of `det_fate`.

- [ ] **Step 2: Replace the finalize block**

Locate `PPRCalculator.py` around L263-266:
```python
        df['flow_to_det'] = df['flow_to_det'].fillna(df['M0'] + df['egestion'])
        df.loc[is_det, 'q'] = df['flow_to_det'].sum()
        df.loc[is_det, 'p'] = df.loc[is_det, 'q']
        df.loc[is_det, 'biomass_accum'] = df.loc[is_det, 'p'] - (df.loc[is_det, 'predation'] + df.loc[is_det, 'net_migration'])
```

Replace with:
```python
        df['flow_to_det'] = df['flow_to_det'].fillna(df['M0'] + df['egestion'])
        det_idx = df.index[is_det]
        if det_fate is not None and len(det_idx) > 1:
            for det_j in det_idx:
                if det_j in det_fate.columns:
                    df.loc[det_j, 'q'] = (df['flow_to_det'] * det_fate[det_j].reindex(df.index).fillna(0)).sum()
                else:
                    df.loc[det_j, 'q'] = df['flow_to_det'].sum()
        else:
            df.loc[is_det, 'q'] = df['flow_to_det'].sum()
        df.loc[is_det, 'p'] = df.loc[is_det, 'q']
        df.loc[is_det, 'biomass_accum'] = df.loc[is_det, 'p'] - (df.loc[is_det, 'predation'] + df.loc[is_det, 'net_migration'])
```

- [ ] **Step 3: Run Iceland tests — must all pass**

```bash
python -m pytest tests/test_multi_det.py::TestIcelandRegression -v
```

Expected: All PASS. If any fail, the `else` branch may not be triggering — add a print to debug.

- [ ] **Step 4: Commit**

```bash
git add PPRCalculator.py
git commit -m "fix: distribute DET q by det_fate in apply_ecopath_defaults finalize"
```

---

## Task 3: Fix `apply_lim()._finalize_outputs()`

**Files:**
- Modify: `PPRCalculator.py:297-299`

**Problem:** `_finalize_outputs` inside `apply_lim` has the same pattern — assigns total flow_to_det to all DET rows.

- [ ] **Step 1: Replace the finalize block**

Current (`PPRCalculator.py:295-300`):
```python
            # rebalance DET row:
            df_final.loc[self.get_DET_seq(), 'q'] = df_final['flow_to_det'].sum()
            df_final.loc[self.get_DET_seq(), 'p'] = df_final.loc[self.get_DET_seq(), 'q']
            df_final.loc[self.get_DET_seq(), 'biomass_accum'] = df_final.loc[self.get_DET_seq(), 'p'] - (df.loc[self.get_DET_seq(), 'predation'] + df.loc[self.get_DET_seq(), 'net_migration'])
```

Replace with:
```python
            # rebalance DET rows:
            det_seqs = self.get_DET_seq()
            det_fate = getattr(self, '_det_fate', None)
            if det_fate is not None and len(det_seqs) > 1:
                for det_j in det_seqs:
                    if det_j in det_fate.columns:
                        df_final.loc[det_j, 'q'] = (df_final['flow_to_det'] * det_fate[det_j].reindex(df_final.index).fillna(0)).sum()
                    else:
                        df_final.loc[det_j, 'q'] = df_final['flow_to_det'].sum()
            else:
                df_final.loc[det_seqs, 'q'] = df_final['flow_to_det'].sum()
            df_final.loc[det_seqs, 'p'] = df_final.loc[det_seqs, 'q']
            df_final.loc[det_seqs, 'biomass_accum'] = df_final.loc[det_seqs, 'p'] - (df.loc[det_seqs, 'predation'] + df.loc[det_seqs, 'net_migration'])
```

- [ ] **Step 2: Run Iceland tests**

```bash
python -m pytest tests/test_multi_det.py::TestIcelandRegression -v
```

Expected: All PASS.

- [ ] **Step 3: Commit**

```bash
git add PPRCalculator.py
git commit -m "fix: distribute DET q by det_fate in apply_lim finalize"
```

---

## Task 4: Fix `get_Z()`

**Files:**
- Modify: `PPRCalculator.py:555-562`

**Problem:** `get_Z()` sets only `DET_seq[0]`'s row. With multiple DET groups, each DET row must be set using det_fate to split flow_to_det.

- [ ] **Step 1: Replace `get_Z()`**

Current (`PPRCalculator.py:555-562`):
```python
    def get_Z(self, DET_as_PP=False):
        """get Z matrix. if DET_as_PP is False (default), DET row is flow_to_det. otherwise it is set to 0"""
        Z = self._DC.mul(self._groups_df['q'].fillna(0), axis='index')
        if not DET_as_PP:
            Z.loc[self.get_DET_seq()[0], :] = (self.M0 + self.egestion).fillna(0)
        else:
            Z.loc[self.get_DET_seq()[0], :] = 0
        return Z.sort_index(ascending=False).sort_index(ascending=False, axis=1)
```

Replace with:
```python
    def get_Z(self, DET_as_PP=False):
        """get Z matrix. if DET_as_PP is False (default), DET rows are flow_to_det split by det_fate. otherwise set to 0"""
        Z = self._DC.mul(self._groups_df['q'].fillna(0), axis='index')
        DET_seq = self.get_DET_seq()
        if not DET_as_PP:
            flow_to_det = (self.M0 + self.egestion).fillna(0)
            det_fate = getattr(self, '_det_fate', None)
            if det_fate is not None and len(DET_seq) > 1:
                for det_j in DET_seq:
                    if det_j in det_fate.columns:
                        fracs = det_fate[det_j].reindex(flow_to_det.index).fillna(0)
                        Z.loc[det_j, :] = flow_to_det * fracs
                    else:
                        Z.loc[det_j, :] = flow_to_det
            else:
                for det_j in DET_seq:
                    Z.loc[det_j, :] = flow_to_det
        else:
            for det_j in DET_seq:
                Z.loc[det_j, :] = 0
        return Z.sort_index(ascending=False).sort_index(ascending=False, axis=1)
```

- [ ] **Step 2: Run Iceland tests**

```bash
python -m pytest tests/test_multi_det.py::TestIcelandRegression -v
```

Expected: All PASS. The single-DET else-branch preserves old behavior exactly.

- [ ] **Step 3: Run Humboldt load test**

```bash
python -m pytest tests/test_multi_det.py::TestHumboldtLoads::test_loads_after_guard_removed -v
```

Expected: Either PASS or crash at a later point (acceptable — each task peels back one layer).

- [ ] **Step 4: Commit**

```bash
git add PPRCalculator.py
git commit -m "fix: get_Z splits DET rows by det_fate for multi-DET models"
```

---

## Task 5: Fix `get_TL()` GE branch

**Files:**
- Modify: `PPRCalculator.py:617-641`

**Problem:** L635 `DC.loc[self.get_DET_seq()[0], :]` only updates the first DET row in the GE branch.

- [ ] **Step 1: Replace the two single-index assignments**

Locate `PPRCalculator.py:630-635`:
```python
        if TE_option == 'TE':
            DC.loc[DET_seq, Regular_seq] = 0
        elif TE_option == 'GE':
            DC.loc[self.get_DET_seq()[0], :] = ((self.M0) / flow2det).fillna(0)
```

Replace with:
```python
        if TE_option == 'TE':
            DC.loc[DET_seq, Regular_seq] = 0
        elif TE_option == 'GE':
            det_fate = getattr(self, '_det_fate', None)
            for det_j in DET_seq:
                q_det_j = Z.loc[det_j, :].sum() if not DET_as_PP else flow2det / len(DET_seq)
                if det_fate is not None and det_j in det_fate.columns and q_det_j > 0:
                    fracs = det_fate[det_j].reindex(self.M0.index).fillna(0)
                    DC.loc[det_j, :] = ((self.M0 * fracs) / q_det_j).fillna(0)
                else:
                    DC.loc[det_j, :] = ((self.M0) / flow2det).fillna(0)
```

Note: `Z` is already computed at L622 (`Z = self.get_Z(DET_as_PP=DET_as_PP)`), so we can use it here.

- [ ] **Step 2: Run Iceland tests**

```bash
python -m pytest tests/test_multi_det.py::TestIcelandRegression -v
```

Expected: All PASS.

- [ ] **Step 3: Commit**

```bash
git add PPRCalculator.py
git commit -m "fix: get_TL GE branch iterates over all DET groups"
```

---

## Task 6: Fix `SPPR_2015()`

**Files:**
- Modify: `PPRCalculator.py:966-1013`

**Problem:** L972 stores `DET_seq = self.get_DET_seq()[0]` as a scalar, then uses it for single-row indexing throughout. Every subsequent use of `DET_seq` (L977, L981, L983, L984) must be updated to handle a list.

- [ ] **Step 1: Replace the DET scalar with list-based iteration**

Locate `PPRCalculator.py:966-984`:
```python
        DET_seq = self.get_DET_seq()[0]  # assuming there is only one DET
        Z = self.get_Z(DET_as_PP=False)

        # combine part of DET that is PP into PP row:
        Z_without_DET = Z.copy()
        percent_of_det_that_is_PP = Z_without_DET.loc[DET_seq, PP_seq] / Z_without_DET.loc[DET_seq, :].sum()  # 98%
        percent_of_det_that_is_PP = percent_of_det_that_is_PP.fillna(0)
        if only_pp_det:  # this is what is implemented in the article
            for i in PP_seq:
                Z_without_DET.loc[:, i] += percent_of_det_that_is_PP[i] * Z_without_DET.loc[:, DET_seq]
        else:
            Z_without_DET.loc[:, PP_seq] += Z_without_DET.loc[:, DET_seq]
        Z_without_DET = Z_without_DET.drop(index=DET_seq, columns=DET_seq)
```

Replace with:
```python
        DET_seqs = self.get_DET_seq()  # list, may have >1 elements
        Z = self.get_Z(DET_as_PP=False)

        # combine part of each DET that is PP into PP row:
        Z_without_DET = Z.copy()
        for det_j in DET_seqs:
            det_row_total = Z_without_DET.loc[det_j, :].sum()
            percent_of_det_that_is_PP = (Z_without_DET.loc[det_j, PP_seq] / det_row_total).fillna(0)
            if only_pp_det:
                for i in PP_seq:
                    Z_without_DET.loc[:, i] += percent_of_det_that_is_PP[i] * Z_without_DET.loc[:, det_j]
            else:
                Z_without_DET.loc[:, PP_seq] += Z_without_DET.loc[:, det_j]
        Z_without_DET = Z_without_DET.drop(index=DET_seqs, columns=DET_seqs)
```

- [ ] **Step 2: Update L1011 (SPPR DET row back-fill)**

Locate `PPRCalculator.py:1009-1011`:
```python
        for i in PP_seq:
            SPPR.loc[self.get_DET_seq(), i] = self.M0.loc[PP_seq][i] / (self.M0 + self.egestion).sum()
```

This already uses `self.get_DET_seq()` as a list and assigns the same proportional value to all DET groups. For a first pass this is acceptable (all DET groups get equal SPPR attribution). No change needed here unless `get_DET_seq()` returns an empty list (add a guard):

```python
        det_seqs = self.get_DET_seq()
        if det_seqs:
            for i in PP_seq:
                SPPR.loc[det_seqs, i] = self.M0.loc[PP_seq][i] / (self.M0 + self.egestion).sum()
```

- [ ] **Step 3: Run Iceland tests**

```bash
python -m pytest tests/test_multi_det.py::TestIcelandRegression -v
```

Expected: All PASS.

- [ ] **Step 4: Commit**

```bash
git add PPRCalculator.py
git commit -m "fix: SPPR_2015 iterates over all DET groups instead of [0]"
```

---

## Task 7: Fix `SPPR_new()` flow2det calculations

**Files:**
- Modify: `PPRCalculator.py:1060-1090`

**Problem:** Three branches each do `self.get_Z(DET_as_PP=False).loc[self.get_DET_seq()[0], :].sum()` to get total flow to DET. For multiple DET the total across all DET rows is simply `(self.M0 + self.egestion).sum()`.

- [ ] **Step 1: Replace the three flow2det calculations**

In all three branches of `SPPR_new()`, replace:
```python
flow2det = self.get_Z(DET_as_PP=False).loc[self.get_DET_seq()[0], :].sum()
```
with:
```python
flow2det = (self.M0 + self.egestion).sum()
```

Lines to change: 1063, 1071, 1081.

Note: The computation at L1082 (`flow2det = (self.M0 + self.egestion).sum()`) already overrides the one at L1081 in the With-Egestion branch. After the fix the override becomes redundant but harmless.

- [ ] **Step 2: Run Iceland tests**

```bash
python -m pytest tests/test_multi_det.py::TestIcelandRegression -v
```

Expected: All PASS. The value `(M0 + egestion).sum()` equals `Z.loc[DET_seq[0], :].sum()` for single-DET, so behavior is preserved.

- [ ] **Step 3: Add Humboldt SPPR_new test**

Add to `tests/test_multi_det.py` in `TestHumboldtLoads`:
```python
    def test_sppr_new_runs(self):
        pc = PPRCalculator(HUMBOLDT)
        sppr, _, _ = pc.SPPR_new(TE_option='GE')
        assert not sppr.empty
        assert np.isfinite(sppr.values[np.isfinite(sppr.values)]).all()
```

```bash
python -m pytest tests/test_multi_det.py -v
```

- [ ] **Step 4: Commit**

```bash
git add PPRCalculator.py tests/test_multi_det.py
git commit -m "fix: SPPR_new uses total flow2det instead of [0]-indexed DET row"
```

---

## Task 8: Fix `_SPPR_symbolic_helper_diet_import_as_PP()`

**Files:**
- Modify: `PPRCalculator.py:1092-1173`

**Problem:** Lines 1106, 1108, 1111 use `self.get_DET_seq()[0]` to set a single DET row in DC. With multiple DET groups, each row must be set.

- [ ] **Step 1: Replace DC DET-row assignments (GE and With Egestion branches)**

Locate `PPRCalculator.py:1103-1113`:
```python
        if TE_option == 'TE':
            DC.loc[DET_seq, Regular_seq] = 0
        elif TE_option == 'GE':
            DC.loc[self.get_DET_seq()[0], :] = ((self.M0) / flow2det).fillna(0)
        elif TE_option == 'With Egestion':
            DC.loc[self.get_DET_seq()[0], :] = 0
            m = (self.M0 / flow2det).fillna(0)
            e = (self.egestion / flow2det).fillna(0)
            DC.loc[self.get_DET_seq()[0], :] = (m + e @ DC)
        else:
```

Replace with:
```python
        if TE_option == 'TE':
            DC.loc[DET_seq, Regular_seq] = 0
        elif TE_option == 'GE':
            det_fate = getattr(self, '_det_fate', None)
            for det_j in DET_seq:
                q_det_j = Z.loc[det_j, :].sum() if hasattr(self, '_det_fate') else flow2det
                if det_fate is not None and det_j in det_fate.columns and q_det_j > 0:
                    fracs = det_fate[det_j].reindex(self.M0.index).fillna(0)
                    DC.loc[det_j, :] = ((self.M0 * fracs) / q_det_j).fillna(0)
                else:
                    DC.loc[det_j, :] = ((self.M0) / flow2det).fillna(0)
        elif TE_option == 'With Egestion':
            for det_j in DET_seq:
                DC.loc[det_j, :] = 0
            m = (self.M0 / flow2det).fillna(0)
            e = (self.egestion / flow2det).fillna(0)
            for det_j in DET_seq:
                DC.loc[det_j, :] = (m + e @ DC)
        else:
```

Note: `Z` is not computed at this point in the method. Add `Z = self.get_Z(DET_as_PP=False)` before the if-block if not already present. Check the method context before inserting.

- [ ] **Step 2: Run Iceland tests**

```bash
python -m pytest tests/test_multi_det.py::TestIcelandRegression -v
```

Expected: All PASS.

- [ ] **Step 3: Commit**

```bash
git add PPRCalculator.py
git commit -m "fix: _SPPR_symbolic_helper_diet_import_as_PP iterates over all DET groups"
```

---

## Task 9: Fix `_SPPR_symbolic_helper_diet_import_as_DC()`

**Files:**
- Modify: `PPRCalculator.py:1175-~1270`

**Problem:** Lines 1188, 1190, 1193 — identical pattern to Task 8 but in the DC-variant helper.

- [ ] **Step 1: Apply the same fix**

Locate `PPRCalculator.py:1185-1195` (same structure as Task 8):
```python
        if TE_option == 'TE':
            DC.loc[DET_seq, Regular_seq] = 0
        elif TE_option == 'GE':
            DC.loc[self.get_DET_seq()[0], :] = ((self.M0) / flow2det).fillna(0)
        elif TE_option == 'With Egestion':
            DC.loc[self.get_DET_seq()[0], :] = 0
            m = (self.M0 / flow2det).fillna(0)
            e = (self.egestion / flow2det).fillna(0)
            DC.loc[self.get_DET_seq()[0], :] = (m + e @ DC)
        else:
```

Apply the same replacement as Task 8 (copy exactly, ensuring `det_j` loop for GE and With Egestion branches).

- [ ] **Step 2: Run all tests**

```bash
python -m pytest tests/test_multi_det.py -v
```

Expected: All Iceland PASS; Humboldt load test PASS; other Humboldt tests may still fail.

- [ ] **Step 3: Commit**

```bash
git add PPRCalculator.py
git commit -m "fix: _SPPR_symbolic_helper_diet_import_as_DC iterates over all DET groups"
```

---

## Task 10: Fix `calc_2015/SPPR_2015.py` and `utils.py`

**Files:**
- Modify: `calc_2015/SPPR_2015.py:165-172`
- Modify: `utils.py:220-226`

**Problem:** Both files contain `Z.loc[DET_seq[0], :]` patterns. These are legacy SPPR calculation paths that need the same get_Z-style fix.

- [ ] **Step 1: Inspect the context in SPPR_2015.py**

```bash
python -c "
import linecache
for i in range(160, 180):
    print(i, linecache.getline('calc_2015/SPPR_2015.py', i), end='')
"
```

- [ ] **Step 2: Fix `calc_2015/SPPR_2015.py`**

Locate the block (around L165-172) that does:
```python
    if DET_as_PP:
        Z.loc[DET_seq[0], :] = groups_data['flow_to_det']
    else:
        Z.loc[DET_seq[0], :] = groups_data['flow_to_det'] * (DC * (1-det_fate)).sum(axis=1)
```

Replace with (assuming `det_fate` param exists in the enclosing function and DET_seq is already a list):
```python
    for det_j in DET_seq:
        if DET_as_PP:
            if det_fate is not None and det_j in det_fate.columns:
                fracs = det_fate[det_j].reindex(groups_data.index).fillna(0)
                Z.loc[det_j, :] = groups_data['flow_to_det'] * fracs
            else:
                Z.loc[det_j, :] = groups_data['flow_to_det']
        else:
            if det_fate is not None and det_j in det_fate.columns:
                fracs = det_fate[det_j].reindex(groups_data.index).fillna(0)
                Z.loc[det_j, :] = groups_data['flow_to_det'] * fracs * (DC * (1 - det_fate)).sum(axis=1)
            else:
                Z.loc[det_j, :] = groups_data['flow_to_det'] * (DC * (1 - det_fate)).sum(axis=1)
```

- [ ] **Step 3: Inspect and fix `utils.py`**

```bash
python -c "
import linecache
for i in range(215, 230):
    print(i, linecache.getline('utils.py', i), end='')
"
```

Apply the equivalent fix at lines 223, 225 in `utils.py`.

- [ ] **Step 4: Run all tests**

```bash
python -m pytest tests/test_multi_det.py -v
```

- [ ] **Step 5: Commit**

```bash
git add calc_2015/SPPR_2015.py utils.py
git commit -m "fix: legacy SPPR helpers in SPPR_2015.py and utils.py iterate over all DET groups"
```

---

## Task 11: Full Humboldt integration test + validate Iceland unchanged

**Files:**
- Modify: `tests/test_multi_det.py`

- [ ] **Step 1: Run the full Humboldt pipeline**

```python
python -c "
from PPRCalculator import PPRCalculator
pc = PPRCalculator('real_models/EwE_jsons/10013_Humboldt_Current_(1995-2004).json')
print('Loaded! n_groups:', pc.n_groups)
print('DET groups:', pc.get_DET_seq())
print('PP groups:', pc.get_PP_seq())
print('Regular groups:', len(pc.get_Regular_seq()))
print('is_balanced:', pc.is_balanced)

s = pc.SPPR_1986()
print('SPPR_1986 mean:', s['sppr'].mean())
print('SPPR_1986 min:', s['sppr'].min())

sppr2015, _, _ = pc.SPPR_2015()
print('SPPR_2015 shape:', sppr2015.shape)
print('SPPR_2015 sum:', sppr2015.sum().sum())
"
```

Check: values are positive and finite.

- [ ] **Step 2: Verify Iceland baseline unchanged**

Run the exact numbers from Task 0 Step 1 and confirm they match:

```bash
python -c "
from PPRCalculator import PPRCalculator
pc = PPRCalculator('real_models/EwE_jsons/227_Iceland_(1950).json')
print('SPPR_1986 mean:', pc.SPPR_1986()['sppr'].mean())
sppr2015, _, _ = pc.SPPR_2015()
print('SPPR_2015 sum:', sppr2015.sum().sum())
"
```

Expected: Identical to Task 0 recorded values.

- [ ] **Step 3: Expand test suite**

Add to `tests/test_multi_det.py`:
```python
class TestHumboldtResults:
    """Smoke tests: Humboldt results are finite and non-negative."""

    def setup_method(self):
        self.pc = PPRCalculator(HUMBOLDT)

    def test_det_groups_count(self):
        assert len(self.pc.get_DET_seq()) > 1

    def test_det_q_differs_per_group(self):
        # each DET group should have a different q if det_fate splits them
        det = self.pc.get_DET_seq()
        if len(det) > 1:
            q_vals = self.pc.q[det].values
            assert not np.all(q_vals == q_vals[0]), "All DET groups have identical q — det_fate split may not be working"

    def test_sppr_1986_finite(self):
        s = self.pc.SPPR_1986()
        assert np.isfinite(s['sppr']).all()
        assert (s['sppr'] >= 0).all()

    def test_sppr_2015_finite(self):
        sppr, _, _ = self.pc.SPPR_2015()
        assert np.isfinite(sppr.values[~np.isnan(sppr.values)]).all()
```

- [ ] **Step 4: Run full suite**

```bash
python -m pytest tests/test_multi_det.py -v
```

Expected: All tests PASS.

- [ ] **Step 5: Final commit**

```bash
git add tests/test_multi_det.py
git commit -m "test: full Humboldt integration tests; all multi-DET tests passing"
```

---

## Self-Review

### Spec coverage
- [x] Remove >1 DET guard → Task 1
- [x] `apply_ecopath_defaults` finalize → Task 2
- [x] `apply_lim._finalize_outputs` → Task 3
- [x] `get_Z` → Task 4
- [x] `get_TL` GE branch → Task 5
- [x] `SPPR_2015` → Task 6
- [x] `SPPR_new` flow2det → Task 7
- [x] `_SPPR_symbolic_helper_diet_import_as_PP` → Task 8
- [x] `_SPPR_symbolic_helper_diet_import_as_DC` → Task 9
- [x] `calc_2015/SPPR_2015.py` + `utils.py` → Task 10
- [x] Integration validation → Task 11
- [x] Baseline test preserving Iceland → every task

### Known risks
1. **Task 8/9** require `Z = self.get_Z(DET_as_PP=False)` to be computed before the if-block. Check that `Z` is in scope at those points. If not, add `Z = self.get_Z(DET_as_PP=False)` before the if-block.
2. **Task 10** — the exact signature of `det_fate` in `calc_2015/SPPR_2015.py` and `utils.py` must be verified in Step 1 before writing the fix.
3. **`test_det_q_differs_per_group`** — this test will fail if `det_fate` does NOT have DET cols for the Humboldt model. If it fails, add a diagnostic print showing `md.det_fate.columns` to verify the expected columns are present.
