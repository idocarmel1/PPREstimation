# `SPPR_2015` bug: dropped-group block assignment leaks SPPR=1 into a zero-production basal column

**Status:** open, pre-existing. Documented while verifying that `SPPR_new(TE, fix_EE_0_cases=False)`
reproduces `SPPR_2015` (it does, *except* where this bug fires). `SPPR_2015` is the
independent reference and must **not** be edited casually (per the user) — this note records
the bug and the fix for a deliberate future session.

**Scope:** affects only the SPPR of **EE=0 groups** (transfer-efficiency-0 dead-ends, `M0==p`)
in models where a **`PP_seq` member also has zero production** and is therefore dropped — in
practice an **Import pseudo-group with `p==0`**. The mis-assigned value is cosmetic: EE=0 groups
have zero export (catch+growth+net_migration), so their SPPR never enters PPR or the balance
identity. But it makes `SPPR_2015` internally **inconsistent** across models and impossible to
reproduce exactly.

---

## 1. Symptom

`SPPR_2015` assigns EE=0 groups (identical `ee==0`) **different** SPPR in different models:

| model | EE=0 group | `SPPR_2015` | Import (`PP_seq`) production |
|---|---|---|---|
| `47_10047_East_China_Sea_(1997)` | Marine mammals (23), Sharks (22) | **1.0** | Import (25) `p = -0.0` -> dropped |
| `47_11047_East_China_Sea_(2018)` | Marine mammals, Sharks | **1.0** | Import `p = 0` -> dropped |
| `48_12048_Southwestern_Yellow_Sea_(2008)` | Charadriiformes, Gruiformes, Anseriformes | **1.0** | Import `p = 0` -> dropped |
| `48_10048_North_Yellow_Sea_(2019)` | Hexagrammos otakii (4), Piscivores (2) | **0.0** | Import (19) `p = 1.39` -> kept |

The "clean" value (no zero-production Import) is **0.0**; the **1.0** is spurious.

---

## 2. Root cause

In `PPRCalculator.SPPR_2015`, after building the reduced production-requirement matrix
`L = (I - A)^-1` (with zero-production groups dropped), the dropped groups are re-added:

```python
seq_to_drop = P.index[P == 0]          # EE=0 groups AND any p==0 group (incl. a p==0 Import)
...
L[seq_to_drop] = 0                      # dropped columns -> 0
L = L.reindex(new_index).fillna(0)      # dropped rows    -> 0
L.loc[seq_to_drop, seq_to_drop] = 1     # <-- BUG: sets the ENTIRE dropped x dropped block to 1
SPPR = L.loc[PP_seq, :].T               # SPPR[g, src] = L[src, g]
```

`L.loc[seq_to_drop, seq_to_drop] = 1` is intended to put dropped groups on a unit diagonal
(self-basal), but it assigns the **whole submatrix** (every pair) to 1, so `L[d_i, d_j] = 1`
for *all* dropped `d_i, d_j`.

When a **`PP_seq` member is itself dropped** (a zero-production Import group, `p==0`), that
Import seq is in `seq_to_drop`. `SPPR = L.loc[PP_seq, :].T` then reads
`SPPR[ee0_group, Import] = L[Import, ee0_group] = 1` (off-diagonal block entry), so every other
dropped EE=0 group inherits a spurious **1.0** in the Import column.

Confirmed: in `47_10047`, `SPPR_2015.loc[23]` (Marine mammals) is `1.0` entirely in column `25`
(the dropped, zero-production Import group).

When the Import group has nonzero production it is not dropped, the off-diagonal leak lands in
no `PP_seq` column, and the EE=0 groups correctly read `0.0`.

---

## 3. Suggested fix

Set only the diagonal, not the whole block:

```python
for d in seq_to_drop:
    L.loc[d, d] = 1
# or: a = L.values; idx = [L.index.get_loc(d) for d in seq_to_drop]; a[idx, [L.columns.get_loc(d) for d in seq_to_drop]] = 1
```

Then EE=0 groups read `0.0` consistently regardless of whether an Import group has zero
production. (Whether dropped groups should be `0` or self-basal `1` is a separate convention
question; the bug is that the current code is *neither* consistently — it depends on Import `p`.)

Independently, a zero-production Import group probably should not be in `PP_seq`/basis at all;
worth deciding alongside the fix.

---

## 4. Reproduction

```python
import numpy as np, pandas as pd
from ModelData import ModelData
from PPRCalculator import PPRCalculator

def tot(s, self):
    s = PPRCalculator.rename_results(s, self.name2seq)
    return s.sum(axis=1) if isinstance(s, pd.DataFrame) else s

for f, ee0 in [('47_10047_East_China_Sea_(1997)', 23),      # Import p==0 -> EE0 gets spurious 1.0
               ('48_10048_North_Yellow_Sea_(2019)', 4)]:    # Import p!=0 -> EE0 gets clean 0.0
    m = ModelData(f'real_models/new_EwE_jsons/{f}.json')
    self = PPRCalculator.from_modeldata(m, underdetermined=True, zero_biomass_accum=False)
    s15 = PPRCalculator.rename_results(self.SPPR_2015()[0], self.name2seq)
    print(f, '-> EE0 group', ee0, 'SPPR_2015 =', float(s15.loc[ee0].sum()),
          '| Import p =', float(self.p[self.get_Import_seq()[0]]))
```

---

## 5. Relevant code

- `PPRCalculator.SPPR_2015` — the `L.loc[seq_to_drop, seq_to_drop] = 1` block assignment and the
  `SPPR = L.loc[PP_seq, :].T` read.

## 6. Related context

- Found while adding `SPPR_new(..., fix_EE_0_cases=True)`, which re-credits EE=0 dead-ends'
  consumed PP to detritus so `TE_option='TE'` balances (inflow == outflow) on single-DET models.
  With `fix_EE_0_cases=False`, `SPPR_new` reproduces *clean* `SPPR_2015` on all single-DET,
  non-EE=0 groups exactly; the only single-DET mismatches are this `SPPR_2015` EE=0 bug.
- Multi-DET TE mismatches are the separate `MULTIDET_TE_DISCREPANCY.md` issue.
- `SPPR_new` deliberately leaves EE=0 groups at SPPR=0 (clean-2015 behavior); it does **not**
  reproduce the spurious 1.0.
