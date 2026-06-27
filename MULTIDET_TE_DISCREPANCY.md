# Multi-DET TE discrepancy: `SPPR_new` vs `SPPR_symbolic(as_PP)`

**Status:** open, pre-existing (NOT introduced by the open-system / `det_fate` work of
2026-06-26). Noted for a future session.

**Scope:** affects only `TE_option='TE'` on **multi-detritus models that contain secondary
detritus pools** — pools fed by regular (consumer) groups rather than directly by
primary producers / import. Single-DET models and multi-DET models whose pools are all
fed directly by PP+Import are unaffected.

---

## 1. Symptom

`SPPR_symbolic(diet_import_option='as_PP')` is meant to reproduce `SPPR_new` for matching
`TE_option`. After the 2026-06-26 open-system fixes it does so **exactly** for:

| case | TE | GE |
|---|---|---|
| single-DET, closed (e.g. `36_11036`, `227_Iceland`) | ✅ | ✅ |
| single-DET, export (e.g. `105_Grand_Banks`, `650_Mauritania`) | ✅ | ✅ |
| multi-DET, all pools PP-fed (e.g. `705_Georges_Bank` nDET=2, toy `900`) | ✅ | ✅ |

But it diverges on **`13_11013_Humboldt_Current_(1995-2004)`** (nDET=4):

- **TE**: relative error ≈ **0.71** (huge) — FAIL
- **GE**: relative error ≈ **0.001** (matches within tolerance)

So the problem is specific to the **TE** method on this model class. GE already agrees.

---

## 2. Root cause

`13_11013` has 4 detritus pools:

| seq | name | q[det] | fed by |
|---|---|---|---|
| 36 | Anchovy eggs | 82.1 | anchovy (regular) |
| 37 | Fishery offal | 0.0 | fishery/discards |
| 38 | Pelagic detritus | 5081.4 | PP + many groups |
| 39 | Benthic detritus | 1749.0 | benthic groups (regular) |

`SPPR_new`'s **TE branch** scales each detritus column by its **direct PP+Import inflow
share only**:

```python
# PPRCalculator.SPPR_new, TE_option == 'TE' branch
PP_Import_seq = list(self.get_PP_seq() + self.get_Import_seq())
for i, det_l in enumerate(DET_seq):
    q_l = self.q[det_l] if self.q[det_l] > 0 else flow_to_det.sum()
    if det_fate is not None and det_l in det_fate.columns:
        fr = det_fate[det_l].reindex(flow_to_det.index).fillna(0)
        inflow = (flow_to_det[PP_Import_seq] * fr[PP_Import_seq]).sum()   # <-- PP+Import ONLY
    ...
    SPPR[det_l] *= (inflow / q_l) * theta[i]
```

For a **secondary pool** (anchovy eggs, offal, benthic detritus) **no PP or Import group
feeds it directly**, so `inflow = 0` and the whole column is zeroed. Observed DET column
sums for `13_11013` TE: `[0, 0, 0.31, 0]` — only "Pelagic detritus" (PP-fed) survives.

Consequence: living groups that eat the zeroed pools lose that PP-requirement pathway, so
their `SPPR_new` values come out **~1.7× lower** than `SPPR_symbolic(as_PP)`, which solves
the detritus pools inside the food-web linear system and therefore traces the indirect PP
that reaches a secondary pool via the consumers feeding it.

**Why single-DET matches:** with one pool, *all* `flow_to_det` is that pool's, and the
recycling of consumer-derived detritus is captured by the nullspace basis; the direct
PP+Import share is the correct scalar. The mismatch only appears once detritus is split
into pools at different "distances" from PP.

**Why GE matches:** the GE path (`_build_det_BC`) builds the recycling matrix from **all
groups'** `M0 * det_fate / q[det]` (not just PP+Import) and solves the coupled
`(I - B) x = c` system — i.e. it already traces the indirect inflow, like the symbolic
solver. So GE and `as_PP`-GE agree.

---

## 3. The design question (needs a decision before coding)

Which detritus-pool SPPR is *correct* for the TE method?

- **(a) Direct-PP-inflow view** (current `SPPR_new` TE): a pool's SPPR = fraction of its
  inflow arriving **directly** from PP+Import. A pool fed only by consumers → SPPR 0.
- **(b) Food-web-trace view** (current symbolic + GE): a pool's SPPR traces **all** the PP
  that reaches it, including PP that passed through consumers first. Secondary pools get
  nonzero SPPR.

(b) is the more physically defensible (benthic detritus *does* embody PP — the benthos ate
PP), and it is what GE and the 2015 method already do. (a) is the literal TE-method scalar
that happens to coincide with (b) only for single-DET.

**Decision needed:** keep (a) as the intended TE semantics and make `as_PP` reproduce it
(by zeroing secondary pools in the symbolic TE path too), OR change `SPPR_new`'s TE branch
to (b) so it traces secondary pools (then `as_PP` already matches).

---

## 4. Suggested approach (if we pursue (b) — recommended)

Make `SPPR_new`'s TE detritus scaling trace indirect inflow instead of using only the
direct PP+Import share. Concretely, replace the per-pool `inflow = PP+Import direct` scalar
with the same **coupled recycling solve the GE path already uses**, but with the TE
definition of the per-group detritus contribution:

1. In `SPPR_new`, route the `TE_option='TE'` detritus step through `_build_det_BC` +
   `_solve_det_scaling` (as GE does), using `flow_to_det = M0 + egestion` for the TE
   contribution (GE uses `M0` only; "With Egestion" uses `M0 + DC@egestion`).
2. This makes secondary pools inherit SPPR from the consumers feeding them, matching the
   symbolic `as_PP` TE path and the GE results.
3. Re-verify single-DET still matches exactly (it must reduce to the current scalar when
   there is one pool), and that `is_sppr_balanced` is unchanged for single-DET.

If instead we keep (a): change the symbolic `as_PP` TE path to zero any detritus pool with
zero direct PP+Import inflow, so it reproduces `SPPR_new`. (Less physically meaningful;
only do this if (a) is deliberately the intended TE definition.)

**Do NOT** "fix" this by special-casing in only one method — the whole point of the
2026-06-26 work was to keep `SPPR_new`, the symbolic helpers, `get_Z`, and `apply_lim`
mutually consistent. Whatever semantics is chosen must hold in both `SPPR_new` and the
symbolic `as_PP` helper.

---

## 5. Reproduction

```python
import numpy as np, pandas as pd
from ModelData import ModelData
from PPRCalculator import PPRCalculator

m = ModelData('real_models/new_EwE_jsons/13_11013_Humboldt_Current_(1995-2004).json')
self = PPRCalculator.from_modeldata(m, underdetermined=True, zero_biomass_accum=False)

def total(s):
    s = PPRCalculator.rename_results(s, self.name2seq)
    return s.sum(axis=1) if isinstance(s, pd.DataFrame) else s

for opt in ['TE', 'GE']:
    sn = total(self.SPPR_new(TE=None, TE_option=opt, collapse_det=False)[0])
    ss = total(self.SPPR_symbolic(TE=None, TE_option=opt, diet_import_option='as_PP')[1])
    idx = sn.index.union(ss.index); sn = sn.reindex(idx).fillna(0); ss = ss.reindex(idx).fillna(0)
    rel = float((sn - ss).abs().sum() / (sn.abs().sum() + 1e-12))
    print(opt, 'relErr=', round(rel, 4))

# Inspect the zeroed pools:
sppr_new = self.SPPR_new(TE=None, TE_option='TE', collapse_det=False)[0]
print('DET column sums:', sppr_new[self.get_DET_seq()].sum().values)  # ~ [0, 0, 0.31, 0]
```

A clean multi-DET control that **does** match (use to confirm a fix preserves correct cases):
`real_models/EwE_jsons/705_705_Georges_Bank_(1996).json` (nDET=2) and the
`real_models/ToyModels/900_900_Multi_DET_Toy_(2026).json`.

---

## 6. Relevant code (function names, not line numbers — lines drift)

- `PPRCalculator.SPPR_new` — the `TE_option == 'TE'` detritus branch (the direct-PP-inflow
  scalar that zeroes secondary pools). This is where (b) would be implemented.
- `PPRCalculator._build_det_BC` / `_solve_det_scaling` — the coupled recycling solver the
  GE path uses; reuse for the TE path under approach (b).
- `PPRCalculator._SPPR_symbolic_helper_diet_import_as_PP` — the symbolic `as_PP` helper;
  its TE branch builds the detritus diet row from `get_DC(DET_as_PP=False)` (now det_fate-
  consistent via `get_Z`).
- `PPRCalculator.get_Z(DET_as_PP=False)` — detritus row = `flow_to_det * det_fate` per pool
  (fixed 2026-06-26 to apply for single-DET too).

## 7. Related context

- The 2026-06-26 open-system work (commit "Open-system detritus via det_fate…") fixed the
  single-DET `det_fate` regression and added `ModelData.validate_DC` / `validate_det_fate`
  and per-group `det_export`. It deliberately left `is_sppr_balanced` on the classic
  outflow — a `det_export` term overshoots and the residual balance gaps are a
  method-independent property of the balance identity (confirmed identical under
  `SPPR_2015`), not an `SPPR_new` defect.
- `SPPR_2015` was used as an independent reference and must **not** be edited (per the
  user); it agrees with `SPPR_new` to ~0.1% across models, so it is a good control here too.
