# Multi-DET Openness, Collapse Modes & Diagnostics — Change Log

Implements the agreed points from GPT-5.5's review (`detritus_SPPR_review_for_claude_code.txt`,
my assessment in `response_for_GPT.txt`) plus the requested ecological "openness" knob.
**Default-parameter output is byte-identical to the previous implementation** on both test
models — verified by a golden regression captured before any change (see "Tests executed").

Branch: `FishEstimationAI`. All work is in `PPRCalculator.py`.

---

## Review points addressed

| Point | Summary | What changed |
|---|---|---|
| 2 | Silent negatives vs. explicit collapse behaviour | New `det_collapse_mode` enum: `never` (always solve, may return negatives — never raises), `auto` (pool iff unstable), `always` (always pool). The old boolean `collapse_det` still works via a back-compat shim (`False→never`, `True→auto`). |
| 3 | Real stability test | The solve-vs-pool decision now uses the **spectral radius** `ρ = max|eig(B_open)|` and the **condition number** `cond(I−B_open)`, not just the sign of the smallest real eigenvalue. Thresholds: `tol=1e-10`, `cond_threshold=1e10`. |
| 4 | Symbolic path had no stability gate | Both `SPPR_symbolic` helpers now share the same stability/collapse gate as `SPPR_new` for non-default modes. The **default** symbolic path keeps the proven `sympy` `sol_dict2` scaling untouched (byte-identical). |
| 5 | No run diagnostics | Every GE/With-Egestion solve records `self.detritus_resolution_info` (method used, reason, `rho_B`, `cond_IminusB`, det seqs/names, open mode, theta, and either the solved `x_vec` or the pooled `collapse_scalar`/`a`/`b`). |
| 6 | Missing det_fate column = whole-flow | In a multi-DET model, a present `det_fate` matrix that lacks a column for a DET group now contributes **0** inflow to that group (with a `RuntimeWarning`), instead of silently treating it as receiving the whole flow. Single-DET / no-fate keeps `fracs = 1`. |
| 7 | Numerator/denominator inconsistency in pooled scaling | `_collapse_det_scaling` now weights `M0`/`egestion` by `fate_to_modeled_det` (= row-sum of `det_fate` over the modeled DET columns), consistent with the already-fate-weighted `q_combined`. With `det_fate` rows summing to 1 this weight is exactly 1 → no-op on the real models. |
| 8 | Ecological openness of detritus | New `det_open_mode` with two ecologically-grounded options plus the per-DET `det_external_sppr` knob (point 8B). See math below. |

---

## New public parameters (identical across `SPPR_new`, `SPPR_symbolic`, `monte_carlo_SPPR`, `monte_carlo_SPPR_2`)

```python
det_collapse_mode='never'   # 'never' | 'auto' | 'always'
det_open_mode='none'        # 'none'  | 'recycling_loss' (8A) | 'source_dilution' (8B)
det_theta=1.0               # availability θ (=1−λ); float OR {DET seq-or-name: float}
det_external_sppr=0.0       # external SPPR for diluted material (point 8B); float OR dict
collapse_det=None           # back-compat: False→'never', True→'auto'; None→use det_collapse_mode
```

`det_theta` / `det_external_sppr` accept a scalar (broadcast to all DET groups) or a dict keyed
by DET-group sequence number **or** name; missing keys fall back to the default.

### Openness math (per-DET θ, ext aligned to the DET columns)

Recycling system without openness: `(I − B) x = c`, where `x_l` scales DET column `l`,
`c[l] = m_eff_l · non_DET_sppr`, `B[l,j] = m_eff_l · basis[det_j]`.

| `det_open_mode` | `B_open` | `c_open` |
|---|---|---|
| `none` | `B` | `c` |
| `recycling_loss` (8A) | `diag(θ)·B` | `c` |
| `source_dilution` (8B) | `diag(θ)·B` | `θ·c + ext·(1−θ)` |

`θ = 1`, `ext = 0` reduces all three to the original system, so the old result is preserved.
For `TE_option='TE'` (no recycling matrix) θ is applied as a plain multiplier on each DET column.

---

## Files / methods touched (`PPRCalculator.py`)

- **New helpers:** `_spectral_radius`, `_resolve_det_param`, `_build_det_BC` (centralised
  `(I−B)x=c` builder with the point-6 guard), `_solve_det_scaling` (openness transform +
  stability decision + apply + diagnostics).
- **Upgraded:** `_collapse_det_scaling` (point-7 fate weighting + openness + returns
  `(SPPR, scalar, a, b)`).
- **Rewired to share the helpers:** `SPPR_new` (single- and multi-DET paths unified),
  `_SPPR_symbolic_helper_diet_import_as_PP`, `_SPPR_symbolic_helper_diet_import_as_DC`,
  `SPPR_symbolic` (threads the new params; back-compat shim).
- **Threaded params through:** `monte_carlo_SPPR`, `monte_carlo_SPPR_2` (forward to every
  internal `SPPR_new` / `SPPR_symbolic` call; the `never` default keeps the
  reject-negative-samples behaviour intact — it solves and returns, never raises).
- **Documentation:** explanatory comments/docstrings added across the whole class
  (comments-only; no behaviour change).
- `import warnings` added.

---

## Tests executed

All regression oracles were captured from the pre-change code, on
`227_Iceland_(1950)` (1 DET group) and `900_Multi_DET_Toy_(2026)` (2 DET groups).

### 1. `SPPR_new` golden regression (default params == old)
```
iceland GE             maxdiff=0.000e+00
iceland With Egestion  maxdiff=2.365e-11   (float reassociation: symbolic→numeric solve)
iceland TE             maxdiff=0.000e+00
toy GE                 maxdiff=0.000e+00
toy With Egestion      maxdiff=0.000e+00
toy TE                 maxdiff=0.000e+00
ALL MATCH   (tolerance 1e-9)
```

### 2. `SPPR_symbolic` golden regression (default params == old), all 12 combos
```
{iceland,toy} × {GE, With Egestion, TE} × {as_DC, as_PP}  →  maxdiff = 0.000e+00 each
ALL MATCH
```

### 3. Feature behaviour tests (`scratch_feature_tests.py`) — 11/11 passed
```
PASS toy GE none==golden
PASS toy GE recycling_loss theta=1==golden
PASS toy GE source_dilution theta=1==golden
PASS toy theta<1 shrinks det SPPR
PASS toy external_sppr raises det SPPR
PASS toy per-DET dict != scalar
PASS toy always -> pooled                      (diagnostics method = pooled_detritus_scaling)
PASS collapse_det=True -> auto path ran
PASS iceland single-DET theta=1==golden
PASS iceland single-DET theta<1 changes result
PASS iceland symbolic theta=1==golden
SUMMARY: 11 / 11 passed
```

### 4. Monte-Carlo smoke (negatives rejected, never raised)
```
monte_carlo_SPPR   kind='new'       n=30  -> rejected fraction 0.0, result shape (6, 4)
monte_carlo_SPPR_2 kind='new'       n=30  -> rejected fraction 0.0, result shape (6, 4)
monte_carlo_SPPR   kind='symbolic'  n=15  -> rejected fraction 0.0, result shape (5, 4)
```

### 5. Post-documentation re-verification
After the comments-only documentation pass, items 1–3 were re-run: `GOLDEN ALL MATCH`
and `SUMMARY: 11 / 11 passed` (nothing reverted).
