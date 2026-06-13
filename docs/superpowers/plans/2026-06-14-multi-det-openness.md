# Multi-DET Openness, Collapse Modes & Diagnostics — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add ecologically-grounded detritus "openness" (GPT review point 8), an explicit `never/auto/always` collapse mode with a real spectral-radius/condition stability check, fate-weighted pooled scaling (point 7), a missing-fate guard (point 6), and run diagnostics (point 5) to all SPPR methods — while keeping default-parameter output byte-identical to today.

**Architecture:** Centralize the detritus recycling system in two new helpers — `_build_det_BC` (constructs the numeric `(I−B)x=c` system from `M0`, `q`, `egestion`, `det_fate`, `DC`) and `_solve_det_scaling` (applies openness transform, checks stability, then solves or pools, stores diagnostics, and scales the DET columns of the SPPR basis). `SPPR_new` and both `SPPR_symbolic` helpers call these instead of their own inline solve. The default path (`det_open_mode='none'`, `det_collapse_mode='never'`) reproduces the existing math exactly; the symbolic default path keeps its proven `sympy` solve untouched and only routes through the numeric helper when a non-default mode is requested.

**Tech Stack:** Python, NumPy, pandas, SymPy. Models loaded via `PPRCalculator(json_path)`.

**Backward-compat invariants (HARD requirements):**
- `monte_carlo_SPPR`/`_2` depend on unstable solves returning **negative** SPPR (they reject those samples). `det_collapse_mode='never'` must therefore **solve and return negatives, never raise**.
- Default outputs for `227_Iceland_(1950)` and `900_Multi_DET_Toy_(2026)` must be identical pre/post change.

**Golden baselines already captured (default params, this repo state):**
| Model | DET groups | SPPR_new GE det-sum | With Egestion det-sum |
|---|---|---|---|
| Iceland 227 | 1 (`Detritus`) | 61456.2338 | 20027.0144 |
| Toy 900 | 2 (`Fast`/`Slow detritus`) | 5.2127 | 5.2920 |

Model paths: `real_models/EwE_jsons/227_Iceland_(1950).json`, `real_models/ToyModels/900_Multi_DET_Toy_(2026).json`.

---

## New public parameter surface (identical across all methods)

```python
det_collapse_mode='never'   # 'never'|'auto'|'always'
det_open_mode='none'        # 'none'|'recycling_loss'|'source_dilution'
det_theta=1.0               # availability θ (or 1-λ); float OR {DET seq-or-name: float}; default 1.0
det_external_sppr=0.0       # external_sppr_l (point 8B knob); float OR dict; default 0.0
collapse_det=None           # back-compat: False->'never', True->'auto'; None means "use det_collapse_mode"
```

Resolution rule (top of every method): `if collapse_det is not None: det_collapse_mode = 'auto' if collapse_det else 'never'`.

Openness math (per-DET vectors θ, ext aligned to `DET_seq`):
- `none`: `B_open=B`, `c_open=c`
- `recycling_loss` (8A): `B_open=diag(θ)·B`, `c_open=c`
- `source_dilution` (8B): `B_open=diag(θ)·B`, `c_open=θ·c + ext·(1−θ)`

Decision: `ρ=max|eig(B_open)|`, `cond=cond(I−B_open)`.
- `always` → pool. `auto` → pool iff `ρ≥1−tol or cond>cond_threshold`, else solve. `never` → always solve (no raise).
- `tol=1e-10`, `cond_threshold=1e10`.

---

## Task 0: Capture full golden baseline (regression oracle)

**Files:** Create `scratch_baseline.py` (temporary, deleted at end), writes `scratch_golden.pkl`.

- [ ] **Step 1: Write the capture script**

```python
# scratch_baseline.py — capture pre-change outputs as the regression oracle
import pickle, numpy as np
from PPRCalculator import PPRCalculator

MODELS = {
    'iceland': 'real_models/EwE_jsons/227_Iceland_(1950).json',
    'toy': 'real_models/ToyModels/900_Multi_DET_Toy_(2026).json',
}
golden = {}
for name, path in MODELS.items():
    c = PPRCalculator(path)
    for opt in ['GE', 'With Egestion', 'TE']:
        s, A, L = c.SPPR_new(TE_option=opt)
        golden[(name, 'new', opt)] = s.sort_index().sort_index(axis=1).values.copy()
        for dio in ['as_DC', 'as_PP']:
            try:
                _, sm_, _, _ = c.SPPR_symbolic(TE_option=opt, diet_import_option=dio)
                golden[(name, 'sym', opt, dio)] = sm_.sort_index().sort_index(axis=1).values.copy()
            except Exception as e:
                golden[(name, 'sym', opt, dio)] = ('ERROR', repr(e))
with open('scratch_golden.pkl', 'wb') as f:
    pickle.dump(golden, f)
print('captured', len(golden), 'baseline arrays')
for k, v in golden.items():
    print(k, 'ERROR' if isinstance(v, tuple) else np.asarray(v).shape)
```

- [ ] **Step 2: Run it BEFORE any source change**

Run: `python scratch_baseline.py`
Expected: prints "captured N baseline arrays" and a shape per key. Note which symbolic combos error today (those errors are part of the baseline — they must keep erroring identically, not newly crash differently).

- [ ] **Step 3: Commit the plan only** (not the scratch file)

```bash
git add docs/superpowers/plans/2026-06-14-multi-det-openness.md
git commit -m "docs: plan for multi-DET openness + collapse modes"
```

---

## Task 1: Add `_resolve_det_param` and `_spectral_radius` helpers

**Files:** Modify `PPRCalculator.py` (add two methods near the other helpers, e.g. just above `_collapse_det_scaling` at L1055).

- [ ] **Step 1: Add the helpers**

```python
    @staticmethod
    def _spectral_radius(M):
        """Largest absolute eigenvalue of M (0 for empty). Used to test whether the
        detritus recycling matrix B is subcritical (rho < 1 => finite solution)."""
        M = np.asarray(M, dtype=float)
        if M.size == 0:
            return 0.0
        return float(np.max(np.abs(np.linalg.eigvals(M))))

    def _resolve_det_param(self, param, DET_seq, default):
        """Resolve a per-DET parameter into an np.array aligned with DET_seq.
        Accepts a scalar (broadcast to all DET groups) or a dict keyed by DET group
        seq (int) or DET group name (str). Missing keys fall back to `default`."""
        if param is None:
            param = default
        if np.isscalar(param):
            return np.full(len(DET_seq), float(param), dtype=float)
        out = np.full(len(DET_seq), float(default), dtype=float)
        for i, d in enumerate(DET_seq):
            if d in param:
                out[i] = float(param[d])
            elif self.seq2name.get(d) in param:
                out[i] = float(param[self.seq2name[d]])
        return out
```

- [ ] **Step 2: Smoke test in a REPL**

Run:
```bash
python -c "
from PPRCalculator import PPRCalculator
c = PPRCalculator('real_models/ToyModels/900_Multi_DET_Toy_(2026).json')
ds = c.get_DET_seq()
print(c._resolve_det_param(0.5, ds, 1.0))
print(c._resolve_det_param({'Fast detritus':0.8,'Slow detritus':0.4}, ds, 1.0))
import numpy as np
print(c._spectral_radius(np.array([[0.2,0.1],[0.0,0.3]])))
"
```
Expected: `[0.5 0.5]`, a length-2 array with 0.8/0.4 mapped to the right DET seqs, and `0.3`.

- [ ] **Step 3: Commit**

```bash
git add PPRCalculator.py && git commit -m "feat: add _spectral_radius and _resolve_det_param helpers"
```

---

## Task 2: Add `_build_det_BC` (numeric recycling system builder)

**Files:** Modify `PPRCalculator.py` (add method near the helpers). Add `import warnings` at top if absent.

This centralizes the `B`,`c` construction currently duplicated in `SPPR_new`'s GE and With-Egestion branches, and implements point 6 (missing-fate guard).

- [ ] **Step 1: Add the method**

```python
    def _build_det_BC(self, sppr_basis, non_DET_sppr, DET_seq, DC, TE_option):
        """Build the detritus recycling system (I - B) x = c for GE / With Egestion.

        x_l        = scaling factor for detritus group l
        c[l]       = m_eff_l . non_DET_sppr      (production from non-detritus sources)
        B[l, j]    = m_eff_l . sppr_basis[det_j] (recursive dependence on detritus basis j)

        where m_eff_l = M0*fracs/q_l                 (GE)
                      = M0*fracs/q_l + DC.T @ (egestion*fracs/q_l)  (With Egestion)
        and fracs = det_fate[:, det_l] (fraction of each group's flow_to_det reaching det_l).

        Point 6: in a MULTI-DET model, a det_fate matrix that is present but missing a
        column for det_l means nothing feeds det_l -> fracs = 0 (warn), NOT whole-flow.
        For single-DET (or no det_fate) fracs defaults to 1 (old behavior)."""
        k = len(DET_seq)
        det_fate = getattr(self, '_det_fate', None)
        B = np.zeros((k, k))
        c = np.zeros(k)
        for li, det_l in enumerate(DET_seq):
            q_l = self.q[det_l] if self.q[det_l] > 0 else 1.0
            if det_fate is not None and det_l in det_fate.columns:
                fracs = det_fate[det_l].reindex(self.M0.index).fillna(0)
            elif det_fate is not None and k > 1:
                warnings.warn(
                    f"det_fate has no column for DET group {det_l} "
                    f"({self.seq2name.get(det_l)}); treating its detritus inflow as 0.",
                    RuntimeWarning,
                )
                fracs = pd.Series(0.0, index=self.M0.index)
            else:
                fracs = pd.Series(1.0, index=self.M0.index)
            m_l = (self.M0 * fracs / q_l).fillna(0)
            if TE_option == 'With Egestion':
                e_l = (self.egestion * fracs / q_l).fillna(0)
                m_eff_l = m_l + DC.T @ e_l
            else:  # GE
                m_eff_l = m_l
            c[li] = float(m_eff_l @ non_DET_sppr)
            for ji, det_j in enumerate(DET_seq):
                B[li, ji] = float(m_eff_l @ sppr_basis[det_j])
        return B, c
```

- [ ] **Step 2: Verify it reproduces the OLD B/c on the toy model (GE + With Egestion)**

Run:
```bash
python -c "
import numpy as np
from PPRCalculator import PPRCalculator
c = PPRCalculator('real_models/ToyModels/900_Multi_DET_Toy_(2026).json')
DET=c.get_DET_seq()
DC=c.get_DC(DET_as_PP=True, normalize=False)
GE=c.get_TE(TE_option='GE', DET_values=1, as_matrix=True)
# reconstruct the SPPR basis exactly as SPPR_new does up to the DET step is complex;
# instead assert _build_det_BC GE on a known basis matches manual M0/q formula for k=1 sanity:
print('builds without error; DET=',DET)
"
```
Expected: prints DET groups, no error. (Full equivalence is proven in Task 5's golden check.)

- [ ] **Step 3: Commit**

```bash
git add PPRCalculator.py && git commit -m "feat: add _build_det_BC recycling-system builder (point 6 guard)"
```

---

## Task 3: Upgrade `_collapse_det_scaling` (points 7 + openness, return diagnostics)

**Files:** Modify `PPRCalculator.py:1055-1091`.

- [ ] **Step 1: Replace the staticmethod**

```python
    @staticmethod
    def _collapse_det_scaling(SPPR, DET_seq, non_DET_sppr, M0, egestion, q,
                              flow_to_det, DC, TE_option, det_fate=None,
                              theta=None, ext=None, det_open_mode='none'):
        """Fallback DET scaling: treat all DET groups as one pooled pool, solve the 1-D
        self-consistency x = a + b*x, and multiply every DET column by that scalar.

        Point 7: weight M0/egestion by `fate_to_modeled_det` = sum over modeled DET columns
        of det_fate, so the numerator only counts material that actually enters the modeled
        detritus system (consistent with q_combined, which is already fate-weighted).
        With det_fate rows summing to 1 this weight is 1 (no-op).

        Openness is applied to the pooled scalar using the mean θ / ext across DET groups.
        Returns (SPPR, sppr_det, a, b)."""
        q_combined = float(sum(float(q[d]) for d in DET_seq if float(q[d]) > 0))
        if q_combined <= 0:
            q_combined = float(flow_to_det.sum())

        if det_fate is not None:
            fate_to_modeled = (det_fate.reindex(index=M0.index, columns=list(DET_seq))
                               .fillna(0).sum(axis=1).clip(lower=0.0, upper=1.0))
        else:
            fate_to_modeled = pd.Series(1.0, index=M0.index)

        m_combined = (M0 * fate_to_modeled / q_combined).fillna(0)
        if TE_option == 'With Egestion':
            e_combined = (egestion * fate_to_modeled / q_combined).fillna(0)
            m_eff = m_combined + DC.T @ e_combined
        else:  # GE
            m_eff = m_combined

        SPPR_combined_raw = SPPR[list(DET_seq)].sum(axis=1)
        a = float(m_eff @ non_DET_sppr)
        b = float(m_eff @ SPPR_combined_raw)

        if det_open_mode != 'none' and theta is not None:
            th = float(np.mean(theta))
            if det_open_mode == 'recycling_loss':
                b = th * b
            elif det_open_mode == 'source_dilution':
                ex = float(np.mean(ext)) if ext is not None else 0.0
                a = th * a + ex * (1.0 - th)
                b = th * b

        if b >= 1.0:
            raise ValueError(
                f"collapse_det fallback also diverges (b={b:.4f} >= 1). "
                "Detrital cycling is too strong for SPPR_new."
            )
        sppr_det = a / (1.0 - b)
        for det_j in DET_seq:
            SPPR[det_j] *= sppr_det
        return SPPR, sppr_det, a, b
```

- [ ] **Step 2: Update the symbolic call sites that currently expect a bare SPPR return**

At `PPRCalculator.py` ~L1325 and ~L1445 the calls are `sppr_mat = self._collapse_det_scaling(...)`. These will be replaced in Task 6; for now (so the file still runs) wrap them:
`sppr_mat = self._collapse_det_scaling(...)[0]`

Run: `python -c "import PPRCalculator"` — Expected: imports clean.

- [ ] **Step 3: Commit**

```bash
git add PPRCalculator.py && git commit -m "feat: fate-weighted pooled scaling (point 7) + openness + diagnostics return"
```

---

## Task 4: Add `_solve_det_scaling` (openness + stability decision + apply + diagnostics)

**Files:** Modify `PPRCalculator.py` (add method near the others).

- [ ] **Step 1: Add the method**

```python
    def _solve_det_scaling(self, B, c_vec, DET_seq, SPPR, non_DET_sppr, DC, TE_option,
                           det_collapse_mode='never', det_open_mode='none',
                           det_theta=1.0, det_external_sppr=0.0,
                           tol=1e-10, cond_threshold=1e10):
        """Apply openness, decide solve-vs-pool by spectral radius / conditioning, scale the
        DET columns of SPPR in place, and record self.detritus_resolution_info.

        det_collapse_mode: 'never' (always solve; may return negatives, never raises),
                           'auto'  (pool iff rho>=1-tol or cond>cond_threshold),
                           'always'(always pool).
        det_open_mode:     'none' | 'recycling_loss' | 'source_dilution'."""
        k = len(DET_seq)
        theta = self._resolve_det_param(det_theta, DET_seq, 1.0)
        ext = self._resolve_det_param(det_external_sppr, DET_seq, 0.0)

        if det_open_mode == 'none':
            B_open, c_open = B.copy(), c_vec.copy()
        elif det_open_mode == 'recycling_loss':
            B_open, c_open = np.diag(theta) @ B, c_vec.copy()
        elif det_open_mode == 'source_dilution':
            B_open = np.diag(theta) @ B
            c_open = theta * c_vec + ext * (1.0 - theta)
        else:
            raise ValueError("det_open_mode must be 'none', 'recycling_loss', or 'source_dilution'")

        IminusB = np.eye(k) - B_open
        rho = self._spectral_radius(B_open)
        try:
            cond = float(np.linalg.cond(IminusB))
        except np.linalg.LinAlgError:
            cond = np.inf

        if det_collapse_mode == 'always':
            use_collapse, reason = True, 'mode_always'
        elif det_collapse_mode == 'auto':
            if rho >= 1.0 - tol:
                use_collapse, reason = True, 'spectral_radius_ge_1'
            elif cond > cond_threshold:
                use_collapse, reason = True, 'ill_conditioned'
            else:
                use_collapse, reason = False, None
        elif det_collapse_mode == 'never':
            use_collapse, reason = False, None
        else:
            raise ValueError("det_collapse_mode must be 'never', 'auto', or 'always'")

        det_fate = getattr(self, '_det_fate', None)
        flow_to_det = (self.M0 + self.egestion).fillna(0)

        if use_collapse:
            SPPR, scalar, a_p, b_p = self._collapse_det_scaling(
                SPPR, DET_seq, non_DET_sppr, self.M0, self.egestion, self.q,
                flow_to_det, DC, TE_option, det_fate=det_fate,
                theta=theta, ext=ext, det_open_mode=det_open_mode)
            self.detritus_resolution_info = {
                'method': 'pooled_detritus_scaling', 'reason': reason,
                'rho_B': rho, 'cond_IminusB': cond,
                'det_seq': list(DET_seq), 'det_names': [self.seq2name[d] for d in DET_seq],
                'open_mode': det_open_mode, 'theta': theta.tolist(),
                'collapse_scalar': scalar, 'collapse_a': a_p, 'collapse_b': b_p,
            }
        else:
            x_vec = np.linalg.solve(IminusB, c_open)
            for i, det_j in enumerate(DET_seq):
                SPPR[det_j] *= x_vec[i]
            self.detritus_resolution_info = {
                'method': 'multi_detritus' if k > 1 else 'single_detritus', 'reason': reason,
                'rho_B': rho, 'cond_IminusB': cond,
                'det_seq': list(DET_seq), 'det_names': [self.seq2name[d] for d in DET_seq],
                'open_mode': det_open_mode, 'theta': theta.tolist(),
                'x_vec': x_vec.tolist(),
            }
        return SPPR
```

- [ ] **Step 2: Import check** — `python -c "import PPRCalculator"` → clean.

- [ ] **Step 3: Commit**

```bash
git add PPRCalculator.py && git commit -m "feat: add _solve_det_scaling (openness + stability decision + diagnostics)"
```

---

## Task 5: Rewire `SPPR_new` to use the helpers

**Files:** Modify `PPRCalculator.py:1093-1236`.

- [ ] **Step 1: Change the signature**

```python
    def SPPR_new(self, TE=None, TE_option='GE', DET_TE_vals=1, collapse_det=None,
                 det_collapse_mode='never', det_open_mode='none',
                 det_theta=1.0, det_external_sppr=0.0):
```

Immediately after the docstring/first line add the back-compat shim:
```python
        if collapse_det is not None:
            det_collapse_mode = 'auto' if collapse_det else 'never'
```

- [ ] **Step 2: Replace the entire `if len(DET_seq) == 1: ... else: ...` block (current L1144-1234)** with a unified path:

```python
        DET_seq = self.get_DET_seq()
        flow_to_det = (self.M0 + self.egestion).fillna(0)
        non_DET_sppr = SPPR.drop(columns=DET_seq).sum(axis=1)

        if TE_option == 'TE':
            # No detrital recycling matrix in the TE option: scale each DET column by the
            # direct PP+Import inflow share, then apply availability theta as a plain multiplier.
            PP_Import_seq = list(self.get_PP_seq() + self.get_Import_seq())
            det_fate = getattr(self, '_det_fate', None)
            theta = self._resolve_det_param(det_theta, DET_seq, 1.0)
            for i, det_l in enumerate(DET_seq):
                q_l = self.q[det_l] if self.q[det_l] > 0 else flow_to_det.sum()
                if det_fate is not None and det_l in det_fate.columns:
                    fr = det_fate[det_l].reindex(flow_to_det.index).fillna(0)
                    inflow = (flow_to_det[PP_Import_seq] * fr[PP_Import_seq]).sum()
                elif det_fate is not None and len(DET_seq) > 1:
                    inflow = 0.0
                else:
                    inflow = flow_to_det[PP_Import_seq].sum()
                SPPR[det_l] *= (inflow / q_l) * theta[i]
        elif TE_option in ('GE', 'With Egestion'):
            B, c_vec = self._build_det_BC(SPPR, non_DET_sppr, DET_seq, DC, TE_option)
            SPPR = self._solve_det_scaling(
                B, c_vec, DET_seq, SPPR, non_DET_sppr, DC, TE_option,
                det_collapse_mode=det_collapse_mode, det_open_mode=det_open_mode,
                det_theta=det_theta, det_external_sppr=det_external_sppr)
        else:
            raise Exception("TE_option should be in ['GE', 'TE', 'With Egestion', 'global']")

        return SPPR, A, L
```

Note: the single-DET TE branch previously divided by `flow2det = q[det_j] or flow_to_det.sum()` and used PP+Import. The unified TE branch above preserves that for single-DET (det_fate present with the one column => `fr`≈1) and adds the multi-DET point-6 guard.

- [ ] **Step 2b: Run golden regression for SPPR_new**

```bash
python -c "
import pickle, numpy as np
from PPRCalculator import PPRCalculator
g = pickle.load(open('scratch_golden.pkl','rb'))
P = {'iceland':'real_models/EwE_jsons/227_Iceland_(1950).json','toy':'real_models/ToyModels/900_Multi_DET_Toy_(2026).json'}
ok=True
for name,path in P.items():
    c=PPRCalculator(path)
    for opt in ['GE','With Egestion','TE']:
        s,_,_=c.SPPR_new(TE_option=opt)
        new=s.sort_index().sort_index(axis=1).values
        old=g[(name,'new',opt)]
        d=np.max(np.abs(new-old))
        print(name,opt,'maxdiff=%.3e'%d)
        ok = ok and d < 1e-9
print('ALL MATCH' if ok else 'MISMATCH!!')
"
```
Expected: every maxdiff < 1e-9 and `ALL MATCH`. If MISMATCH, do NOT proceed — diff the construction; the most likely culprit is the single-DET fracs vs old `flow2det`. Fix until identical.

- [ ] **Step 3: Commit**

```bash
git add PPRCalculator.py && git commit -m "refactor: SPPR_new uses shared det-scaling helpers; defaults unchanged"
```

---

## Task 6: Rewire both `SPPR_symbolic` helpers + `SPPR_symbolic`

**Files:** Modify `PPRCalculator.py` `_SPPR_symbolic_helper_diet_import_as_PP` (~L1238) and `_SPPR_symbolic_helper_diet_import_as_DC` (~L1340) and `SPPR_symbolic` (~L1465).

Principle (lowest-risk): the DEFAULT path (`det_open_mode='none'` and effective mode `'never'`) keeps the existing `sympy` `sol_dict2` scaling untouched — byte identical. Any non-default mode routes through the numeric helpers.

- [ ] **Step 1: Thread params through `SPPR_symbolic`**

```python
    def SPPR_symbolic(self, TE=None, TE_option='GE', diet_import_option='as_DC', DET_TE_vals=1,
                      sppr_det_value=None, collapse_det=None,
                      det_collapse_mode='never', det_open_mode='none',
                      det_theta=1.0, det_external_sppr=0.0):
        if collapse_det is not None:
            det_collapse_mode = 'auto' if collapse_det else 'never'
        kwargs = dict(TE=TE, TE_option=TE_option, DET_TE_vals=DET_TE_vals,
                      sppr_det_value=sppr_det_value, det_collapse_mode=det_collapse_mode,
                      det_open_mode=det_open_mode, det_theta=det_theta,
                      det_external_sppr=det_external_sppr)
        if diet_import_option == 'as_PP':
            return self._SPPR_symbolic_helper_diet_import_as_PP(**kwargs)
        elif diet_import_option == 'as_DC':
            return self._SPPR_symbolic_helper_diet_import_as_DC(**kwargs)
```

- [ ] **Step 2: Update both helper signatures** to accept the same new kwargs (replace `collapse_det=False` with the four new params):

```python
    def _SPPR_symbolic_helper_diet_import_as_PP(self, TE, TE_option, DET_TE_vals, sppr_det_value,
                                                det_collapse_mode='never', det_open_mode='none',
                                                det_theta=1.0, det_external_sppr=0.0):
```
(identical for the `_as_DC` helper).

- [ ] **Step 3: Replace the collapse decision block in EACH helper**

In `_as_PP` the current block is `PPRCalculator.py` ~L1320-1333; in `_as_DC` ~L1440-1453. Replace each with:

```python
        if sppr_det_value is not None:
            sppr_mat[DET_seq] *= float(sppr_det_value)
        elif det_open_mode == 'none' and det_collapse_mode == 'never':
            # DEFAULT: keep the proven exact symbolic per-DET scaling (byte-identical).
            for det_j in DET_seq:
                det_sym = sppr_vec.loc[det_j].values.ravel()[0]
                sppr_mat[det_j] *= float(sol_dict2[det_sym].evalf(subs=subs_dict))
        else:
            # Non-default: build the numeric (I-B)x=c from the symbolic basis sppr_mat and
            # route through the shared openness/stability/collapse solver.
            non_DET_sppr = sppr_mat.drop(columns=list(DET_seq), errors='ignore').sum(axis=1)
            B, c_vec = self._build_det_BC(sppr_mat, non_DET_sppr, DET_seq, DC, TE_option)
            sppr_mat = self._solve_det_scaling(
                B, c_vec, DET_seq, sppr_mat, non_DET_sppr, DC, TE_option,
                det_collapse_mode=det_collapse_mode, det_open_mode=det_open_mode,
                det_theta=det_theta, det_external_sppr=det_external_sppr)
```

Note for `_as_DC`: the DET-column scaling there used `sol_dict2[det_sym].evalf(subs=subs_dict)` with its own `subs_dict`; keep that exact expression in the default branch (copy from the existing line). Leave the subsequent `DIET_*` handling untouched.

- [ ] **Step 4: Golden regression for symbolic (default params)**

```bash
python -c "
import pickle, numpy as np
from PPRCalculator import PPRCalculator
g = pickle.load(open('scratch_golden.pkl','rb'))
P = {'iceland':'real_models/EwE_jsons/227_Iceland_(1950).json','toy':'real_models/ToyModels/900_Multi_DET_Toy_(2026).json'}
ok=True
for name,path in P.items():
    c=PPRCalculator(path)
    for opt in ['GE','With Egestion','TE']:
        for dio in ['as_DC','as_PP']:
            old=g[(name,'sym',opt,dio)]
            try:
                _,s,_,_=c.SPPR_symbolic(TE_option=opt, diet_import_option=dio)
                new=s.sort_index().sort_index(axis=1).values
            except Exception as e:
                new=('ERROR',repr(e))
            if isinstance(old,tuple):
                match = isinstance(new,tuple)  # both errored
                print(name,opt,dio,'OLD-ERROR new-error?',match)
            else:
                d=np.max(np.abs(new-old)); print(name,opt,dio,'maxdiff=%.3e'%d); match = d<1e-9
            ok = ok and match
print('ALL MATCH' if ok else 'MISMATCH!!')
"
```
Expected: `ALL MATCH`. Combos that errored in the baseline must still error (not newly crash). If a previously-working combo mismatches, fix before proceeding.

- [ ] **Step 5: Commit**

```bash
git add PPRCalculator.py && git commit -m "fix: symbolic path shares stability gate (point 4); default byte-identical"
```

---

## Task 7: Thread new params through Monte-Carlo

**Files:** Modify `PPRCalculator.py` `monte_carlo_SPPR` (~L1483) and `monte_carlo_SPPR_2` (~L1574).

- [ ] **Step 1: Add params to both signatures** (after existing args):

```python
                 det_collapse_mode='never', det_open_mode='none',
                 det_theta=1.0, det_external_sppr=0.0,
```

- [ ] **Step 2: Forward to every internal `SPPR_new` / `SPPR_symbolic` call** inside both methods (initial collector call + the in-loop calls). Example for the `SPPR_new` calls:

```python
                sppr, _, _ = self.SPPR_new(TE=TE_sample, TE_option=TE_option, DET_TE_vals=DET_TE_vals,
                                           det_collapse_mode=det_collapse_mode, det_open_mode=det_open_mode,
                                           det_theta=det_theta, det_external_sppr=det_external_sppr)
```
and analogously for the `SPPR_symbolic` call in `monte_carlo_SPPR` (pass the same four kwargs plus its existing `diet_import_option`).

- [ ] **Step 3: Run a tiny Monte-Carlo to prove negatives are still rejected, not raised**

```bash
python -c "
from PPRCalculator import PPRCalculator
c = PPRCalculator('real_models/ToyModels/900_Multi_DET_Toy_(2026).json')
s, arr, frac, *_ = c.monte_carlo_SPPR(n_samples=30, TE_option='GE', kind='new')
print('ran MC; rejected fraction =', frac, '; result shape', s.shape)
"
```
Expected: runs without raising; prints a rejected fraction and a result shape.

- [ ] **Step 4: Commit**

```bash
git add PPRCalculator.py && git commit -m "feat: thread det openness/collapse params through Monte-Carlo"
```

---

## Task 8: Feature behavior tests (openness + collapse modes)

**Files:** Create `scratch_feature_tests.py` (temporary). These assert the NEW knobs behave as designed and that no-dampening == old.

- [ ] **Step 1: Write the checks**

```python
# scratch_feature_tests.py
import numpy as np, pickle
from PPRCalculator import PPRCalculator
g = pickle.load(open('scratch_golden.pkl','rb'))
toy='real_models/ToyModels/900_Multi_DET_Toy_(2026).json'
ice='real_models/EwE_jsons/227_Iceland_(1950).json'
results=[]
def check(name, cond):
    results.append((name,bool(cond))); print(('PASS' if cond else 'FAIL'), name)

c=PPRCalculator(toy); DET=c.get_DET_seq()
# 1) open_mode none == old
s0,_,_=c.SPPR_new(TE_option='GE', det_open_mode='none')
check('toy GE none==golden', np.max(np.abs(s0.sort_index().sort_index(axis=1).values-g[('toy','new','GE')]))<1e-9)
# 2) theta=1 == old (recycling_loss and source_dilution)
for mode in ['recycling_loss','source_dilution']:
    s=c.SPPR_new(TE_option='GE', det_open_mode=mode, det_theta=1.0, det_external_sppr=0.0)[0]
    check(f'toy GE {mode} theta=1==golden', np.max(np.abs(s.sort_index().sort_index(axis=1).values-g[('toy','new','GE')]))<1e-9)
# 3) theta<1 shrinks detritus SPPR (recycling_loss)
base=c.SPPR_new(TE_option='GE')[0][DET].sum().sum()
shrunk=c.SPPR_new(TE_option='GE', det_open_mode='recycling_loss', det_theta=0.5)[0][DET].sum().sum()
check('toy theta<1 shrinks det SPPR', shrunk < base)
# 4) external_sppr>0 raises det SPPR back up (source_dilution)
low=c.SPPR_new(TE_option='GE', det_open_mode='source_dilution', det_theta=0.5, det_external_sppr=0.0)[0][DET].sum().sum()
high=c.SPPR_new(TE_option='GE', det_open_mode='source_dilution', det_theta=0.5, det_external_sppr=5.0)[0][DET].sum().sum()
check('toy external_sppr raises det SPPR', high > low)
# 5) per-DET dict differs from scalar
d_scalar=c.SPPR_new(TE_option='GE', det_open_mode='recycling_loss', det_theta=0.5)[0][DET].sum().sum()
d_dict=c.SPPR_new(TE_option='GE', det_open_mode='recycling_loss', det_theta={'Fast detritus':0.5,'Slow detritus':0.9})[0][DET].sum().sum()
check('toy per-DET dict != scalar', abs(d_dict-d_scalar)>1e-9)
# 6) always-collapse gives equal scalar across DET columns
sa=c.SPPR_new(TE_option='GE', det_collapse_mode='always')[0]
# ratio of each DET column to its golden basis should be equal across DET groups
print('  diagnostics after always:', c.detritus_resolution_info['method'])
check('toy always -> pooled', c.detritus_resolution_info['method']=='pooled_detritus_scaling')
# 7) back-compat: collapse_det=True maps to auto
c.SPPR_new(TE_option='GE', collapse_det=True)
check('collapse_det=True -> auto path ran', 'method' in c.detritus_resolution_info)
# 8) single-DET openness works (Iceland) and theta=1 == old
ci=PPRCalculator(ice)
si=ci.SPPR_new(TE_option='GE', det_open_mode='recycling_loss', det_theta=1.0)[0]
check('iceland single-DET theta=1==golden', np.max(np.abs(si.sort_index().sort_index(axis=1).values-g[('iceland','new','GE')]))<1e-9)
si2=ci.SPPR_new(TE_option='GE', det_open_mode='recycling_loss', det_theta=0.5)[0]
check('iceland single-DET theta<1 changes result', np.max(np.abs(si2.values-si.values))>1e-9)
# 9) symbolic non-default routes through solver (no crash) and theta=1 default==golden
_,sm1,_,_=ci.SPPR_symbolic(TE_option='GE', diet_import_option='as_DC', det_open_mode='recycling_loss', det_theta=1.0)
check('iceland symbolic theta=1==golden', np.max(np.abs(sm1.sort_index().sort_index(axis=1).values-g[('iceland','sym','GE','as_DC')]))<1e-9)

print('\\nSUMMARY:', sum(1 for _,ok in results if ok), '/', len(results), 'passed')
assert all(ok for _,ok in results), 'feature tests failed'
```

- [ ] **Step 2: Run**

Run: `python scratch_feature_tests.py`
Expected: all PASS, `SUMMARY: N / N passed`. Capture this output for the change log. If a behavioral assertion fails, fix the implementation (not the test, unless the test's expectation is wrong).

- [ ] **Step 3: Commit** (source only; scratch files are not committed)

```bash
git add PPRCalculator.py && git commit -m "test: verify openness/collapse behavior + no-dampening==old"
```

---

## Task 9: Documentation pass over all of `PPRCalculator.py` (subagent)

**Files:** Modify `PPRCalculator.py` — comments only, no behavior change.

- [ ] **Step 1: Dispatch a subagent** (general-purpose) with explicit constraints:
  - Add clear explanatory comments to every method in `PPRCalculator.py` (constructors, getters, `apply_ecopath_defaults`, `apply_lim`, `balance_model`, all `SPPR_*`, helpers, Monte-Carlo, class methods).
  - MUST NOT change any executable code, signatures, defaults, string literals, or whitespace that affects execution — comments and docstrings ONLY.
  - Must run `graphify query` before exploring, per project rules.
  - After editing, run `python -c "import PPRCalculator"` and the Task 5 + Task 6 golden scripts; all must still say ALL MATCH.

- [ ] **Step 2: Re-run BOTH golden scripts (Task 5 Step 2b and Task 6 Step 4) after the doc pass**

Expected: `ALL MATCH` for both. If anything moved, the doc pass changed code — revert offending hunk.

- [ ] **Step 3: Commit**

```bash
git add PPRCalculator.py && git commit -m "docs: explanatory comments across PPRCalculator.py (no behavior change)"
```

---

## Task 10: Change log, cleanup, git, graphify

- [ ] **Step 1: Write `CHANGELOG_multi_det_openness.md`** at repo root summarizing: each point addressed (2,3,4,5,6,7,8), new params + defaults, the openness math, files/methods touched, and a **"Tests executed"** section pasting the golden-regression results (Tasks 5/6) and the feature-test summary (Task 8), plus the Monte-Carlo smoke result.

- [ ] **Step 2: Delete scratch files**

```bash
rm -f scratch_baseline.py scratch_feature_tests.py scratch_golden.pkl
```
(Keep `response_for_GPT.txt` and `detritus_SPPR_review_for_claude_code.txt`.)

- [ ] **Step 3: Update knowledge graph**

```bash
python -m graphify update .
```
(Use the saved interpreter from `graphify-out/.graphify_python` if `graphify` is not on PATH.)

- [ ] **Step 4: Final commit**

```bash
git add -A && git commit -m "docs: change log for multi-DET openness; update knowledge graph"
```

- [ ] **Step 5: Report** the branch state and a summary of all commits to the user.

---

## Self-review notes
- **Spec coverage:** point 2 (Task 5/6 mode + shim), point 3 (Task 4 rho+cond), point 4 (Task 6 shared gate), point 5 (Task 4 diagnostics), point 6 (Task 2 guard), point 7 (Task 3 fate weighting), point 8 (Task 4 openness, both modes + external_sppr). Monte-Carlo (Task 7). Single-DET openness (Task 5 unified path + Task 8 check 8). All covered.
- **No-dampening == old:** enforced by golden checks in Tasks 5, 6, 9 and assertions 1/2/8/9 in Task 8.
- **Monte-Carlo no-raise invariant:** `never` never raises (Task 4); proven in Task 7 Step 3.
- **Naming consistency:** `_build_det_BC`, `_solve_det_scaling`, `_collapse_det_scaling`, `_resolve_det_param`, `_spectral_radius`, `detritus_resolution_info` used consistently across tasks.
