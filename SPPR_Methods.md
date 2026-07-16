# SPPR Calculation Methods in `PPRCalculator`

This document explains every SPPR (**Specific Primary Production Required**) method exposed by
`PPRCalculator.py`, the flags and input options each one accepts, and — crucially — the
biological / ecological meaning behind each choice. It also documents the `get_PPR`,
`get_NPP` and `get_PPR2NPP_ratio` functions that turn a per-group SPPR into an
ecosystem-level footprint.

---

## 1. What SPPR means

**PPR** (Primary Production Required) is the amount of primary production (plant/algal
production at the base of the food web) that must be fixed to ultimately support a given
harvest. **SPPR** is the *specific* version: the primary production required **per unit of
production of a given group**. If a group sits at trophic level 4 and energy is lost at every
transfer up the chain, then producing 1 tonne of that group requires many tonnes of primary
production far below it.

Concretely, every SPPR method returns, for each consumer group, **how many units of basal
production** (primary production, and sometimes detritus / imported production) are required
per unit of that group's own production. Multiplying by the catch and summing gives the total
primary production required to sustain the fishery — the ecological "cost" of the catch.

Two broad families of methods exist:

- **Trophic-level methods** (`SPPR_1986`, `SPPR_1995`, `SPPR_1995_TL_fix`) — collapse the
  whole food web into a single number per group (its trophic level) and apply a fixed
  transfer efficiency. Fast, classic, coarse.
- **Flow-network methods** (`SPPR_EwE`, `SPPR_EwE_Ido`, `SPPR_2015`, `SPPR_new`,
  `SPPR_symbolic`) — trace production back through the actual diet matrix, resolving how much
  of *each* basal source (each primary producer, each detritus pool, imported food) is drawn
  down. These give a **per-basal-source breakdown**, not just a single number, and respect the
  real topology of the modelled ecosystem.

---

## 2. Shared building blocks

All methods draw on the same underlying quantities, unpacked from the Ecopath model in
`_fill_properties`:

| Symbol | Attribute | Ecological meaning |
|--------|-----------|--------------------|
| `p` | production | total production of a group (biomass produced per unit time) |
| `q` | consumption | total food eaten by a group |
| `M0` | natural mortality flow | production dying of causes other than predation/fishing → goes to detritus |
| `egestion` | egestion | unassimilated food (faeces) → goes to detritus |
| `EE` | ecotrophic efficiency | fraction of production used *within* the system (eaten or exported), the rest is `M0` |
| `GE` | gross efficiency = `p/q` | how efficiently consumed food becomes production |
| `catch` | catch | the harvest we are attributing a PPR cost to |
| `DC` | diet composition | `DC[i,j]` = fraction of predator *i*'s diet made up of prey *j* |
| `Z` | flow matrix | `Z = DC · q`, absolute flow of biomass from prey to predator |

### Trophic categories (`trophic_info`)
Groups are classified into four kinds, each with its own `get_*_seq()` accessor:

- **PP** (primary producers) — basal autotrophs; the "source" the whole calculation traces back to.
- **DET** (detritus) — dead organic matter pools; recycled basal sources fed by `M0` + egestion.
- **Import** — external production entering the system (diet imported from outside the model boundary).
- **Regular** — consumers (heterotrophs).

PP, DET and Import are **basal sources**: SPPR is expressed as production required *from these*.

### `get_TE` — transfer efficiency (`TE_option`)
Transfer efficiency is the fraction of production at one level that becomes production at the
next. It is the single most important ecological assumption in every flow-network method.
`get_TE` builds a per-group TE from one of four options:

| `TE_option` | Formula | Ecological meaning |
|-------------|---------|--------------------|
| `'GE'` | `p/q` | **Gross growth efficiency** — production per unit consumed. Ignores that not all production is passed on. This is what the EwE user guide uses. |
| `'TE'` | `(p/q)·(1 − M0/p)` | Gross efficiency times the **ecotrophic fraction**: only the part of production that is actually consumed/exported counts as "transferred". Production that dies naturally (`M0`) is treated as lost from the up-web pathway. Equivalent to `GE·EE`. |
| `'With Egestion'` | `(p/q)·(q/(q−egestion))` | Efficiency computed on **assimilated** intake rather than gross intake, i.e. faeces are removed from the denominator so they are accounted separately (they flow to detritus, not up the chain). |
| `'global'` | one scalar for all groups | A single system-wide TE broadcast to every group — the classic "10% rule" assumption. Controlled by `global_TE`. |

- **`global_TE`** (only used when `TE_option='global'`): either a literal float (e.g. `0.1`)
  or `'mean'`, which computes the catch-weighted (or biomass-weighted if there is no catch)
  mean of the per-group `'TE'` efficiency. Biologically, `'mean'` lets the data set its own
  effective efficiency instead of imposing 10%.
- **`DET_values`** (default 1): the TE assigned to detritus rows. TE=1 means detritus passes
  its content on without a further "trophic loss", because detritus is itself a basal source,
  not a trophic step.
- **`as_matrix`**: return an *n×n* matrix (per-edge weights) vs a length-*n* vector.

### `get_TL` — trophic level
Solves the standard linear-algebra definition `TL = (I − DC)⁻¹ · 1`: each group's trophic
level is 1 plus the diet-weighted mean TL of its prey. Basal sources sit at TL 1.
- **`break_cycles`**: remove cycles (Ulanowicz algorithm) before inverting, so recursive
  loops (e.g. detritus ↔ bacteria) don't make the inversion ill-behaved.
- **`DET_as_PP`**: treat detritus as a basal production source at TL 1.
- **`TE_option`**: how detritus diet rows are redefined before inversion.

### `DET_as_PP` and `normalize` (in `get_DC` / `get_Z`)
These control whether **detritus is a source or a recycling loop**:
- `DET_as_PP=True` — detritus is a basal source (its diet row is treated as terminal). This is
  what the flow-network SPPR methods use so that detritus becomes one of the "columns" you get
  an SPPR breakdown against.
- `DET_as_PP=False` — the detritus row is rebuilt from where its material actually came from
  (`M0` + egestion shares, routed by `det_fate`), i.e. detritus recycles back into the web.

---

## 3. Trophic-level SPPR methods

### `SPPR_1986()` — Pauly & Christensen (1986)
```python
SPPR_1986() -> pd.DataFrame   # single 'sppr' column
```
Computes **one** catch-weighted mean trophic level for the entire catch, then
`SPPR = TE^(1−TL)` with a **fixed TE = 0.1**. Every group receives the same number.

- **No user flags.** TE is hard-wired to the classic 10% rule; TL is broken-cycle and
  DET-as-PP forced.
- **Ecological meaning:** the coarsest estimate — "the average fish in this catch is at TL X,
  so at 10% efficiency per level it costs `10^(TL−1)` units of primary production". Good for
  back-of-envelope comparisons across fisheries; blind to who eats whom.
- Returns all zeros if there is no catch.

### `SPPR_1995(global_TE=0.1)` — Christensen & Pauly (1995)
```python
SPPR_1995(global_TE: str | float = 0.1) -> pd.DataFrame
```
Same `TE^(1−TL)` form, but uses **each group's own continuous trophic level** rather than one
catch-average, with a single global TE.

- **`global_TE`**: the system-wide transfer efficiency — a float (default `0.1`, the 10% rule)
  or `'mean'` (let the model's own efficiencies set it).
- **Ecological meaning:** resolves the fact that different harvested groups sit at different
  trophic levels, so a TL-2 shellfish and a TL-4.5 tuna get very different PPR costs. Still
  assumes a uniform efficiency at every step.

### `SPPR_1995_TL_fix(global_TE=0.1)`
```python
SPPR_1995_TL_fix(global_TE: str | float = 0.1) -> pd.DataFrame
```
A numerical refinement of `SPPR_1995`. Instead of raising `1/TE` to a fractional power
directly, it **linearly interpolates between the two bracketing integer trophic levels**:
`SPPR = (1−frac)·(1/TE)^(TLint−1) + frac·(1/TE)^TLint`.

- **`global_TE`**: as above.
- **Ecological meaning:** a group at TL 3.4 is treated as a mixture of "40% of a TL-4 feeder
  and 60% of a TL-3 feeder", which is closer to how a real diet spans discrete prey levels than
  a smooth exponential. Avoids the discontinuity/curvature artefacts of the direct fractional
  exponent.

---

## 4. Flow-network SPPR methods

These trace production back through the diet matrix and return a **matrix**: rows = groups,
columns = basal sources (each PP group, each detritus pool, and/or Import). Cell `(i, s)` = the
production required from basal source *s* per unit of group *i*'s production.

### `SPPR_EwE(TE_option, use_EE=True, return_paths=True, silent=True)`
```python
SPPR_EwE(TE_option, use_EE=True, return_paths=True, silent=True)
    -> (SPPR, A, paths_dict)
```
The classic **path-enumeration** approach of Ecopath with Ecosim. It enumerates *every simple
path* from each group down to a basal terminal node and sums the product of edge weights along
each path. The per-edge weight matrix is `A = DC/TE` (diet fraction divided by transfer
efficiency = production required per unit passed along that link).

- **`TE_option`** — `'GE'`, `'TE'`, `'With Egestion'`, or `'global'`. Governs how much basal
  production each trophic link implies (see §2). The choice of TE is the dominant ecological
  assumption.
- **`use_EE`** (default True) — scale each group's row by its ecotrophic efficiency. Biologically
  this discounts production that never gets eaten (dies to `M0`) so it isn't double-counted as
  supporting the harvest.
- **`return_paths`** (default True) — if True, use the slower implementation that also returns
  the explicit list of food-chain paths (useful for tracing *which* chains dominate the PPR); if
  False, use the fast vectorized implementation (returns an empty paths dict).
- **`silent`** — suppress progress bars.
- **Returns** `(SPPR, A, paths_dict)`.
- **Ecological meaning:** the most literal reading of "trace the energy back". Its weakness is
  cycles: recycling loops create infinitely many paths, so a depth safety-valve caps
  enumeration. This is why the matrix reformulations below exist.

### `SPPR_EwE_Ido(TE_option, global_TE='mean', use_EE=True)`
```python
SPPR_EwE_Ido(TE_option, global_TE='mean', use_EE=True) -> (SPPR, A, L)
```
A **matrix (nullspace) reformulation** of `SPPR_EwE` that avoids path enumeration. It builds
`A = DC/TE` with cycles removed, replaces each basal (producer) row with an identity row, and
finds the steady-state SPPR as the **nullspace of `L = A − I`** (i.e. `A·x = x`). The nullspace
is RREF-normalized so each output column is anchored to exactly one basal source.

- **`TE_option`** — as above.
- **`global_TE`** — only used when `TE_option='global'`; `'mean'` or a literal float.
- **`use_EE`** — scale rows by EE.
- **Returns** `(SPPR, A, L)`.
- **Ecological meaning:** mathematically equivalent to summing all paths (including cycles,
  handled by cycle removal) but exact and fast. This is the bridge between the EwE path method
  and the newer numeric solver.

### `SPPR_2015(only_pp_det=True)` — the 2015 method
```python
SPPR_2015(only_pp_det=True) -> (SPPR, A, L)
```
A **Leontief input–output** formulation (matrix inversion). Detritus columns are *dissolved*:
the PP-derived fraction of each detritus flow is reassigned back onto the PP groups, leaving
only living compartments. The production-normalized transaction matrix `A` then yields the
production-requirement matrix `L = (I − A)⁻¹`, whose PP columns give the per-group SPPR. A
balancing detritus SPPR is added back at the end.

- **`only_pp_det`** — nominally whether to reassign only the PP-derived fraction of detritus
  back onto PP (the article's choice). **Note: the body forces this to `True`**, so it is
  effectively always on.
- **Returns** `(SPPR, A, L)`.
- **Ecological meaning:** treats the ecosystem like an economy where each group's production
  "requires" inputs from the groups it eats; `(I − A)⁻¹` sums the full direct + indirect
  requirement chain. Consistent, closed-form, and the reference implementation of the 2015
  paper. **Caveat (see project memory):** `SPPR_2015` should *not* be used as an oracle to
  validate multi-detritus methods — argue correctness from the equations instead.

### `SPPR_new(...)` — primary numeric solver
```python
SPPR_new(TE=None, TE_option='GE', DET_TE_vals=1, collapse_det=None,
         det_collapse_mode='never', det_open_mode='none',
         det_theta=1.0, det_external_sppr=0.0, fix_EE_0_cases=True)
    -> (SPPR, A, L)
```
The main, most general numeric SPPR method. Like `SPPR_EwE_Ido` it builds `A = DC/TE`, replaces
basal rows with identity rows, and solves the nullspace of `L = A − I`, RREF-normalized to one
column per basal source. Its novelty is **explicit, tunable detritus handling** — how recycled
dead organic matter is credited as a basal source.

**Core inputs:**
- **`TE`** — supply an explicit TE matrix (e.g. a Monte-Carlo sample). If `None`, built from
  `TE_option`.
- **`TE_option`** — `'GE'` (default), `'TE'`, `'With Egestion'`, `'global'` (see §2). This also
  selects how detritus is resolved:
  - `'TE'`: each detritus column is scaled by its **direct PP+Import inflow share** (detritus
    treated as a pass-through of the primary production that fell into it).
  - `'GE'` / `'With Egestion'`: build the **coupled recycling system** `(I − B)x = c` (via
    `_build_det_BC`) and solve it (via `_solve_det_scaling`). This captures that detritus feeds
    consumers whose mortality feeds detritus again — a genuine recycling loop.
- **`DET_TE_vals`** (default 1) — TE for detritus rows when building the TE matrix.

**Detritus recycling knobs** (the ecological heart of the method):
- **`collapse_det`** — legacy boolean, mapped onto `det_collapse_mode`: `False`→`'never'`,
  `True`→`'auto'`, `None`→leave as given.
- **`det_collapse_mode`** — solve-vs-pool strategy when recycling is strong:
  - `'never'` (default): always solve the coupled system; may return negative SPPR but never raises.
  - `'auto'`: pool all detritus into one compartment only if the coupled system is unstable
    (spectral radius ≥ 1, i.e. recycling would "blow up", or ill-conditioned).
  - `'always'`: always pool.
  - *Ecological meaning:* if detritus recycling is too strong (each unit regenerates ≥1 unit),
    the closed system has no finite solution; pooling merges detritus pools so the combined
    denominator tames the loop.
- **`det_open_mode`** — how "open" detritus recycling is:
  - `'none'` (default): closed recycling (all dead matter is reprocessed within the system).
  - `'recycling_loss'`: damp the recycling matrix `B` by `diag(theta)` — a fraction of detritus
    is lost (buried, exported) rather than recycled.
  - `'source_dilution'`: damp `B` **and** dilute the source toward an external SPPR — detritus
    is partly supplied from outside the modelled system.
- **`det_theta`** (default 1.0) — detritus **availability / retention fraction**. `1.0` = the
  closed system; lower values mean less detritus is actually available to consumers (the rest is
  buried/exported). A float (all pools) or a dict keyed by DET seq or name.
- **`det_external_sppr`** (default 0.0) — the SPPR assigned to externally-sourced detritus under
  `'source_dilution'`. Float or dict.
- **`fix_EE_0_cases`** (default True) — a mass-balance correction for the `'TE'` option. Groups
  with `EE=0` (all their production dies naturally, `M0=p`) get transfer efficiency 0 and are
  severed from the nullspace, which **leaks the primary production they consumed** (it goes
  neither up the web nor back to detritus). When True, that consumed PP is **re-credited to the
  detritus pool**, closing the global PP balance. Only active for **single-detritus** models
  under `'TE'`; a `RuntimeWarning` is emitted whenever EE=0 groups are present. It does *not* fix
  near-singular `0 < EE ≪ 1` groups, whose `SPPR ~ 1/te` blows up (an inherent singularity of the
  TE method — a separate warning is raised).
- **Returns** `(SPPR, A, L)`.

### `SPPR_symbolic(...)` — symbolic solver
```python
SPPR_symbolic(TE=None, TE_option='GE', diet_import_option='as_DC', DET_TE_vals=1,
              sppr_det_value=None, collapse_det=None, det_collapse_mode='never',
              det_open_mode='none', det_theta=1.0, det_external_sppr=0.0,
              fix_EE_0_cases=True) -> (sppr_symbolic, sppr_mat, equations, variables)
```
A **symbolic (SymPy)** counterpart to `SPPR_new`. It writes the per-group SPPR balance
equations `A·x − x = 0` symbolically and solves them exactly, returning both the symbolic
solution and a numeric basis matrix. It shares all of `SPPR_new`'s detritus knobs
(`collapse_det`, `det_collapse_mode`, `det_open_mode`, `det_theta`, `det_external_sppr`,
`fix_EE_0_cases`) with identical meaning.

Its distinctive input is **how imported diet is treated**:
- **`diet_import_option`**:
  - `'as_DC'` (default): imported diet is kept as a **separate production source** with its own
    `DIET_SPPR` solved from a second linear system. Each imported group carries the
    diet-composition-weighted production it requires — i.e. imported food is costed by what it,
    in turn, was made of. Ecologically the most faithful treatment of cross-boundary subsidies.
  - `'as_PP'`: imported diet is treated as **just another primary-production source** (its own
    free SPPR symbol fixed to 1). Simpler; treats external food as if it were free basal
    production entering the system.
- **`sppr_det_value`** — if set, every detritus column is scaled by this fixed value instead of
  being solved. Useful to impose an externally-determined detritus SPPR.
- Other args (`TE`, `TE_option`, `DET_TE_vals`, and the detritus knobs) match `SPPR_new`.
- **`fix_EE_0_cases`** — the EE=0 re-credit; only consumed by the `'as_PP'` + `'TE'` path (a
  no-op for `'as_DC'`, for GE / With Egestion, and for multi-DET). Mirrors `SPPR_new`.
- **Returns** `(sppr_symbolic, sppr_mat, equations, variables)` — the symbolic solution, the
  numeric basis matrix, the full equation system, and the ordered symbol list.
- **Default path** (`det_open_mode='none'`, `det_collapse_mode='never'`) keeps the exact SymPy
  per-detritus scaling and is designed to reproduce `SPPR_new` for matching `TE_option`s.

### Helpers used by the numeric/symbolic solvers (for reference)
- **`_build_det_BC`** — assembles the detritus recycling system `(I − B)x = c`, where `B[l,j]`
  is how much detritus pool *l*'s SPPR depends on pool *j*'s (via `M0`/egestion routed by
  `det_fate`), and `c` is the contribution from non-detritus sources. This encodes the
  "dead matter → consumer → dead matter" loop.
- **`_spectral_radius`** — largest eigenvalue of `B`; `< 1` means recycling converges (a finite,
  non-negative SPPR exists), `≥ 1` means it diverges.
- **`_solve_det_scaling`** — applies the openness transform, tests stability, decides
  solve-vs-pool per `det_collapse_mode`, and scales the detritus columns in place.
- **`_collapse_det_scaling`** — the fallback that pools all detritus into one compartment and
  solves a scalar `x = a + b·x`, used when the coupled system is unstable.
- **`_resolve_det_param`** — expands the `det_theta` / `det_external_sppr` scalar-or-dict knobs
  into per-detritus arrays.

### `_sample_SPPR_new_forced_balance(TE_option='TE', sppr_det=None)`
Runs `SPPR_new` and then solves the single detritus scaling value that forces **exact global
mass balance** (PP inflow = export outflow of catch + growth + net migration weighted by SPPR).
Returns `(sppr, sppr_det)`. Ecologically, it pins the one free recycling degree of freedom so
that total primary production in equals total exported production out.

### `monte_carlo_SPPR(...)` / `monte_carlo_SPPR_2(...)`
Uncertainty propagation: resample the TE matrix from gamma distributions centred on the model's
TEs, recompute SPPR (via `SPPR_new` or `SPPR_symbolic`) many times, discard non-physical
(negative-SPPR) draws, and average. Key knobs: `n_samples`, `TE_error_percent` (the CV of the
TE prior), `TE_error_cut_percent` (clip band), plus all the detritus knobs. Ecologically this
propagates the well-known uncertainty in transfer efficiency into an uncertainty band on PPR.

---

## 5. From SPPR to ecosystem footprint

### `get_PPR(sppr, only_inner=False, only_pp=False)`
```python
get_PPR(sppr, only_inner=False, only_pp=False) -> pd.DataFrame | pd.Series
```
Converts a **per-group SPPR** into the **total primary production required by the catch**:
`PPR = catch · SPPR`. Each group's SPPR is weighted by how much of it we actually harvest and
summed. The input SPPR is relabelled to seq, reindexed onto the catch, and infinities zeroed.

- **`only_inner`** — drop the Import columns, so only production generated *inside* the system is
  counted (exclude subsidies imported across the boundary).
- **`only_pp`** — drop **both** Import and Detritus columns, counting only genuine within-system
  **primary** production. (`only_pp` is stronger than `only_inner` and overrides it.)
- **Returns** a 1-row DataFrame if `sppr` is a DataFrame (a PPR per basal source), else a Series.
- **Ecological meaning:** this is the headline number — the tonnes of primary production the
  fishery ultimately appropriates. The `only_*` flags let you choose the accounting boundary:
  all basal sources, in-system only, or strictly primary producers.

### `get_NPP(only_inner=True)`
```python
get_NPP(only_inner=True) -> float
```
Returns the ecosystem's **net primary production** — the summed production `p` of the PP groups.
Only the within-system case is implemented (`only_inner=True`; `False` raises). This is the
total plant/algal production available as the denominator for the footprint ratio below.

### `get_PPR2NPP_ratio(sppr, only_pp=False)`
```python
get_PPR2NPP_ratio(sppr, only_pp=False) -> float
```
The **fraction of available net primary production appropriated by the catch**: within-system
`PPR / NPP`. This is the classic "%PPR" indicator (à la Pauly & Christensen) — what share of the
sea's primary production the fishery consumes.

- **`only_pp`** — if True, count only within-system primary production in the numerator (drop
  Import and Detritus columns) for a like-for-like PP/PP comparison.
- **Ecological meaning:** a low ratio means the fishery is a small draw on the ecosystem's
  productive base; a ratio approaching (or exceeding) 1 signals that the harvest is
  appropriating an unsustainable share of primary production.

---

## 6. Quick chooser

| Want… | Use |
|-------|-----|
| A one-line classic estimate | `SPPR_1986` |
| Per-group classic estimate | `SPPR_1995` / `SPPR_1995_TL_fix` |
| Explicit food-chain paths | `SPPR_EwE` (`return_paths=True`) |
| Fast EwE-equivalent | `SPPR_EwE_Ido` |
| The 2015 input–output method | `SPPR_2015` |
| General numeric solver with detritus control | `SPPR_new` |
| Exact symbolic solution / import-cost detail | `SPPR_symbolic` |
| Uncertainty bands | `monte_carlo_SPPR` |
| Total PPR of the catch | `get_PPR` |
| Share of NPP appropriated | `get_PPR2NPP_ratio` |
