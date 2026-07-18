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
- **Flow-network methods** (`SPPR_EwE`, `SPPR_EwE_Ulanowicz`, `SPPR_2015`, `SPPR_new`,
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
| `M0` | non-predatory mortality flow | production dying of causes other than predation/fishing → routed to detritus |
| `egestion` | egestion `U` | unassimilated food (faeces) → routed to detritus |
| `EE` | ecotrophic efficiency | the fraction of production **not routed to detritus** — i.e. used within the system or exported: eaten by predators, caught, accumulated as biomass, or net-migrated. `EE = (predation + catch + BA + Nm)/p = (p − M0)/p`. Its complement `M0/p` is the fraction that dies non-predatorily to detritus. |
| `GE` | gross efficiency = `p/q` | how efficiently consumed food becomes production |
| `catch` | catch `C` | the harvest we are attributing a PPR cost to |
| `BA` | biomass accumulation (`growth`) | net change in standing biomass over the period |
| `Nm` | net migration | net biomass gained by immigration − emigration |
| `DC` | diet composition | `DC[i,j]` = fraction of predator *i*'s diet made up of prey *j* (consumer rows sum to 1; PP/Import rows are all zeros) |
| `Z` | biotic transaction matrix | `Z[i,j] = q_i · DC[i,j]`, absolute predatory flow from prey *j* into predator *i* |

The two Ecopath mass-balance identities tie these together, per group:

$$ q = \text{respiration} + U + p \qquad\qquad p = BA + \text{predation} + C + N_m + M0 $$

so that `EE·p = p − M0 = predation + C + BA + Nm` is exactly the production that leaves a
group through a route *other* than non-predatory death. This is the quantity `EE` measures, and
it is what makes `TE = GE·EE` the *fraction transferred up the web* (see below).

### Trophic categories (`trophic_info`)
Groups are classified into four kinds, each with its own `get_*_seq()` accessor:

- **PP** (primary producers) — basal autotrophs; the "source" the whole calculation traces back to.
- **DET** (detritus) — dead organic matter pools; recycled basal sources fed by `M0` + `egestion`.
- **Import** — external production entering the system (diet imported from outside the model boundary).
- **Regular** — consumers (heterotrophs).

PP, DET and Import are **basal sources**: SPPR is expressed as production required *from these*.

### `get_TE` — transfer efficiency (`TE_option`)
Transfer efficiency is the fraction of production at one level that becomes production at the
next. It is the single most important ecological assumption in every flow-network method.
`get_TE` builds a per-group TE from one of four options:

| `TE_option` | Formula | Ecological meaning |
|-------------|---------|--------------------|
| `'GE'` | `p/q` | **Gross growth efficiency** — production per unit consumed. Ignores that not all production is passed on. |
| `'TE'` | `(p/q)·(1 − M0/p)` | Gross efficiency times the **ecotrophic fraction**: only the part of production that is actually consumed/exported/accumulated/migrated counts as "transferred". Production that dies naturally (`M0`) is treated as lost from the up-web pathway. Equivalent to `GE·EE`. |
| `'With Egestion'` | `(p/q)·(q/(q−egestion))` | Efficiency computed on **assimilated** intake rather than gross intake, i.e. feces are removed from the denominator so they are accounted separately (they flow to detritus, not up the chain). |
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
Trophic level follows the standard recursive definition: a group's TL is one more than the
diet-weighted mean TL of its prey,

$$ TL_i = 1 + \sum_j DC_{ij} \cdot TL_j \quad\Longleftrightarrow\quad \mathbf{TL} = \mathbf{1} + \mathbf{DC} \cdot \mathbf{TL} \quad\Longleftrightarrow\quad \mathbf{TL} = (I - DC)^{-1}\mathbf{1} $$

Basal sources (with an all-zero diet row) sit at TL 1; the matrix inverse resolves the
recursive dependencies (including cycles) in one shot.

- **`break_cycles`**: remove cycles (Ulanowicz algorithm) before inverting, so recursive loops
  (e.g. detritus ↔ bacteria) don't distort the inversion.
- **`DET_as_PP`** — sets the **base** detritus diet row, which `TE_option` (below) may then
  overwrite:
  - `DET_as_PP=True` — detritus has an all-zero diet row (it "eats" nothing), so `TL_DET = 1`.
    This is the **Christensen & Pauly (1995) convention** — zero the detritus rows so detritus
    sits at the base — and it is what the code actually uses: `SPPR_1986`/`SPPR_1995`/
    `SPPR_1995_TL_fix` all call `get_TL(break_cycles=True, DET_as_PP=True)`.
  - `DET_as_PP=False` — the detritus row is rebuilt from the groups whose death fed it, so
    detritus rises above TL 1 (single-detritus row `= (M0 + egestion)/flow2det`).
- **`TE_option`** — rewrites the detritus diet row before the inversion; since detritus feeds
  many consumers, this shifts every TL above it. The rows use two flow quantities:
  `flow2det = Σ_k (M0_k + egestion_k)` (total dead matter reaching detritus) and, per detritus
  pool, its own inflow `q_DET` routed by the `det_fate` fractions `fracs` (with one detritus pool,
  `fracs = 1` and `q_DET = flow2det`).
  - `'TE'` — zero the detritus→consumer entries, leaving detritus with only its primary-producer
    ancestry (no consumer-derived contribution). Combined with `DET_as_PP=True` the row is already
    zero, so detritus stays at TL 1.
  - `'GE'` — rebuild the detritus row from **mortality** provenance: `DC[DET, :] = (M0·fracs)/q_DET`,
    i.e. single-detritus `DC[DET, k] = M0_k/flow2det`. Detritus "eats" each group in proportion to
    its non-predatory mortality, and inherits a fractional TL just above the mean TL of that dead
    matter.
  - `'With Egestion'` (**default**) — leave the row as `DET_as_PP` built it. With `DET_as_PP=False`
    (single detritus) that row is `(M0_k + egestion_k)/flow2det`, so both mortality and egested
    material carry trophic ancestry into detritus.

### `DET_as_PP` and `normalize` (in `get_DC` / `get_Z`)
`DET_as_PP` controls whether **detritus is a source or a recycling loop**:
- `DET_as_PP=True` — detritus is a basal source: its diet row is the stored `DC` (in practice
  terminal / zero outgoing diet in `get_Z`), so no flow is traced *out* of detritus. This is what
  the flow-network SPPR methods use, so that detritus becomes one of the basal-source "columns"
  you get an SPPR breakdown against.
- `DET_as_PP=False` — the detritus row is rebuilt from where its material actually came from
  (`M0` + egestion shares, routed by `det_fate`), i.e. detritus recycles back into the web; the
  non-detritus rows are rescaled to preserve their original row sums.

`normalize` controls **whether each diet row is forced to sum to 1**:
- `normalize=False` (default) — rows keep their *original* sums. A consumer whose reported diet
  does not sum exactly to 1 (rounding, or a diet-import fraction stripped out) keeps that sum, so
  `Z = DC·q` and `A = DC/TE` preserve the true absolute flows. The SPPR solvers rely on this,
  because they build `A` from the real diet fractions rather than renormalized ones.
- `normalize=True` — every row is divided by its own sum (`DC = DC / DC.sum(axis=1)`), forcing
  each consumer's diet fractions to add to exactly 1. Use this when a downstream calculation
  needs a strict probability distribution per row (e.g. a clean `(I − DC)⁻¹` TL inversion) and
  the small deviations from unity would otherwise bias the result. It is *not* used by the SPPR
  flow solvers, precisely because renormalizing would silently rescale genuine mass flows.

---

## 3. Trophic-level SPPR methods

### `SPPR_1986()` — Pauly & Christensen (1986)
```python
SPPR_1986() -> pd.DataFrame   # single 'sppr' column
```
Computes **one** catch-weighted mean trophic level for the entire catch, then
`SPPR = TE^(1−TL)` with a **fixed TE = 0.1**. Every group receives the same number.

**Equation.** With a single catch-weighted mean trophic level
$\overline{TL} = \big(\sum_i C_i TL_i\big)/\sum_i C_i$ and fixed `TE = 0.1`,

$$ \mathrm{SPPR}_i = TE^{\,1-\overline{TL}} = \left(\tfrac{1}{TE}\right)^{\overline{TL}-1} = 10^{\overline{TL}-1}\quad\text{(same for every }i). $$

This is the original Pauly & Christensen pyramid: each trophic step multiplies the requirement
by `1/TE = 10`, so a mean-TL-3.5 catch costs `10^{2.5}` units of primary production per unit
caught.

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

**Equation.** Using each group's own fractional trophic level `TL_i` (from `get_TL`) and one
global `TE`,

$$ \mathrm{SPPR}_i = TE^{1-TL_i} = \left(\tfrac{1}{TE}\right)^{TL_i-1}. $$

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
directly, it **linearly interpolates between the two bracketing integer trophic levels**.

**Equation.** Writing `TL_i = n + f` with integer part `n = ⌊TL_i⌋` and fraction
`f = TL_i mod 1`,

$$ \mathrm{SPPR}_i = (1-f)\left(\tfrac{1}{TE}\right)^{n-1} + f\left(\tfrac{1}{TE}\right)^{n} $$

So a group at `TL = 3.4` is scored as `0.6` of a pure TL-3 feeder plus `0.4` of a pure TL-4
feeder, instead of `(1/TE)^{2.4}`. The two agree at integer TL but differ in between, because
`x^{TL}` is convex — the direct fractional exponent sits *below* the straight-line blend.

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

### Mathematical foundation (shared by all flow-network methods)

**The per-edge weight matrix `A`.** Every flow-network method starts from the
**production-normalized transaction matrix**

$$ A_{ik} = \frac{DC_{ik}}{TE_i} \qquad\text{(equivalently } A_{ik}=\tfrac{Z_{ik}}{P_i\cdot EE_i}=\tfrac{Z_{ik}}{P_i-M0_i}\text{).} $$

`A_{ik}` is the number of units of prey *k*'s production directly required to make `DC_{ik}`
units of consumer *i*'s production: the diet fraction `DC_{ik}` says how much of *i*'s intake is
*k*, and dividing by the transfer efficiency `TE_i` converts "intake" into "production required"
(you need `1/TE` units in for one unit out). The `TE_option` (§2) chooses which efficiency sits
in that denominator.

**Why `A·sppr = sppr` for the sppr vector.** Let `x_i = SPPR_i` be the source requirement of one unit of group *i*. If the
system is closed with respect to the chosen basal sources, then the requirement of *i* is just
the sum of the requirements of everything it directly needs:

$$ x_i = \sum_k A_{ik}\cdot x_k \qquad\Longleftrightarrow\qquad \mathbf{x} = A\cdot\mathbf{x} \qquad\Longleftrightarrow\qquad (A - I)\cdot\mathbf{x} = 0. $$

So the SPPR vector is a fixed point of `A` — an eigenvector with eigenvalue 1 — i.e. it lives in
the **nullspace of `L = A − I`**. Intuitively, `x = Ax` says "the cost of a group equals the
summed cost of its diet"; solving it in closed form automatically sums **all** pathways and
**all** cycles (the geometric series `I + A + A² + …`), which is exactly what path enumeration
struggles with.

**Terminal sources and the basis columns.** Basal sources (PP, detritus, imports) have no diet
of their own to trace, so their rows of `A` are replaced by identity rows
(`Ã_{ii}=1`, `Ã_{ij}=0`). The nullspace of `I − Ã` then has one basis vector per terminal
source, normalized (via RREF) so that basis vector `s^b` carries a `1` on source `b` and `0` on
the other sources. Column `b` of the returned SPPR matrix is `s^b`: the units of source `b`
required, directly and indirectly, to produce one unit of each group. Summing the columns gives
the aggregated SPPR.

**The balance identity (used to sanity-check every method).** Under mass balance, primary
production entering the system equals the production leaving it, both measured in SPPR units:

$$ \text{Inflow} = NPP + \text{DietImport} = (\mathbf{N_m} + \mathbf{C} + \mathbf{BA})\cdot\mathbf{SPPR} = \text{Outflow}, $$

which is why `get_PPR2NPP_ratio` computes `C·SPPR / [(Nm+C+BA)·SPPR]`.

### `SPPR_EwE(TE_option, use_EE=True, return_paths=True, silent=True)`
```python
SPPR_EwE(TE_option, use_EE=True, return_paths=True, silent=True)
    -> (SPPR, A, paths_dict)
```
The classic **path-enumeration** approach of Ecopath with Ecosim. It enumerates *every simple
path* from each group down to a basal terminal node and sums the product of the per-edge weights
`A = DC/TE` along each path.

**Equation.** For focal group *x*,

$$ \mathrm{SPPR}_x = EE_x \!\! \sum_{\mathrm{path}\in\mathcal{P}_x}\ \prod_{(\mathrm{pred},\mathrm{prey})\in\mathrm{path}} \frac{DC_{\mathrm{pred},\mathrm{prey}}}{TE_{\mathrm{pred}}}, $$

where `P_x` is the set of **simple** paths (no node repeated) from *x* down to a basal terminal
(PP or detritus). Because paths are simple, cannibalism and cycles are **not** expanded into
repeated loops — this is exactly what the matrix methods below fix.

- **`TE_option`** — `'GE'`, `'TE'`, `'With Egestion'`, or `'global'`. Governs the per-edge weight
  `DC/TE`, i.e. how much basal production each trophic link implies (see §2). The dominant
  ecological assumption. EwE's software uses `'TE'`.
- **`use_EE`** (default True) — this is the **leading `EE_x` factor** in the equation above, and
  it multiplies **the whole row of the focal group `x`** (the group whose SPPR is being computed):
  in code, `SPPR = SPPR.mul(EE, axis='index')`, i.e. row *i* is scaled by `EE_i`.

  Why only the focal group, and why it matters: with `TE = GE·EE`, every edge weight `DC/TE`
  already contains a `1/EE_pred`. Tracing a chain backward from *x*, the **first** step uses
  `1/EE_x`; the leading `EE_x` **cancels it**, so the focal group's own step effectively uses
  gross efficiency `GE_x`, while every downstream predator keeps the full `TE = GE·EE`.
  Ecologically: a group's production splits into the ecotrophically-used
  part (eaten/caught/accumulated/migrated) and the part that dies to `M0` and drops to detritus. The
  leading `EE_x` charges the harvest only for the focal group's **useful** production, excluding
  the M0-to-detritus fraction — whereas for the intermediate predators along the chain, the
  `M0` loss *is* counted, because supporting them required feeding the fraction that later died.
  Setting `use_EE=False` drops this factor, charging the focal group's full production
  (including its `M0`) to the requirement.
- **`return_paths`** (default True) — if True, use the slower implementation that also returns
  the explicit list of food-chain paths (useful for tracing *which* chains dominate the PPR); if
  False, use the fast vectorized implementation (returns an empty paths dict).
- **`silent`** — suppress progress bars.
- **Returns** `(SPPR, A, paths_dict)`.
- **Ecological meaning:** the most literal reading of "trace the energy back". Its weakness is
  cycles: recycling loops create infinitely many paths, so only simple paths are kept (and a
  depth safety-valve caps enumeration). This *undercounts* cyclic contributions — a cannibal
  group at `TE=0.1` eating 99% PP + 1% itself gets `SPPR = 0.99/TE = 9.9` from the one simple
  path, whereas the true cycle-summed value is `(0.99/TE)/(1 − 0.01/TE) = 11`. This is why the
  matrix reformulations below exist.

### `SPPR_EwE_Ulanowicz(TE_option, global_TE='mean', use_EE=True)`
```python
SPPR_EwE_Ulanowicz(TE_option, global_TE='mean', use_EE=True) -> (SPPR, A, L)
```
A **matrix (nullspace) reformulation** in the spirit of `SPPR_EwE`. It builds `A = DC/TE`,
**removes cycles first**, replaces each basal (producer) row with an identity row, and finds the
SPPR as the **nullspace of `L = A − I`** (i.e. `A·x = x`), RREF-normalized so each output column
is anchored to one basal source.

The cycle handling is the crucial subtlety. It computes `DCNoCyc = remove_cycles(DC)` (Ulanowicz
weakest-link removal) and then **zeros in the original `DC` exactly the edges that
`remove_cycles` drove to zero** (`DC[DCNoCyc == 0] = 0`) — *without renormalizing the surviving
diet fractions*. It then solves the nullspace on that pruned matrix.

- **`TE_option`** — as above.
- **`global_TE`** — only used when `TE_option='global'`; `'mean'` or a literal float.
- **`use_EE`** — scale each focal group's row by its `EE` (same meaning as in `SPPR_EwE`).
- **Returns** `(SPPR, A, L)`.

> **`SPPR_EwE_Ulanowicz` is *not* mathematically equivalent to `SPPR_EwE` in general — they agree only
> when the food web has no cycles.**
>
> **Why.** When the graph is acyclic, `remove_cycles` is a no-op, so both methods operate on the
> same `A`; and on an acyclic graph the nullspace sum equals the simple-path sum — hence
> identical. When cycles exist the two diverge because they handle them differently and neither
> equals the true cycle-summed (Leontief) answer:
> - `SPPR_EwE` keeps **all simple paths at their full diet weights**, but drops any path that
>   would repeat a node (so it undercounts cyclic recycling).
> - `SPPR_EwE_Ulanowicz` **deletes the weakest-link edges** of each cycle and keeps the survivors at
>   full weight *without renormalizing*, so the pruned consumer's diet no longer sums to 1 — it
>   loses mass, and its nullspace SPPR is generally *lower* than the simple-path sum.
>
> So this method is best read as a fast, cycle-pruned matrix cousin of `SPPR_EwE`, **not** as the
> exact all-cycles nullspace. For that, use `SPPR_new` / `SPPR_2015`, which count cycles fully.

### `SPPR_2015()` — the 2015 method
```python
SPPR_2015() -> (SPPR, A, L)
```
A **Leontief input–output** formulation (matrix inversion). Detritus columns are *dissolved*:
the PP-derived fraction of each detritus flow is reassigned back onto the PP groups, leaving
only living compartments. The production-normalized transaction matrix `A` then yields the
production-requirement matrix `L = (I − A)⁻¹`, whose PP columns give the per-group SPPR. A
balancing detritus SPPR is added back at the end.

**Equation.** In this method the code builds `A` directly from flows and living production,

$$ A_{ij} = \frac{Z_{ij}}{P_i}, \qquad L = (I - A)^{-1}, \qquad \mathrm{SPPR}_i = L_{i,\mathrm{PP}}, $$

where `P_i` here is the **living/useful production** `P_i − M0_i` (production net of
non-predatory mortality). The code computes this as `P[non_PP] = p·EE` for consumers and then
`A = Zᵀ / P`.

**Why `A = Z/P` is the same matrix as `A = DC/TE`.** These two constructions are algebraically
identical, which is what makes `SPPR_2015` the same underlying calculation as the other
flow-network methods. Substituting the definitions `Z_{ij} = q_i·DC_{ij}`,
`P_i − M0_i = p_i·EE_i`, `GE_i = p_i/q_i` and `TE_i = GE_i·EE_i`:

$$ A_{ij} = \frac{Z_{ij}}{P_i - M0_i} = \frac{q_i\cdot DC_{ij}}{p_i\cdot EE_i} = \frac{DC_{ij}}{(p_i/q_i)\cdot EE_i} = \frac{DC_{ij}}{GE_i\cdot EE_i} = \frac{DC_{ij}}{TE_i}. $$

So dividing the absolute prey flow by the group's *useful* production is exactly the same as
dividing the diet fraction by the ecotrophic transfer efficiency `TE = GE·EE`. The powers of `A`
then enumerate paths: `[A^k]_{i,\mathrm{PP}}` is the PP required through all length-`k` chains,
and `L = (I−A)⁻¹ = I + A + A² + …` sums every path length — **cycles included** — provided the
spectral radius of `A` is `< 1` (the standard Leontief convergence condition).

- **Returns** `(SPPR, A, L)`.
- **Ecological meaning:** treats the ecosystem like an economy where each group's production
  "requires" inputs from the groups it eats; `(I − A)⁻¹` sums the full direct + indirect
  requirement chain. Consistent, closed-form, and the reference implementation of the 2015
  paper

### `SPPR_new(...)` — primary numeric solver
```python
SPPR_new(TE=None, TE_option='GE', DET_TE_vals=1,
         det_collapse_mode='never', det_open_mode='none',
         det_theta=1.0, det_external_sppr=0.0, fix_EE_0_cases=True)
    -> (SPPR, A, L)
```
The main, most general numeric SPPR method. It builds `A = DC/TE` (with detritus treated as a
basal source, `DET_as_PP=True`), replaces basal rows with identity rows, and solves the
**nullspace of `L = A − I`** (the shared foundation above), RREF-normalized to one column per
basal source. Unlike `SPPR_EwE_Ulanowicz` it does **not** prune cycles, so the living-network solution
counts all cycles exactly. Its novelty is **explicit, tunable detritus handling** — how recycled
dead organic matter is credited as a basal source.

#### The core calculation, step by step (single detritus, `det_collapse_mode='never'`)

This is the simplest and default case; start here. Take a model with one detritus pool `DET`.

1. **Solve the living network as if detritus were a free basal source.** The nullspace of
   `L = A − I` gives, for every group *i*, a raw SPPR split into one column per basal source. Two
   kinds of column matter here: the primary-producer / import columns, and the detritus column.
   Write `nonDET_sppr_i` for the summed PP+Import part and `basis_i[DET]` for the raw detritus
   column (the units of detritus required per unit of *i*, treating one unit of detritus as
   "1" for now). At this stage detritus is just another source; one unit of it is worth one unit
   of itself.

2. **Find what one unit of detritus is actually worth in primary-production units.** Detritus is
   not a true primary source — its value is the flow-weighted average SPPR of everything dying
   into it: 

   $$SPPR_{DET} = \frac{\sum_{k} F_{k \to DET} \cdot SPPR_{k}}{\sum_{k} F_{k \to DET}}$$

   Let `q_DET` be the total inflow to the pool and, for `TE_option='GE'`, define the per-group inflow share

$$ m_k = \frac{M0_k}{q_{DET}} \quad\text{(the fraction of the detritus pool supplied by group }k\text{'s non-predatory death).} $$

   By definition, one unit of detritus is worth the inflow-weighted sum of the SPPR of the material
   flowing into it, using the per-group inflow shares `m_k`:

$$ \mathrm{sppr\_det} = \sum_k m_k\cdot\mathrm{SPPR}_k. $$

   A subtlety on those weights: `q_DET` is the pool's *total* inflow — `(M0 + egestion)` summed —
   but the GE weight `m_k = M0_k/q_DET` counts only the non-predatory-mortality inflow. So the `m_k`
   sum to `ΣM0/q_DET ≤ 1`, **not** to 1: the egestion share of the inflow carries no SPPR under the
   GE convention, so it dilutes `q_DET` without contributing to the sum. (Under `'With Egestion'`
   the egestion inflow is instead credited with the SPPR of the food it came from — the extra `DCᵀ`
   term noted below — which restores the effective weights to ≈ 1)

   The catch is that each `SPPR_k` on the right is not a fixed number: from step 1 it splits into a
   **known** primary-producer/import part and an **unknown** detritus part,

$$ \mathrm{SPPR}_k = \mathrm{nonDET\_sppr}_k + \mathrm{sppr\_det}\cdot\mathrm{basis}_k[DET], $$

   whose only unknown is `sppr_det` itself — the value of one detritus unit we are solving for. So
   detritus's value depends on the groups dying into it, but their SPPRs depend on detritus's value
   (they ate detritus while alive). Substituting this split into the average and factoring the
   common `sppr_det` out of the detritus term collapses all the groups `k` into one scalar equation
   in which `sppr_det` appears on **both** sides — a self-consistency (fixed-point) equation:

$$ \mathrm{sppr\_det} = \underbrace{\sum_k m_k\cdot\mathrm{nonDET\_sppr}_k}_{a\ \text{(PP-origin material entering DET)}} + \underbrace{\Big(\sum_k m_k\cdot\mathrm{basis}_k[DET]\Big)}_{b\ \text{(recycled DET-origin material)}}\cdot \mathrm{sppr\_det}. $$

   This is exactly the cannibal-cycle logic, now applied to the whole
   detritus pool: detritus feeds consumers, whose death feeds detritus again, so its value
   depends on itself. Solving the scalar fixed point,

$$ \boxed{\ \mathrm{sppr\_det} = \dfrac{a}{1-b}\ } \qquad (b<1\text{ required for a finite, positive value).} $$

3. **Scale the detritus column** of the SPPR matrix by `sppr_det` and add it to the PP/Import
   columns. `'never'` means step 2 is always solved directly (never pooled), so if `b ≥ 1`
   (recycling so strong it diverges) the result may go negative rather than raise.

   The other two `TE_option`s change what seeds the detritus value.

   **`'With Egestion'`.** Now *both* routes into detritus carry SPPR, so the weight on group `k`
   gains a second term:

$$ m_k = \underbrace{\frac{M0_k}{q_{DET}}}_{k\text{'s own dead body}} \;+\; \underbrace{\sum_j DC_{jk}\cdot\frac{\mathrm{egestion}_j}{q_{DET}}}_{\text{prey }k\text{ egested undigested by consumers }j} \;=\; \frac{M0_k}{q_{DET}} + \Big(DC^{\mathsf{T}}\big(\mathrm{egestion}/q_{DET}\big)\Big)_k. $$

   The first term is the mortality route from GE: `k`'s carcass carries `k`'s own SPPR. The second
   term handles faeces, which are *not* the egesting consumer's production — they are food that
   passed through `j` undigested — so they must **not** be charged `SPPR_j`. Instead the diet row
   `DC[j, ·]` says what `j` ate: a fraction `DC[j,k]` of `j`'s intake (and hence of `j`'s faeces)
   was prey `k`, so that faecal flow carries `SPPR_k`. Summing each egesting consumer `j`'s faeces
   over the prey that composed them is exactly the matrix–vector product `DCᵀ · (egestion/q_DET)`:
   the transpose "un-mixes" every consumer's faeces back into its prey species and credits detritus
   the *prey's* SPPR, not the consumer's — matching the physical fact that faeces are undigested
   food. As a bonus this closes the accounting gap noted above: the consumer diet rows sum to 1, so
   the egestion term sums to `Σegestion/q_DET`; added to the mortality term's `ΣM0/q_DET` the
   weights now total `(ΣM0 + Σegestion)/q_DET = 1`, the fully normalized flow-weighted average that
   GE was missing.

   **`'TE'`.** Here the per-edge efficiency already folds in the ecotrophic factor `EE = (p−M0)/p`,
   so production lost to non-predatory death `M0` is treated as *gone* — it carries no SPPR onward.
   Dead consumer bodies therefore do **not** re-seed detritus with recycled SPPR: there is no
   detritus→consumer→detritus loop, hence **no recycling matrix and no fixed point** to solve. Each
   detritus column is simply scaled by the share of the pool's inflow arriving *directly* from
   primary producers and imports:

$$ \mathrm{sppr_{det}} = \theta\cdot\frac{\sum_{k\in PP\cup Import} F_{k\to DET}}{q_{DET}}, $$

   i.e. detritus is credited only with the genuinely primary material that fell into it, while all
   consumer-derived inflow is uncounted because the `'TE'` convention has already written it off as
   lost. Note that in this convention, `SPPR_det` is always $\leq$ 1. Here `θ = det_theta` is the availability/retention damping. (The one exception is `EE=0`
   dead-end groups, whose `TE=0` severs them from the nullspace and leaks the PP they consumed;
   `fix_EE_0_cases` re-credits that leak to detritus with a small linear correction.)

#### The coupled recycling system `(I − B)x = c` (multiple detritus pools)

With more than one detritus pool the single scalar becomes a **vector** `x = (sppr_det_1, …,
sppr_det_k)`, because pools feed each other: a consumer eating pool *j* can die into pool *l*, so
pool *l*'s value depends on pool *j*'s value. Repeating step 2 per pool `l`:

$$ x_l = \underbrace{\sum_k m^{(l)}_k\cdot\mathrm{nonDET\_sppr}_k}_{c_l} + \sum_{j} \underbrace{\Big(\sum_k m^{(l)}_k\cdot\mathrm{basis}_k[DET_j]\Big)}_{B_{lj}}\cdot x_j, $$

which in matrix form is the linear system `_build_det_BC` assembles and `_solve_det_scaling`
solves:

$$ \mathbf{x} = \mathbf{c} + B\cdot\mathbf{x} \qquad\Longleftrightarrow\qquad (I - B)\cdot\mathbf{x} = \mathbf{c}. $$

The inflow share `m^{(l)}_k` — the fraction of pool `l`'s inflow supplied by group `k`, routed by
`det_fate` — is the multi-pool version of the single-detritus weight from step 2:

$$ m^{(l)}_k = \frac{M0_k\cdot\mathrm{fracs}^{(l)}_k}{q_l} \;+\; \underbrace{\Big(DC^{\mathsf{T}}\big(\mathrm{egestion}\cdot\mathrm{fracs}^{(l)}/q_l\big)\Big)_k}_{\text{egestion route ('With Egestion' only)}}, $$

where `q_l` is pool `l`'s total inflow and the routing weight `fracs^{(l)}_k = det_fate[k, l]` is the
fraction of group `k`'s flow-to-detritus that reaches pool `l` (with a single pool `fracs^{(l)} = 1`,
recovering the scalar `m_k`). The egestion term is present only for `TE_option='With Egestion'`;
`'GE'` keeps just the mortality term. Reading the remaining pieces (all defined per pool `l`):

- **`c_l`** = the primary-production-origin SPPR entering pool `l` (mortality/egestion of PP and
  of the PP-derived part of consumers) — the "new" material.
- **`B_{lj}`** = how much of pool `l`'s value comes from pool `j`'s value, via consumers that eat
  pool `j` and then die into pool `l`. This is the **recycling coupling**, and it is where `DC`
  and `A` enter: `basis_k[DET_j]` comes straight from the living-network nullspace of `A = DC/TE`
  (how much of pool `j` each group needs), and `m^{(l)}_k` comes from the mortality/egestion
  flows into pool `l`. So `B` is the composition "trace pool `j` up through the diet
  (`A`/nullspace), then back down into pool `l` through death (`M0`, egestion)".
- **`(I − B)⁻¹ = I + B + B² + …`** is again a Leontief sum: `I` is the direct PP-origin input,
  `B` is one recycling loop through the detritus system, `B²` two loops, and so on. It converges
  when the spectral radius of `B` is `< 1` (recycling loses mass each pass). The single-detritus
  `a/(1−b)` above is exactly this with a 1×1 `B = [b]` and `c = [a]`.

So the detritus resolution is a **second linear solve layered on top of the living-network
nullspace**: the nullspace (from `A = DC/TE`) fixes how much of each detritus pool every group
needs; the `(I − B)x = c` system then converts those pools from "one unit of themselves" into
primary-production-equivalent values by closing the death→detritus→consumption→death loop.

**Core inputs:**
- **`TE`** — supply an explicit TE matrix (e.g. a Monte-Carlo sample). If `None`, built from
  `TE_option`.
- **`TE_option`** — `'GE'` (default), `'TE'`, `'With Egestion'`, `'global'` (see §2). Beyond
  setting the per-edge weights `A = DC/TE`, it also selects the detritus resolution described
  above: `'TE'` scales each detritus column by its direct PP+Import inflow share (no recycling
  matrix), while `'GE'` / `'With Egestion'` build and solve the coupled recycling system
  `(I − B)x = c`.
- **`DET_TE_vals`** (default 1) — TE for detritus rows when building the TE matrix.

**Detritus recycling knobs** (the ecological heart of the method): these all act on the recycling
system from the previous section — the single-pool fixed point `sppr_det = a/(1−b)` and its
multi-pool form `(I − B)x = c`, where `b` (the spectral radius of `B`) measures recycling
strength: how much detritus one unit of detritus regenerates via the death→detritus→consumption→death
loop.

- **`det_collapse_mode`** — what to do when that loop is **too strong to invert**. A finite,
  non-negative solution exists only while `b < 1`; once `b ≥ 1` each unit of detritus regenerates
  ≥ 1 unit and the closed system diverges (`1 − b ≤ 0`). This flag chooses the response:
  - `'never'` (default): always solve the coupled system directly, even at `b ≥ 1`. The answer may
    come out **negative** — a clear signal that recycling diverged — but the call never raises.
    This is the default because the Monte-Carlo samplers solve thousands of TE draws and simply
    discard the non-physical (negative) ones.
  - `'auto'`: solve directly **unless** the system is unstable (`b ≥ 1`, or `(I − B)`
    ill-conditioned), in which case fall back to pooling. The safe, self-correcting choice.
  - `'always'`: always pool.
  - *What "pooling" does, ecologically:* it merges all detritus pools into **one** combined
    compartment and solves a single scalar `sppr_det = a/(1−b)`. Because the combined pool's inflow
    `q_combined` is the sum of every pool's inflow, it is far larger than any single `q_l`, and the
    bigger denominator drives `b` back below 1 — taming an otherwise runaway loop at the cost of
    resolution (every detritus pool then shares one blended value).

- **`det_open_mode`** — whether detritus recycling is **closed** (all dead matter is reprocessed
  inside the system) or **open** (some escapes, or some is supplied from outside). It works by
  transforming `(B, c)` *before* the solve, using the retention fraction `θ = det_theta` and the
  external value `ext = det_external_sppr` (both defined below):
  - `'none'` (default): closed recycling — `B` and `c` are used unchanged, so the equations above
    are solved as-is.
  - `'recycling_loss'`: a fraction `1 − θ` of the dead matter entering detritus is permanently lost
    (buried in sediment, exported off the shelf) instead of recycled. Only the **recycling term** is
    damped: `B → diag(θ)·B`, `c` unchanged, so the scalar becomes `sppr_det = a/(1 − θ·b)`. Each pass
    through detritus loses mass, shrinking the recycled contribution — and rescuing divergent
    (`b ≥ 1`) cases, since `θ·b` can fall below 1.
  - `'source_dilution'`: the same damping **plus** the lost fraction of the inflow is replaced by
    detritus supplied from *outside* the model, carrying a fixed SPPR `ext`. Now **both** terms
    change: `B → diag(θ)·B` and `c → θ·c + ext·(1 − θ)`, giving
    `sppr_det = (θ·a + ext·(1 − θ)) / (1 − θ·b)`. Ecologically, detritus becomes a partly-subsidized
    basal source: a blend of internally-recycled material (weight `θ`) and imported dead matter of
    value `ext` (weight `1 − θ`).

- **`det_theta`** (`θ`, default 1.0) — the detritus **availability / retention fraction**: the share
  of dead matter actually recycled within the system. `1.0` recovers the closed system (making
  `det_open_mode` a no-op); lower values mean more is buried/exported/lost. It is exactly the
  `diag(θ)` in the transforms above. A float (all pools) or a dict keyed by DET seq or name (a
  per-pool `θ`).
- **`det_external_sppr`** (`ext`, default 0.0) — the SPPR charged to the externally-supplied
  detritus under `'source_dilution'` (the `ext` in `c → θ·c + ext·(1 − θ)`). `0.0` treats the
  imported dead matter as free primary production; a positive value gives it a cost. Same float-or-
  dict form as `det_theta`.

- **Returns** `(SPPR, A, L)`.

**The `fix_EE_0_cases` correction (default True).** This is *not* a recycling knob but a
mass-balance fix specific to `TE_option='TE'`. Groups with `EE=0` (all their production dies
naturally, `M0=p`) get transfer efficiency 0 and are severed from the nullspace, which **leaks the
primary production they consumed** — it goes neither up the web nor back to detritus, so the global
PP balance breaks. When True, that consumed PP is **re-credited to the detritus pool**, since a
group whose production is 100% non-predatory mortality physically flows to detritus; this restores
inflow = outflow. It is only active for **single-detritus** models under `'TE'` (a no-op otherwise),
and emits a `RuntimeWarning` whenever `EE=0` groups are present. It does *not* fix near-singular
`0 < EE ≪ 1` groups, whose `SPPR ~ 1/te` blows up — an inherent singularity of the TE method that
raises a separate warning.

### `SPPR_symbolic(...)` — symbolic solver
```python
SPPR_symbolic(TE=None, TE_option='GE', diet_import_option='as_DC', DET_TE_vals=1,
              sppr_det_value=None, det_collapse_mode='never',
              det_open_mode='none', det_theta=1.0, det_external_sppr=0.0,
              fix_EE_0_cases=True) -> (sppr_symbolic, sppr_mat, equations, variables)
```
A **symbolic (SymPy)** counterpart to `SPPR_new`. It writes the per-group SPPR balance
equations `A·x − x = 0` symbolically and solves them exactly, returning both the symbolic
solution and a numeric basis matrix. It shares all of `SPPR_new`'s detritus knobs
(`det_collapse_mode`, `det_open_mode`, `det_theta`, `det_external_sppr`, `fix_EE_0_cases`) with
identical meaning.

Its distinctive input is **how imported diet is treated**. Imported food crosses the model
boundary, so one cannot automatically say one unit of it equals one unit of internal primary
production. The two options answer different accounting questions:
- **`diet_import_option`**:
  - `'as_DC'` (default): imported diet is kept as a **separate production source** with its own
    `DIET_SPPR`, solved from a second linear system. After the internal SPPRs are known, the
    import value for consumer *i* is inferred as the **weighted mean SPPR of that consumer's
    non-import diet**:

$$ (1-DC_{i,DI})\cdot\mathrm{DIET\_SPPR}_i = \sum_{k\in\mathcal{X}} DC_{ik}\cdot\mathrm{SPPR}_k \;\Longrightarrow\; \mathrm{DIET\_SPPR}_i = \frac{\sum_{k\in\mathcal{X}} DC_{ik}\cdot\mathrm{SPPR}_k}{\sum_{k\in\mathcal{X}} DC_{ik}}, $$

    where `𝒳` is the set of internal (regular/detritus/PP) compartments and `DI` the import node.
    So imported food is costed by what the consumer's *internal* diet is made of — the most
    faithful treatment of cross-boundary subsidies.
  - `'as_PP'`: imported diet is treated as **just another primary-production source** — its row
    of `A` becomes an identity row and its own SPPR symbol is fixed to 1, exactly like a primary
    producer. Simpler and transparent, but makes one unit of imported food directly comparable to
    one unit of internal primary production (so the import contribution should be reported
    separately).
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
Runs `SPPR_new`, then treats the detritus scaling as a single unknown `x = sppr_det` and solves
the scalar equation that forces **exact global mass balance** — PP inflow equals the export
outflow (catch + growth + net migration) weighted by SPPR:

$$ \sum_{i\in PP} p_i \;=\; \sum_i (C_i + BA_i + N_{m,i})\big(\mathrm{SPPR}^{PP}_i + x\cdot\mathrm{SPPR}^{DET}_i\big), $$

which is linear in `x`, so it is solved directly; the detritus columns are then scaled by the
solved `x`. Returns `(sppr, sppr_det)`. Ecologically, it pins the one free recycling degree of
freedom so that total primary production in equals total exported production out.

### `monte_carlo_SPPR(...)` / `monte_carlo_SPPR_2(...)`
Uncertainty propagation over transfer efficiency. Each per-group TE is resampled from a gamma
distribution centred on the model value `\overline{TE}_i`, with shape and scale set so the mean
is preserved and the coefficient of variation is `η` (= `TE_error_percent`):

$$ \widetilde{TE}_i \sim \mathrm{Gamma}(\alpha,\theta_i),\qquad \alpha=\frac{1}{\eta^2},\qquad \theta_i=\overline{TE}_i\cdot\eta^2, $$

then clipped to `[\overline{TE}_i(1-\delta),\ \overline{TE}_i(1+\delta)]` (`δ` =
`TE_error_cut_percent`). SPPR is recomputed (via `SPPR_new` or `SPPR_symbolic`) on each sampled
matrix, non-physical (negative-SPPR) draws are discarded, and the accepted draws are averaged.
This matters ecologically because `A = DC/TE` depends on `1/TE`, which is convex, so
`E[1/TE] ≥ 1/E[TE]` (Jensen): sampling *before* the nonlinear solve gives a higher, less biased
expected PPR than plugging in the mean TE. Key knobs: `n_samples`, `TE_error_percent`,
`TE_error_cut_percent`, plus all the detritus knobs.

---

## 5. From SPPR to ecosystem footprint

### `get_PPR(sppr, only_inner=False, only_pp=False)`
```python
get_PPR(sppr, only_inner=False, only_pp=False) -> pd.DataFrame | pd.Series
```
Converts a **per-group SPPR** into the **total primary production required by the catch**. Each
group's SPPR is weighted by how much of it we actually harvest and summed — an inner product
with the catch vector, per basal source *s*:

$$ PPR_s = \mathbf{C}\cdot\mathbf{SPPR}_s = \sum_i C_i\cdot(\mathrm{SPPR}_s)_i, \qquad PPR = \sum_{s\in\text{sources}} PPR_s. $$

The input SPPR is relabelled to seq, reindexed onto the catch, and infinities zeroed.

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
sea's primary production the fishery consumes. Using the balance identity (§4), the denominator
can be read as the total production leaving the system, so

$$ \frac{PPR}{NPP} = \frac{\mathbf{C}\cdot\mathbf{SPPR}}{(\mathbf{N_m}+\mathbf{C}+\mathbf{BA})\cdot\mathbf{SPPR}}. $$

With `Nm = 0` (the usual case), the ratio can exceed 1 when biomass is being depleted
(`BA·SPPR < 0`), and falls below 1 when large outflows (mostly detrital) absorb part of the
primary production.

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
| Fast cycle-pruned matrix cousin of EwE (≠ EwE when cycles exist) | `SPPR_EwE_Ulanowicz` |
| The 2015 input–output method (full cycles) | `SPPR_2015` |
| General numeric solver with detritus control | `SPPR_new` |
| Exact symbolic solution / import-cost detail | `SPPR_symbolic` |
| Uncertainty bands | `monte_carlo_SPPR` |
| Total PPR of the catch | `get_PPR` |
| Share of NPP appropriated | `get_PPR2NPP_ratio` |
