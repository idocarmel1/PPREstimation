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

- **`det_collapse_mode`** — what to do when that loop is **too strong to invert**. `B` is entrywise
  non-negative, so by Perron–Frobenius `(I − B)⁻¹ ≥ 0` **iff** `b < 1`: a finite non-negative
  solution is guaranteed for *every* admissible `c` only while `b < 1`, and at `b ≥ 1` each unit of
  detritus regenerates ≥ 1 unit. This flag chooses the response:
  - `'never'` (default): always solve the coupled system directly, even at `b ≥ 1`. The answer
    usually comes out **negative** — a clear signal that recycling diverged — and the call
    essentially never raises. This is the default because the Monte-Carlo samplers solve thousands
    of TE draws and simply discard the non-physical (negative) ones.
  - `'auto'`: solve directly **unless** the system is unstable (`b ≥ 1`, or `(I − B)`
    ill-conditioned), in which case fall back to pooling. The safe, self-correcting choice.
  - `'always'`: always pool.
  - *What "pooling" does, ecologically:* it merges all detritus pools into **one** combined
    compartment and solves a single scalar `sppr_det = a/(1−b)`. Because the combined pool's inflow
    `q_combined` is the sum of every pool's inflow, it is usually far larger than any single `q_l`,
    and the bigger denominator drives `b` back below 1 — taming an otherwise runaway loop at the
    cost of resolution (every detritus pool then shares one blended value). This is not guaranteed:
    pooling also sums the numerator, so when one pool dominates `q_combined` it can make `b`
    *worse*. 

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

**The `fix_EE_0_cases` correction (default True).** This is a
mass-balance fix specific to `TE_option='TE'`. Groups with `EE=0` (all their production dies
naturally, `M0=p`) get transfer efficiency 0 and are severed from the nullspace (Becasue the denominator in A is 0), which **leaks the primary production they consumed** — it goes neither up the web nor back to detritus, so the global PP balance breaks. When True, that consumed PP is **re-credited to the detritus pool**, since a
group whose production is 100% non-predatory mortality physically flows to detritus; this restores
inflow = outflow. It is only active for **single-detritus** models under `'TE'` (a no-op otherwise),
and emits a `RuntimeWarning` whenever `EE=0` groups are present. It does *not* fix near-singular
`0 < EE ≪ 1` groups, whose `SPPR ~ 1/te` blows up — an inherent singularity of the TE method that
raises a separate warning.

- **Returns** `(SPPR, A, L)`.

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

$$ (1-DC_{i,DI})\cdot\mathrm{DIET\_SPPR}_i = \sum_{k\in\mathcal{X}} DC_{ik}\cdot\mathrm{SPPR}_k \;\Longrightarrow\; \mathrm{DIET\_SPPR}_i = \frac{\sum_{k\in\mathcal{X}} DC_{ik}\cdot\mathrm{SPPR}_k}{\sum_{k\in\mathcal{X}} DC_{ik}} $$

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
  non-negative SPPR is guaranteed for every admissible `c`), `≥ 1` means the Neumann series
  diverges. Note the solve still returns finite numbers at `b ≥ 1`; they are the analytic
  continuation of a divergent series, not a solution.
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

### `monte_carlo_SPPR(...)` — uncertainty propagation over transfer efficiency
```python
monte_carlo_SPPR(n_samples=1000, TE_error_percent=10, TE_error_cut_percent=20,
                 TE_option='GE', DET_TE_vals=1, kind='new', method_kwargs=None,
                 exclude_diverged=False, return_diagnostics=False,
                 diet_import_option='as_DC',
                 silent=True, det_collapse_mode='never', det_open_mode='none',
                 det_theta=1.0, det_external_sppr=0.0)
    -> (mean_sppr, samples, rejection_fraction, equations, variables[, diagnostics])
```
Transfer efficiency is the single most uncertain input in every flow-network method, and SPPR
depends on it nonlinearly. Rather than solve once at the model's mean TE, this method **propagates
the TE uncertainty**: it draws many TE matrices, re-solves SPPR on each, and averages. Each
per-group TE is resampled from a **gamma** distribution centred on the model value
$\overline{TE}_i$, with shape and scale set so the mean is preserved and the coefficient of
variation equals $\eta$ (`TE_error_percent`):

$$ \widetilde{TE}_i \sim \mathrm{Gamma}(\alpha,\theta_i),\qquad \alpha=\frac{1}{\eta^2},\qquad \theta_i=\overline{TE}_i\cdot\eta^2. $$

Each draw is then clipped to $[\overline{TE}_i(1-\delta),\ \overline{TE}_i(1+\delta)]$ (with
$\delta$ = `TE_error_cut_percent`), SPPR is recomputed (via `SPPR_new` or `SPPR_symbolic`) on the
sampled matrix, non-physical (negative-SPPR) draws are discarded, and the accepted draws are
averaged. The gamma is chosen deliberately: it keeps every sampled TE **strictly positive** with
the requested mean and spread — essential because $1/TE$ diverges near zero and a negative TE is
meaningless, so a symmetric normal would be wrong here.

**Non-positive TEs are pinned, not sampled.** `get_TE` legitimately returns $\overline{TE}_i = 0$ for
some groups (and floating-point cancellation can leave such an entry a tiny negative), which the
deterministic solvers handle but gamma rejects outright — its scale $\theta_i$ must be strictly
positive, so drawing those would raise and take the whole Monte-Carlo run down. Every entry with
$\overline{TE}_i \le 0$ is therefore held fixed at its model value for all draws and excluded from
the clip band, so the solver sees exactly what the deterministic method sees for that group; only
the strictly positive TEs are resampled. When all TEs are positive this is a no-op. Note that such
groups usually make the system near-singular ($SPPR \sim 1/TE$), so draws on those models tend to be
rejected as diverged — an honest `n_rejected_diverged` rather than a crash. `diagnose_sppr` reports
them via its "near-zero TE" warning.

**Choosing the solver and passing it parameters.** `kind` selects `SPPR_new` (`'new'`) or
`SPPR_symbolic` (`'symbolic'`) — the only two solvers that accept an injected `TE` matrix, which is
what resampling requires. Anything else that solver takes goes through `method_kwargs`, e.g.
`method_kwargs={'fix_EE_0_cases': True}`. Two guards fire *before* any sampling: a key the chosen
solver does not accept raises `TypeError` naming it, and a key this wrapper sets itself (`TE`,
`TE_option`, `DET_TE_vals`, `diet_import_option`, the four `det_*` knobs) raises `ValueError`
rather than being silently overridden by one side or the other.

**Rejecting diverged draws, not just negative ones.** `exclude_diverged=True` (requires
`kind='new'`) routes each draw through `diagnose_sppr(..., return_sppr=True)` and discards it when
`report['divergence']['status'] == 'FAIL'`, i.e. when the resampled TE pushed the recycling gain
$b$ over 1. This is free: `diagnose_sppr` grades and solves in one pass, so the cost is the same
solve the method was already paying. It matters because **a negative SPPR is an unreliable symptom
of divergence** — a draw with $b>1$ can still return an all-positive SPPR and slip past the
negative-SPPR test, biasing the mean. Rejections are counted separately by cause; pass
`return_diagnostics=True` for a sixth return element holding `n_accepted`,
`n_rejected_negative`, `n_rejected_diverged`, the matching fractions, and `per_sample_reason`
(one of `'accepted'` / `'negative'` / `'diverged'` per draw, in draw order). If *every* draw is
rejected the mean is all-NaN and a `UserWarning` is emitted instead of averaging an empty stack.

**Why resample instead of solving at the mean — the Jensen point (the whole reason this exists).**
The per-edge weight is $A = DC/TE$, so SPPR is a function of $1/TE$, which is **convex**. By
Jensen's inequality, the average of a convex function exceeds the function of the average:

$$ \mathbb{E}\!\left[\tfrac{1}{TE}\right] \;\ge\; \frac{1}{\mathbb{E}[TE]} \qquad\Longrightarrow\qquad \mathbb{E}\big[\mathrm{SPPR}(TE)\big] \;\ge\; \mathrm{SPPR}\big(\mathbb{E}[TE]\big). $$

In words: solving once at the mean TE (a single deterministic call) **systematically
underestimates** the expected PPR, because it ignores the extra requirement contributed by the
low-TE tail — where each inefficient transfer costs disproportionately more basal production.
Averaging over sampled solves recovers the unbiased (higher) expectation. The gap grows with the
assumed uncertainty $\eta$ and with trophic depth (the more `1/TE` factors multiply along a chain,
the stronger the convexity). This is the ecological payoff of the method, not merely an error bar.

**Parameters:**
- **`n_samples`** (default 1000) — number of TE draws / SPPR solves averaged. Larger → tighter,
  less noisy Monte-Carlo estimate (at linear cost).
- **`TE_error_percent`** ($\eta$, default 10) — the assumed **relative uncertainty** on transfer
  efficiency: each TE's standard deviation as a percentage of its mean (the gamma's coefficient of
  variation). `10` means every TE is drawn with a 10% spread around the model value. This is the
  knob that drives the Jensen correction — set it to 0 and the method collapses to the single
  deterministic solve.
- **`TE_error_cut_percent`** ($\delta$, default 20) — a hard **clip band** of $\pm\delta\%$ applied
  to each sampled TE *after* drawing, trimming extreme tail draws that would otherwise produce
  absurd `1/TE` spikes. `20` means no sample departs more than 20% from its mean.
- **`TE_option`** / **`DET_TE_vals`** — select *which* mean TE matrix is sampled around (see §2 /
  `SPPR_new`).
- **`kind`** — which solver to resample: `'new'` (`SPPR_new`) or `'symbolic'` (`SPPR_symbolic`).
- **`diet_import_option`** — `'as_DC'` / `'as_PP'`, forwarded to `SPPR_symbolic` when
  `kind='symbolic'` (ignored for `'new'`).
- **`silent`** — suppress the progress bar / prints.
- **`det_collapse_mode`**, **`det_open_mode`**, **`det_theta`**, **`det_external_sppr`** — the
  detritus recycling knobs (see `SPPR_new`), forwarded unchanged to every SPPR solve in the loop.
- **Returns** `(mean_sppr, samples, rejection_fraction, equations, variables)`: the accepted-sample
  mean SPPR, the stacked array of accepted samples (for building confidence intervals), the fraction
  of draws discarded as non-physical, and the symbolic `equations` / `variables` (both `None` unless
  `kind='symbolic'`). A **negative SPPR** — the rejection trigger — signals the resampled TE drove
  the detritus recycling loop unstable (spectral radius `b ≥ 1`); dropping those draws keeps them
  from biasing the mean, and a high `rejection_fraction` flags a model near that instability.

---

### `diagnose_sppr(...)` — is this result trustworthy?
```python
diagnose_sppr(TE_option='GE', *, short=False, flat=False, thresholds=None,
              return_sppr=False, **sppr_kwargs)
    -> report                      # or (report, SPPR, A, L) when return_sppr=True
```
Every method above returns numbers whether or not those numbers mean anything. `diagnose_sppr`
runs `SPPR_new` once and answers the prior question: **is the solve you just made a convergent
one, on data that closes?** It grades the Ecopath input, the two convergence conditions and the
global PP budget, and reports the catch footprint without grading it.

The verdict describes a **configuration, not a model**. The recycling gain `b` is a function of
`TE_option`, `det_theta` and `det_open_mode`, so one model can be healthy under one configuration
and divergent under another; the configuration evaluated is echoed back under `config`.

#### The two convergence conditions

The detritus condition derived above (`b < 1` for `(I − B)x = c`) has a twin on the living side,
and the report tests both. Recall from *Terminal sources and the basis columns* that `SPPR_new`
solves the fixed point `x = Ã·x`, where `Ã` is the per-edge weight matrix `A_{ik} = DC_{ik}/TE_i`
with every **basal** row replaced by an identity row. Split the group index set accordingly:

- **`L`** — the **living** groups: those with a diet, i.e. the rows of `A` left untouched
  (`Regular` groups).
- **`B`** — the **basal** sources: `PP`, `Import` and — because `SPPR_new` calls
  `get_DC(DET_as_PP=True)` — the detritus pools too. These are the rows replaced by identity rows,
  and the columns of the returned SPPR matrix.

Writing the fixed point `x = Ã·x` in those blocks, with `x_L` the SPPR values of the living groups
and `x_B` those of the basal sources:

$$ \mathbf{x}_L = A_{LL}\cdot\mathbf{x}_L + A_{LB}\cdot\mathbf{x}_B, \qquad\qquad \mathbf{x}_B = \mathbf{x}_B \;\;\text{(the identity rows).} $$

- **`A_LL`** is the living→living block of `A`: `A_LL[i,k]` = units of **living** prey *k*'s
  production required per unit of living consumer *i*'s production. Its non-zeros are exactly the
  predator-eats-predator edges, so this block is where predation cycles and cannibalism live.
- **`A_LB`** is the living→basal block: how much of each basal source (PP, import, detritus) each
  living group requires **directly**, before any indirect path is traced.

Solving the first block for `x_L` gives

$$ \mathbf{x}_L = (I - A_{LL})^{-1}\cdot A_{LB}\cdot\mathbf{x}_B \;=\; \big(I + A_{LL} + A_{LL}^2 + \dots\big)\cdot A_{LB}\cdot\mathbf{x}_B, $$

the same Leontief series as `(I − B)⁻¹` in the detritus solve, but summing paths **through living
groups** instead of through detritus pools: `A_LB·x_B` is the direct basal requirement, `A_LL·A_LB·x_B`
the requirement one predation step further back, and so on. It converges iff
`ρ(A_LL) < 1`, reported as `rho_living`, while `b = ρ(diag(θ)·B) < 1` is the corresponding condition
for the detritus loop.

The two are not additive and not interchangeable. `B` is assembled **from** the living nullspace
basis — `B[l,j]` contains `basis_k[DET_j]`, which is a column of exactly the inverse above — so
`ρ(A_LL) ≥ 1` corrupts the basis and thereby invalidates `b`, rather than contributing a second,
separate failure. This is why a living-side failure is reported as `FAIL` with a note that `b` is
unreliable, and not merely added to the tally.

#### Why convergence is necessary but nowhere near sufficient

`1/(1−b)` is finite for every `b < 1`, but it is not *bounded*: as `b` approaches 1 the recycled
contribution grows without limit while every convergence test still passes. The report therefore
grades two further quantities alongside the convergence booleans — how close `b` is to 1
(`b_warn`, default `0.7`), and the magnitude of the resulting per-pool cost `max_sppr_det`
(`sppr_det_warn`, default `10`).

Ecologically, `sppr_det ≈ 1` is the healthy signature: one unit of detritus costs about one unit of
primary production, because most of it fell in directly from producers. Values up to ~10 remain
plausible for a detritus pool fed largely by consumer mortality, several trophic steps up. Beyond
that the pool's material has been round-tripped so many times that the accounting charges it tens
of times its own mass in upstream production — North Sea (1991) under `'GE'` reaches
`sppr_det = 62.6` at `b = 0.888`, which converges and is still not a usable number.

#### How the reported quantities are obtained

`b` is read from `detritus_resolution_info['rho_B']`, which `_solve_det_scaling` records **before**
the solve-vs-pool decision — so `det_collapse_mode` cannot soften the diagnosis, and only sets the
reported `would_pool`. Each pool's `sppr_det` is read straight off the result as `SPPR.loc[d, d]`:
the detritus column is pivoted at 1 on its own row before scaling (see *Terminal sources and the
basis columns*), so after scaling that entry **is** the pool's resolved SPPR. The identity holds for
every `TE_option` and for the pooled fallback alike, and reproduces `x_vec` exactly when the coupled
system was solved directly.

`n_negative_sources` counts basal-source **columns** containing a negative value, not groups: a
negative detritus column is routinely masked inside a group's row total by its positive PP columns.
`expect_negatives` is the *prediction* `b ≥ 1`; the reducible-`B` exception noted in §4 is not tested
for, so treat `expect_negatives` as the expectation and `n_negative_sources` as the fact.

#### It does not raise on the model it is diagnosing

A singular `b = 1` solve, or a pooled fallback that itself diverges, is caught and reported as
`FAIL`; a diagnostic re-solve with `det_collapse_mode='never'` then recovers `b` and `rho_living`,
so the report can still say *why* the configuration failed instead of returning nothing. In that
case `balance` and `footprint` are `None` — the requested configuration produced no SPPR to grade —
and `would_pool` is predicted from the requested mode rather than observed. A `ValueError` is still
raised for an unrecognized `TE_option`, `det_open_mode` or `det_collapse_mode`, since those are
caller mistakes rather than model pathologies.

**Parameters:**
- **`TE_option`** (default `'GE'`) — as in `SPPR_new`. It moves *both* convergence numbers, since
  both are functions of `A = DC/TE`, and it selects the detritus resolution: `'GE'` /
  `'With Egestion'` / `'global'` build the recycling system, while `'TE'` has no recycling matrix at
  all and reports `b = 0.0` with a note.
- **`short`** (default `False`) — return only `status`, `model_input`, `divergence` and `balance`,
  dropping the ungraded `footprint`, the echoed `config` and the `warnings` list.
- **`flat`** (default `False`) — return one level of `'<section>_<field>'` keys, with `sppr_det`
  expanded per detritus seq and list fields replaced by counts, so many models concatenate into a
  single `pd.DataFrame`.
- **`thresholds`** (default `None`) — per-key overrides of `DEFAULT_DIAGNOSTIC_THRESHOLDS`; missing
  keys keep their defaults. The `model_balance_*` pair was calibrated against
  `real_models/EwE_jsons/`: across 222 loadable models the `is_model_balanced` residuals are
  bimodal — the 176 balanced models sit at `≤ 1e-6` with nothing at all between `1e-6` and `1e-4` —
  so `1e-4` separates cleanly and `0.1` isolates the 32 severely broken ones.
- **`return_sppr`** (default `False`) — also return `(SPPR, A, L)` from the internal solve, so the
  sympy nullspace is not paid for twice when you want both the diagnosis and the result.
- **`**sppr_kwargs`** — forwarded verbatim to `SPPR_new` (`TE`, `DET_TE_vals`,
  `det_collapse_mode`, `det_open_mode`, `det_theta`, `det_external_sppr`, `fix_EE_0_cases`). Two are
  worth flagging here: `TE` is how a single Monte-Carlo draw is screened — low-TE draws are precisely
  the ones that push `b` over 1, the same rejection signal `monte_carlo_SPPR` detects after the fact
  from negative SPPR — and `det_external_sppr` enters only `c`, never `B`, so it cannot change `b` or
  the convergence verdict at all.

#### The report

A `dict` of three graded sections plus, unless `short=True`, the ungraded `footprint`, the echoed
`config` and a `warnings` list. Each graded section carries its own `'OK'` / `'WARN'` / `'FAIL'`,
and the top-level `status` is the worst of the three.

**`status`** — `'OK'` (nothing tripped a threshold), `'WARN'` (usable with stated caveats) or
`'FAIL'` (do not interpret the numbers). It is the worst of `model_input`, `divergence` and
`balance`; `footprint` is deliberately excluded, because a large `PPR/NPP` is a finding about the
ecosystem rather than a defect in the calculation.

**`model_input`** — the Ecopath input data. Independent of `TE_option`: these fields do not change
with the configuration.
- **`is_model_balanced`** — `is_model_balanced()`'s strict boolean (an `np.isclose` test on both
  mass-balance identities of §2). Reported for completeness; the grading uses the residuals below,
  because the boolean is all-or-nothing and roughly a fifth of real models fail it by wildly
  differing amounts.
- **`p_max_rel_residual`**, **`q_max_rel_residual`** — the worst group's relative deviation between
  the recomputed and stored production and consumption,
  `max |recomputed − stored| / max(|stored|, 1e-12)`. Scale-free, so one threshold works across
  models of any total throughput. Every SPPR method reads `p`, `q`, `M0` and `predation` as mutually
  consistent; when they are not, the nullspace is solving a food web that does not close.
- **`dc_rows_sum_to_1`**, **`dc_max_deviation`** — whether every consumer's diet row (including
  `diet_import`) sums to 1 within `dc_row_tol`, and the largest observed deviation. Diet fractions
  must partition a consumer's intake or `A = DC/TE` misweights every path through it. Note
  `ModelData.validate_DC` already **raises at load time** at `tol=1e-3`, so anything that loaded
  passes at that tolerance; the default `dc_row_tol=1e-6` is what makes this field informative.
- **`n_negative_catch`**, **`n_zero_catch`**, **`total_catch`**, **`has_catch`** — `PPR = C·SPPR`,
  so a negative catch yields a negative footprint from a perfectly healthy SPPR, and an all-zero
  catch makes every footprint number identically 0 (vacuous rather than wrong). Zero-catch groups
  are normal and are counted, not penalised. Neither condition invalidates the convergence
  diagnostics, so both grade `WARN`.
- **`has_ee_issues`**, **`n_ee0`**, **`ee0_groups`**, **`n_ee_marginal`**,
  **`ee_marginal_groups`**, **`n_ee_gt_1`** — `EE = 0` means a group's production is entirely
  non-predatory death; under `TE_option='TE'` its whole TE row becomes 0, which severs it from the
  nullspace and leaks the PP it consumed (`fix_EE_0_cases` addresses this for single-detritus
  models). *Marginal* EE is `0 < EE < ee_marginal` (default `1e-3`), where `SPPR ~ 1/te` is
  near-singular. `n_ee_gt_1` counts the classic Ecopath over-consumption flag. These are properties
  of the input data, hence reported here whatever `TE_option` is in use, while the
  TE-matrix-dependent consequence is `near_singular_te` below.

**`divergence`** — the solve itself.
- **`b`**, **`b_converges`** — the detritus recycling gain `ρ(diag(θ)·B)` and whether it is `< 1`.
  Measured before the solve-vs-pool decision. Under `TE_option='TE'` there is no recycling matrix,
  so `b` is `0.0` and a note is added to `warnings`.
- **`rho_living`**, **`living_converges`** — `ρ(A_LL)` and whether it is `< 1`, as derived above.
- **`sppr_det`** — `{detritus_seq: sppr_det}`, one entry per detritus pool: the primary production
  charged per unit of that pool. This is the per-source answer to "what is detritus worth".
- **`max_sppr_det`** — the largest of those, graded against `sppr_det_warn` (default `10`).
- **`max_sppr_group`**, **`max_tl_group`** — two landmarks of the solved SPPR vector, reported
  and never graded: `{'seq', 'tl', 'sppr', 'inv_te'}` for the group carrying the largest total
  SPPR, and the same record for the group at the top of the trophic ordering. `inv_te` is `1/te`
  taken from the TE matrix actually in use, so it follows `TE_option` (`p/q` under `'GE'`,
  `(p/q)(1−M0/p)` under `'TE'`, the draw itself when an explicit `TE` was passed) rather than any
  one definition of transfer efficiency. It is the per-step amplification `A = DC/TE` applies to
  that group's edges, and it usually explains a large `sppr`: a group at `1/te ≈ 1900` is
  expensive because each unit of its production needs 1900 units of prey production, not because
  the detritus loop ran away. `None` where `te = 0`, i.e. where the group is severed from the
  nullspace entirely. Trophic levels are taken
  from `get_TL(break_cycles=True, DET_as_PP=True)` — the convention `SPPR_1986` and
  `get_PPR2NPP_ratio` already use — because the cached `TL` attribute is degenerate on real
  models. Since SPPR grows with trophic depth (`SPPR = TE^{1−TL}` in the chain limit, §3), the
  two records naming the same group is the expected picture; when they disagree, something
  other than trophic depth is dominating the solution — a near-singular `te`, or a detritus
  column close to runaway. Both are `None` when the solve produced no SPPR.
- **`n_negative_sources`** — how many basal-source columns contain a negative value; `≥ 1` is
  `FAIL`, since a negative SPPR is not a physical quantity.
- **`expect_negatives`** — the prediction `b ≥ 1`, i.e. whether negatives *should* be there.
- **`near_singular_te`** — group seqs whose transfer efficiency in the TE matrix **actually in use**
  is within `1e-3` of zero, so their `SPPR ~ 1/te` blows up. Depends on `TE_option` and on any
  explicit `TE` draw, unlike the EE fields above.
- **`solve_error`** — `None` normally; otherwise the exception type and message from the failed
  solve, with the remaining fields recovered by the diagnostic re-solve.

**`balance`** — `is_sppr_balanced` on the result (§2, *The balance identity*).
- **`inflow`**, **`outflow`** — primary production entering the system, and the SPPR-weighted
  export leaving it (`catch + growth + net_migration`).
- **`rel_gap`** — `|outflow − inflow| / |inflow|`, graded against `balance_warn` / `balance_fail`
  (defaults `1%` / `5%`). An SPPR that converged but does not close the PP budget is arithmetically
  fine and physically wrong.
- **`is_balanced`** — `is_sppr_balanced`'s strict boolean, reported alongside the gap it is derived
  from.

**`footprint`** — §5's quantities, reported and never graded. `get_NPP(only_inner=False)` raises
`not implemented yet`, so NPP has a single value and the second view is `only_pp` rather than an
outer/inner pair.
- **`ppr_all`** — total PPR over every basal source (`only_inner=False`).
- **`ppr_inner`** — with the `Import` columns dropped (`only_inner=True`): within-system PPR only.
- **`ppr_pp_only`** — with `Import` **and** detritus dropped (`only_pp=True`): the share of the
  footprint traceable to primary producers directly. The gap between `ppr_all` and `ppr_pp_only` is
  how much of the footprint is routed through detritus.
- **`npp`** — `get_NPP(only_inner=True)`, the summed production of the `PP` groups.
- **`ppr2npp`**, **`ppr2npp_pp_only`** — `ppr_inner / npp` and `ppr_pp_only / npp`, the %PPR
  appropriated by the catch under each convention. A ratio above 1 means the catch requires more
  primary production than the system produces, which is a strong hint that something upstream is
  wrong even when every convergence test passed.

**`config`** — what was actually evaluated: `TE_option`, `det_open_mode`, `det_theta`,
`det_external_sppr`, `det_collapse_mode`, `explicit_TE` (whether a `TE` matrix was supplied),
`method` (`'single_detritus'`, `'multi_detritus'` or `'pooled_detritus_scaling'`, from
`detritus_resolution_info`), `would_pool`, and `model` (name and year, `None` for toy or
`from_dict` instances). Present because the verdict is only meaningful together with the
configuration that produced it.

**`warnings`** — one human-readable string per tripped threshold, naming the quantity, its value and
the threshold it crossed. Where a threshold is a bare ratio, the warning also quotes the magnitude
that ratio is inflating: the worst `sppr_det` on the `b` warning, and the largest group SPPR with its
seq and TL on the `rho_living` warning. A near-divergence number means little on its own — `ρ(A_LL) =
0.9995` reads very differently once you see it comes with a TL-2.16 group charged 10⁵ units of primary
production. The EE warnings likewise name the offending groups rather than counting them.

---

## 5. From SPPR to ecosystem footprint

### `get_PPR(sppr, only_inner=False, only_pp=False)`
```python
get_PPR(sppr, only_inner=False, only_pp=False) -> pd.DataFrame
```
Converts a **per-group SPPR** into the **total primary production required by the catch**. Each
group's SPPR is weighted by how much of it we actually harvest and summed — an inner product
with the catch vector, per basal source *s*:

$$ PPR_s = \mathbf{C}\cdot\mathbf{SPPR}_s = \sum_i C_i\cdot(\mathrm{SPPR}_s)_i, \qquad PPR = \sum_{s\in\text{sources}} PPR_s. $$

The input SPPR is relabelled to seq, reindexed onto the catch, and infinities zeroed.

- **`only_inner`** — drop the Import columns, so only production generated *inside* the system is
  counted (exclude subsidies imported across the boundary).
- **`only_pp`** — drop **both** Import and Detritus columns, counting only genuine within-system
  **primary** production. (`only_pp` is stronger than `only_inner` and overrides it.) Only affects
  a per-source DataFrame input; a Series input is already aggregated over sources.
- **Returns** — **always a 1-row DataFrame**, so the output type never depends on the input. For a
  per-source DataFrame input the columns are the basal sources (a PPR per source); for an
  already-aggregated Series input there is a single `'PPR'` column. Either way, get the scalar total
  with `.sum(axis=1).sum()`. (This unified return type is why every caller downstream — including
  `get_PPR2NPP_ratio` and the notebooks — uses `.sum(axis=1)`.)
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

**What the code computes.** The method is a single line:

```python
get_PPR(sppr, only_inner=True, only_pp=only_pp).sum(axis=1).sum() / get_NPP(only_inner=True)
```

Read piece by piece:
- **Numerator** — `get_PPR(sppr, only_inner=True, ...)` weights each group's SPPR by its catch and
  sums over groups (`C·SPPR`), with `only_inner=True` dropping the Import columns so only
  within-system basal sources remain; `.sum(axis=1).sum()` then collapses the 1-row per-source
  DataFrame to a **single scalar** — the total within-system PPR, $\sum_{s\in\text{inner}} (\mathbf{C}\cdot\mathbf{SPPR})_s = \mathbf{C}\cdot\mathbf{SPPR}_{\text{inner}}$.
- **Denominator** — `get_NPP(only_inner=True)` = $\sum_{i\in PP} p_i$, the summed production of the
  primary-producer groups (the ecosystem's net primary production).

$$ \frac{PPR}{NPP} = \frac{\mathbf{C}\cdot\mathbf{SPPR}_{\text{inner}}}{\sum_{i\in PP} p_i}. $$

**How it relates to the inflow–outflow balance (§4).** Under mass balance, all basal production
entering the system leaves it, both sides measured in SPPR units. The inner form of that identity is:

$$NPP = (\mathbf{N_m}+\mathbf{C}+\mathbf{BA})\cdot\mathbf{SPPR}_{\text{inner}}$$

the net primary production equals the total *inner outflow* (catch + biomass accumulation + net migration, each
weighted by the basal production it required). So the denominator $\sum_{i\in PP} p_i$ **is** that
total outflow, and the numerator $\mathbf{C}\cdot\mathbf{SPPR}_{\text{inner}}$ is precisely the
**catch's slice** of it. The ratio is therefore the share of the primary-production budget the
harvest claims:

$$ \frac{PPR}{NPP} = \frac{\mathbf{C}\cdot\mathbf{SPPR}_{\text{inner}}}{(\mathbf{N_m}+\mathbf{C}+\mathbf{BA})\cdot\mathbf{SPPR}_{\text{inner}}}. $$

With `Nm = 0` (the usual case), the ratio can exceed 1 when biomass is being depleted
(`BA·SPPR < 0`, so the outflow denominator shrinks below the catch term), and falls below 1 when
large outflows (mostly detrital, or biomass accumulation) absorb part of the primary production.

- **`only_pp`** — if True, count only within-system primary production in the numerator (drop
  Import and Detritus columns) for a like-for-like PP/PP comparison.
- **Ecological meaning:** a low ratio means the fishery is a small draw on the ecosystem's
  productive base; a ratio approaching (or exceeding) 1 signals that the harvest is
  appropriating an unsustainable share of primary production.

---

## 6. Quick chooser

### By goal

| Want… | Use | Why |
|-------|-----|-----|
| A one-line classic estimate for a whole fishery | `SPPR_1986` | one catch-averaged TL, fixed 10% TE — the back-of-envelope number |
| Per-group classic estimate | `SPPR_1995` / `SPPR_1995_TL_fix` | each group's own trophic level; `_TL_fix` interpolates between integer TLs |
| Explicit food-chain paths (which chains dominate) | `SPPR_EwE` (`return_paths=True`) | returns the enumerated paths, not just the totals |
| Fast cycle-pruned matrix cousin of EwE (≠ EwE when cycles exist) | `SPPR_EwE_Ulanowicz` | nullspace after Ulanowicz weakest-link removal; not the exact all-cycles answer |
| The 2015 input–output method (full cycles) | `SPPR_2015` | Leontief `(I−A)⁻¹`; detritus folded back onto PP |
| **General numeric solver with full detritus control** (the default workhorse) | `SPPR_new` | all cycles counted, tunable recycling/openness knobs, per-source breakdown |
| Exact symbolic solution / explicit imported-diet cost | `SPPR_symbolic` | SymPy solve; `diet_import_option` keeps import cost explicit |
| Uncertainty bands / unbiased expected PPR | `monte_carlo_SPPR` | resamples TE and averages (the Jensen correction, §4) |
| Total PPR of a catch | `get_PPR` | `C·SPPR`; always a 1-row DataFrame — total via `.sum(axis=1).sum()` |
| Share of NPP appropriated (%PPR) | `get_PPR2NPP_ratio` | within-system `PPR / NPP` |
| **Whether any of the above can be trusted** | `diagnose_sppr` | grades both convergence conditions (`b`, `ρ(A_LL)`), the input data and the PP budget in one call (§4) |

### Method comparison at a glance

| Method | Family | Cycles | Per-source breakdown | Detritus knobs | Notes |
|--------|--------|--------|----------------------|----------------|-------|
| `SPPR_1986` | trophic-level | via TL inversion | no (one number) | — | fixed TE=0.1, catch-averaged TL |
| `SPPR_1995` | trophic-level | via TL inversion | no | — | per-group TL, global TE |
| `SPPR_1995_TL_fix` | trophic-level | via TL inversion | no | — | integer-TL interpolation |
| `SPPR_EwE` | flow, path-enum | **simple paths only** (undercounts) | yes | — | returns explicit paths; `use_EE` |
| `SPPR_EwE_Ulanowicz` | flow, nullspace | pruned (weakest-link removed) | yes | — | matrix cousin of EwE; not exact all-cycles |
| `SPPR_2015` | flow, Leontief | **full** | PP columns (DET folded in) | fixed convention | reference 2015 implementation |
| `SPPR_new` | flow, nullspace | **full** | yes (PP / DET / import) | `det_collapse_mode`, `det_open_mode`, `det_theta`, `det_external_sppr` | primary numeric solver |
| `SPPR_symbolic` | flow, symbolic | **full** | yes | same as `SPPR_new` | `diet_import_option='as_DC'/'as_PP'` |
| `monte_carlo_SPPR` | wrapper | inherits solver | inherits solver | forwarded | resamples TE; returns mean + samples |
| `diagnose_sppr` | diagnostic | tests both `ρ(B)` and `ρ(A_LL)` | reads `sppr_det` per pool | forwarded | grades a *configuration*, not a model; never raises on a sick model |

**Rule of thumb:** start with `SPPR_new` (default `TE_option='GE'`) for a modern, cycle-correct,
per-source estimate; switch to `SPPR_1995` when you only need the classic TL-based number, to
`SPPR_EwE` when you want to *see* the dominant food chains, and wrap any of them in
`monte_carlo_SPPR` when you need uncertainty. Run `diagnose_sppr` on any model you have not
vetted before reading its numbers: convergence is not visible in the output, and a model can
converge cleanly while charging ten times the primary production it should. See §2 for the `TE_option` choice (the single most
important assumption) and §4 for how the detritus knobs reshape the recycling solve.
