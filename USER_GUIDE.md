# User Guide: `ModelData` and `PPRCalculator`

This guide walks a new user through the two core classes of FishEstimationAI:

- **`ModelData`** (in `ModelData.py`) — loads one Ecopath ecosystem model from disk and exposes its
  raw tables (group parameters, diet composition, detritus fate, name lookups).
- **`PPRCalculator`** (in `PPRCalculator.py`) — takes a model, completes its mass balance, and computes
  **SPPR** (Specific Primary Production Required) by a family of published and in-house methods.

It explains every public method, its parameters, its return type, **and the ecological meaning** of the
options, with human-readable equations where they clarify the math.

---

## 1. Background: what is PPR / SPPR?

An Ecopath model describes a marine food web as a set of **functional groups** (species or aggregates of
species) connected by who-eats-whom flows. Each group is in **mass balance**: what it produces equals what
is removed from it (eaten, caught, dies, migrates).

- **PPR** (Primary Production Required) = how much phytoplankton primary production must ultimately be
  fixed to support a given fish catch, after accounting for the energy lost at every trophic step.
- **SPPR** (Specific PPR) = PPR **per unit of production** of a group. It answers: "for one tonne of this
  group's production, how many tonnes of primary production were required upstream?"
- **TE** (Transfer Efficiency) = the fraction of energy that survives one trophic transfer (predator
  production ÷ prey production). Roughly 10% in classic estimates.

The headline relationship for a single linear food chain is:

```
SPPR(group) = TE ^ (1 - TL)
```

where `TL` is the trophic level. With `TE = 0.1` and `TL = 3`, `SPPR = 0.1^(-2) = 100`: 100 units of
primary production support 1 unit of this group's production. For a real (branching, recycling) web the
chain formula is replaced by a matrix solve, which is what most methods here do.

Total PPR for the catch is then:

```
PPR = sum over groups of ( catch(group) * SPPR(group) )
```

---

## 2. Key ecological quantities and group types

Every group carries a **trophic_info** label that controls how it is treated:

| `trophic_info` | Meaning | Role in the math |
|----------------|---------|------------------|
| `Regular`      | A consumer (fish, invertebrate, etc.) | Has predators and prey; full mass balance applies. |
| `PP`           | Primary producer (phytoplankton, algae) | A **basal source**: SPPR = 1 by definition (it *is* primary production). |
| `DET`          | Detritus (dead organic matter) | A recycling pool; its SPPR depends on what dies/egests into it. |
| `Import`       | Imported diet (food entering from outside the modeled area) | An external production source; a synthetic "diet_import" group is added automatically. |

Per-group flow quantities used throughout (all in the same currency, e.g. t/km²/year):

| Symbol | Attribute | Meaning |
|--------|-----------|---------|
| `p`    | `model.p` | Production. |
| `q`    | `model.q` | Consumption (food eaten). |
| `catch`| `model.catch` | Fishery landings. |
| `predation` | `model.predation` | Production lost to predators. |
| `M0`   | `model.M0` | Other (non-predation) mortality. |
| `egestion` | `model.egestion` | Unassimilated food (faeces) → detritus. |
| `respiration` | `model.respiration` | Energy lost to respiration. |
| `growth` | `model.growth` | Biomass accumulation. |
| `net_migration` | `model.net_migration` | Emigration − immigration. |
| `EE`   | `model.EE` | Ecotrophic efficiency = fraction of production used in the system, `1 - M0/p`. |
| `GE`   | `model.GE` | Gross efficiency = `p/q`. |
| `GS`   | `model.GS` | Unassimilated fraction = `egestion/q`. |
| `TL`   | `model.TL` | Trophic level. |

Two mass-balance identities link them (and the calculator enforces them):

```
production:   p = catch + predation + growth + net_migration + M0
consumption:  q = p + egestion + respiration
```

---

## 3. `ModelData` — loading a model

`ModelData` is a thin container: it reads one model from disk and exposes its tables. You usually do not
call it directly — `PPRCalculator` wraps it — but it is useful for inspecting raw data.

### 3.1 Constructor

```python
from ModelData import ModelData

md = ModelData(model_input)
```

**`ModelData(model_input: int | str)`**

| Parameter | Type | Meaning |
|-----------|------|---------|
| `model_input` | `int` | **Legacy API**: a model number; data is pulled from the bundled `real_models/SpeciesGroups.json` / diet data. |
| `model_input` | `str` | **New API**: a path to a per-model JSON file, e.g. `"real_models/EwE_jsons/227_Iceland_(1950).json"`. The filename must follow `{number}_{name}_({year}).json`. |

> **Where the files live:** real EwE models are JSON files under `real_models/EwE_jsons/`; small
> hand-built toy models (used in the test notebooks) are under `real_models/ToyModels/`.

**Returns:** `None` (populates the instance in place).
**Raises:** `TypeError` if `model_input` is neither `int` nor `str`; `ValueError` if no group matches a
given model number or the filename can't be parsed.

Both paths automatically inject a synthetic **`diet_import`** group (an extra `Import` row/column) so that
food imported from outside the modeled area can be handled uniformly downstream.

### 3.2 Attributes available after construction

| Attribute | Type | Meaning |
|-----------|------|---------|
| `md.groups_data` | `pd.DataFrame` | Per-group ecological parameters, indexed by `group_seq` (descending). |
| `md.DC` | `pd.DataFrame` | **Diet composition** matrix. `DC[i, j]` = fraction of predator `i`'s diet made up of prey `j`. Rows sum to ~1 for consumers. |
| `md.det_fate` | `pd.DataFrame` | **Detritus fate** matrix: for each group, the fraction of its dead/egested material routed to each detritus pool. |
| `md.seq2name` | `dict` | `group_seq -> group_name`. |
| `md.name2seq` | `dict` | `group_name -> group_seq`. |
| `md.model_number`, `md.model_name`, `md.model_year`, `md.model_country`, `md.lme` | metadata | Model identity. |

### 3.3 Useful methods

- **`md.get_groups_df(species_groups)`** *(staticmethod)* — convert parsed group records into the
  `(groups_data, diet_import)` pair. Rarely called by users.
- **`md.get_DC(data_json)`** *(staticmethod)* — `tuple[pd.DataFrame, pd.DataFrame]`: `(DC, det_fate)`.
- **`md.get_seq2name(data_json)`** *(staticmethod)* — the `seq -> name` mapping (with a `diet_import`
  entry appended).

> **Ecological note on `DC` vs `det_fate`:** `DC` tells you what each group *eats*; `det_fate` tells you
> where each group's *waste and dead biomass* goes. Both are needed because detritus is recycled back into
> the food web, and the recycling loop is what makes the PPR of detritus-feeders non-trivial.

---

## 4. `PPRCalculator` — building it

`PPRCalculator` loads a model, completes any missing flows so the model is mass-balanced, and then offers
the SPPR methods. There are three constructors.

### 4.1 Primary constructor

```python
from PPRCalculator import PPRCalculator

model = PPRCalculator(model_number)        # e.g. PPRCalculator(227)
model = PPRCalculator("real_models/227_Iceland_(1950).json")
```

**`PPRCalculator(model_number, underdetermined=False, zero_catch=True, zero_biomass_accum=True,
default_gs=True, weight_flow=1.0, weight_guess=1.0)`**

| Parameter | Type | Default | Ecological / numerical meaning |
|-----------|------|---------|-------------------------------|
| `model_number` | `int \| str` | — | Model number or JSON filepath (passed straight to `ModelData`). |
| `underdetermined` | `bool` | `False` | If `True`, missing flows are filled by a **Linear Inverse Model** (`apply_lim`, an optimizer) instead of only the deterministic Ecopath defaults. Use when the published model leaves several cells blank. |
| `zero_catch` | `bool` | `True` | Treat missing `catch` as 0 (most groups are not fished). |
| `zero_biomass_accum` | `bool` | `True` | Treat missing biomass accumulation as 0 (assume steady state). |
| `default_gs` | `bool` | `True` | Assign the textbook **GS = 0.2** (20% of food unassimilated) to regular groups lacking a value. |
| `weight_flow` | `float` | `1.0` | LIM penalty weight favoring the *smallest total flows* (parsimony). |
| `weight_guess` | `float` | `1.0` | LIM penalty weight favoring *staying near biologically sensible guesses*. |

**Returns:** a ready-to-use `PPRCalculator`.
**Raises:** `Exception` if the model has no detritus (`DET`) group.

### 4.2 Alternative constructors

- **`PPRCalculator.from_modeldata(modeldata, ...)`** → `PPRCalculator`. Build from an already-loaded
  `ModelData`. Same keyword options as above. This is the core path `__init__` delegates to:

  ```python
  from ModelData import ModelData
  md = ModelData("real_models/EwE_jsons/10013_Humboldt_Current_(1995-2004).json")
  model = PPRCalculator.from_modeldata(md)
  ```

- **`PPRCalculator.from_dict(data_dict, ...)`** → `PPRCalculator`. Rebuild from a dict of pre-existing
  attributes (e.g. a hand-made toy model or a deserialized state). Re-runs the full
  defaults → LIM → fill → sort pipeline. `data_dict` must at least contain a valid `_groups_df`. This is
  how the test notebooks build minimal toy webs:

  ```python
  # A 3-group toy: PP -> A, with detritus DET.
  model = PPRCalculator.from_dict({
      '_groups_df': groups_df,   # columns: group_name, trophic_info, q, p, catch, predation, ...
      '_DC': DC,                 # diet composition matrix
      'seq2name': seq2name, 'name2seq': name2seq,
  })
  ```

> **What "completing the model" does:** Ecopath models are mass-balanced by construction, but published
> tables often omit cells that are derivable. `apply_ecopath_defaults` fills one missing cell per balance
> equation (e.g. derives `M0 = p·(1−EE)`, `egestion = q·GS`, `predation = column sum of the flow matrix`).
> If too many cells are missing for that to work, set `underdetermined=True` to invoke the optimizer.

---

## 5. `PPRCalculator` — inspecting the model

These getters expose the completed model. All are cheap and return copies / fresh objects.

| Method | Returns | Meaning |
|--------|---------|---------|
| `get_model()` | `ModelData \| None` | The backing `ModelData` (or `None` for `from_dict` toy models). |
| `get_groups_df()` | `pd.DataFrame` | A sorted copy of the per-group parameter table. |
| `get_DC(DET_as_PP=True, normalize=False)` | `pd.DataFrame` | Diet composition matrix (see below). |
| `get_Z(DET_as_PP=False)` | `pd.DataFrame` | Flow matrix `Z = DC · q`: absolute prey→predator flows. |
| `get_DET_seq()` | `list` | Seq IDs of detritus groups. |
| `get_PP_seq()` | `list` | Seq IDs of primary producers. |
| `get_Regular_seq()` | `list` | Seq IDs of regular consumers. |
| `get_Import_seq()` | `list` | Seq IDs of imported-diet groups. |
| `get_TE(TE_option, DET_values=1, as_matrix=True, global_TE='mean')` | `pd.DataFrame \| pd.Series` | Transfer-efficiency vector/matrix (see §5.1). |
| `get_TL(break_cycles, DET_as_PP, TE_option='With Egestion')` | `pd.Series` | Trophic level per group (see §5.2). |
| `get_NPP(only_inner=True)` | `float` | Net primary production = total production of PP groups. |
| `get_PPR(sppr, only_inner=False)` | `pd.DataFrame \| pd.Series` | Total PPR = `catch · SPPR` (see §5.3). |
| `get_PPR2NPP_ratio(sppr)` | `float` | Fraction of NPP appropriated by the catch = PPR / NPP. |

### 5.1 `get_TE` and the `TE_option` choices

`TE_option` selects **how transfer efficiency is defined** — this is the single most important ecological
modeling choice in the package, because it decides what counts as "lost" energy at each step:

| `TE_option` | Formula | Ecological meaning |
|-------------|---------|--------------------|
| `'GE'` | `TE = p / q` (gross efficiency) | Everything not turned into production (respiration **and** egestion) is treated as a transfer loss. The detritus feeds back into the web. |
| `'TE'` | `TE = (p/q) · (1 − M0/p) = GE · EE` | Counts only the production that is actually *used* by the system (predation/export). Energy lost to non-predation mortality is excluded. This reproduces the classic 2015-style accounting. |
| `'With Egestion'` | `TE = (p/q) · (q / (q − egestion))` | Efficiency on *assimilated* intake: egestion is sent to detritus and recycled rather than counted as a dead loss, so assimilation is the denominator. |
| `'global'` | a single scalar broadcast to all groups | Uses one fixed efficiency for the whole web (classic "≈10%" assumption). Set by `global_TE`. |

`global_TE` (only used when `TE_option='global'`): either a literal `float` (e.g. `0.1`) or `'mean'`,
in which case it is the catch-weighted mean of the per-group `'TE'` efficiency (biomass-weighted if total
catch is 0).

`DET_values` sets the TE assigned to detritus rows (default `1`: detritus is fully available as a basal
source). `as_matrix=True` returns an n×n matrix (the vector broadcast across columns); `False` returns the
length-n vector.

### 5.2 `get_TL` — trophic levels

Trophic level is solved from the diet matrix as

```
TL = (I − DC)^(-1) · 1
```

i.e. each group's TL is one plus the diet-weighted mean TL of its prey; basal sources sit at TL 1.

- `break_cycles` (`bool`): if `True`, remove cycles from the flow matrix (Ulanowicz algorithm) before
  inverting. **Why:** food webs contain loops (A eats B eats detritus eats A…); these make the naive
  inversion unstable or push TLs to infinity. Breaking cycles yields finite, well-defined levels.
- `DET_as_PP` (`bool`): whether detritus is treated as a basal source (TL 1) rather than a consumer.
- `TE_option` (`str`): how the detritus rows of `DC` are redefined before inversion.

### 5.3 `get_PPR` — turning SPPR into total PPR

```python
sppr = model.SPPR_new()[0]          # per-group SPPR
ppr  = model.get_PPR(sppr)          # total PPR weighted by catch
```

`only_inner=True` drops the `Import` columns so only **within-system** primary production is counted (food
that came in pre-made from outside is excluded). `get_PPR2NPP_ratio` then expresses that as a fraction of
the system's own NPP — a number > 1 means the catch demands more primary production than the system makes
internally (i.e. it relies on imports).

---

## 6. `PPRCalculator` — SPPR methods

All SPPR methods return a `pd.DataFrame` (or `Series`) of SPPR values. The matrix-based methods return a
**tuple** whose first element is the SPPR table; unpack accordingly.

### 6.1 Quick reference

| Method | Returns | One-line description |
|--------|---------|----------------------|
| `SPPR_1986()` | `DataFrame` | Single catch-weighted mean TL, `TE=0.1`, `SPPR = TE^(1−TL)`. |
| `SPPR_1995(global_TE=0.1)` | `DataFrame` | Per-group `SPPR = TE^(1−TL)` with continuous TL. |
| `SPPR_1995_TL_fix(global_TE=0.1)` | `DataFrame` | 1995 with linear interpolation between integer TLs. |
| `SPPR_2015(only_pp_det=True)` | `(SPPR, A, L)` | Leontief matrix-inversion formulation. |
| `SPPR_EwE(TE_option, use_EE=True, return_paths=True, silent=True)` | `(SPPR, A, paths)` | Explicit path enumeration over the network. |
| `SPPR_EwE_Ulanowicz(TE_option, global_TE='mean', use_EE=True)` | `(SPPR, A, L)` | Nullspace reformulation of the EwE path sum. |
| `SPPR_new(...)` | `(SPPR, A, L)` | Primary numeric solver with full detritus recycling. |
| `SPPR_symbolic(...)` | `(sppr_symbolic, sppr_mat, equations, variables)` | Symbolic solver (keeps imported diet explicit). |
| `monte_carlo_SPPR(...)` | `(mean, samples, reject_frac, eqs, vars)` | Uncertainty propagation over TE. |
| `monte_carlo_SPPR_2(...)` | `(mean, samples, reject_frac)` | Same, `kind='new'` only. |

### 6.2 The classic chain methods

**`SPPR_1986()` → `pd.DataFrame`** (one `'sppr'` column).
Pauly & Christensen (1986). One catch-weighted mean trophic level for the whole catch, fixed `TE = 0.1`:

```
TL_mean = catch-weighted mean of TL
SPPR    = 0.1 ^ (1 − TL_mean)        (same value for every group)
```

Returns all zeros if there is no catch. This is the coarsest estimate.

**`SPPR_1995(global_TE=0.1)` → `pd.DataFrame`.**
Christensen & Pauly (1995). Uses each group's *own* continuous trophic level:

```
SPPR(group) = TE ^ (1 − TL(group))
```

`global_TE` is the single global efficiency (`float`, or `'mean'` for the weighted-mean value).

**`SPPR_1995_TL_fix(global_TE=0.1)` → `pd.DataFrame`.**
Same idea but interpolates between the two integer trophic levels bracketing a fractional TL:

```
SPPR = (1 − frac)·(1/TE)^(TLint − 1) + frac·(1/TE)^TLint
```

**Why:** raising `1/TE` to a non-integer power directly overweights omnivores (Jensen's inequality). The
interpolation keeps the estimate between the integer-level chains and avoids that bias.

### 6.3 The network methods

These replace the single-chain formula with the real branching web, so a group's PPR is summed over **all**
the pathways by which energy reaches it.

**`SPPR_EwE(TE_option, use_EE=True, return_paths=True, silent=True)` → `tuple[DataFrame, DataFrame, dict]`.**
Path-enumeration in the style of Ecopath with Ecosim. It builds the per-edge weight matrix `A = DC / TE`
and, for each (group, basal source) pair, sums the product of edge weights over every simple path:

```
SPPR(group, source) = Σ over paths  Π over edges  A[edge]
```

- `TE_option`: one of `'GE'`, `'TE'`, `'With Egestion'`, `'global'` (see §5.1).
- `use_EE`: scale each row by ecotrophic efficiency `EE` (account for the fraction of production not used
  within the system).
- `return_paths`: `True` returns the explicit path lists (slower) in the third element; `False` uses a fast
  vectorized solver and returns an empty dict.
- `silent`: suppress progress bars.

Returns `(SPPR, A, paths_dict)`.

**`SPPR_EwE_Ulanowicz(TE_option, global_TE='mean', use_EE=True)` → `tuple[DataFrame, DataFrame, DataFrame]`.**
A matrix reformulation that gives the *same* answer as the path sum without enumerating paths. It builds
`A = DC/TE` (cycles removed), replaces basal rows with identity rows, and solves for the steady state as the
**nullspace** of `L = A − I` (so `A·x = x`). The nullspace basis is RREF-normalized so each output column is
anchored to one basal source. Returns `(SPPR, A, L)`.
**Raises** `ValueError` if no steady state exists (disconnected web).

**`SPPR_2015(only_pp_det=True)` → `tuple[DataFrame, DataFrame, DataFrame]`.**
A Leontief input–output formulation. Detritus columns are dissolved by reassigning the PP-derived fraction
of each detritus flow back onto the PP groups, leaving only living compartments. With the
production-normalized transaction matrix `A`:

```
L = (I − A)^(-1)
```

The PP columns of `L` give per-group SPPR; a balancing detritus SPPR is added back at the end. Returns
`(SPPR, A, L)`. (Note: `only_pp_det` is forced to `True` internally, matching the published method.)

### 6.4 `SPPR_new` — the primary solver

```python
SPPR, A, L = model.SPPR_new(TE_option='GE')
```

`SPPR_new` is the recommended general-purpose solver. It builds `A = DC/TE`, replaces basal rows with
identity rows, and finds the steady-state SPPR as the nullspace of `L = A − I` (RREF-normalized per basal
source). Then it **resolves detritus recycling** depending on `TE_option`:

- `'TE'`: each detritus column is scaled by its direct PP+Import inflow share (no feedback loop).
- `'GE'` / `'With Egestion'`: detritus genuinely recycles, so the solver builds and solves a coupled linear
  system (§6.6).

**Full signature:**

```python
SPPR_new(TE=None, TE_option='GE', DET_TE_vals=1, collapse_det=None,
         det_collapse_mode='never', det_open_mode='none',
         det_theta=1.0, det_external_sppr=0.0)
```

| Parameter | Type | Default | Meaning |
|-----------|------|---------|---------|
| `TE` | `pd.DataFrame \| None` | `None` | An explicit TE matrix (e.g. a Monte-Carlo sample). If `None`, built from `TE_option`. |
| `TE_option` | `str` | `'GE'` | `'GE'`, `'TE'`, `'With Egestion'`, `'global'` (see §5.1). |
| `DET_TE_vals` | `float` | `1` | TE assigned to detritus rows when building the TE matrix. |
| `collapse_det` | `bool \| None` | `None` | **Legacy** boolean: `False` → `det_collapse_mode='never'`, `True` → `'auto'`, `None` → leave `det_collapse_mode` as given. |
| `det_collapse_mode` | `str` | `'never'` | Detritus stability strategy (see §6.5). |
| `det_open_mode` | `str` | `'none'` | Detritus openness model (see §6.5). |
| `det_theta` | `float \| dict` | `1.0` | Detritus availability/retention fraction (see §6.5). |
| `det_external_sppr` | `float \| dict` | `0.0` | External SPPR for diluted material under `'source_dilution'` (see §6.5). |

**Returns** `(SPPR, A, L)`.
**Raises** `ValueError` (empty nullspace) or `Exception` (bad `TE_option`).

### 6.5 The detritus knobs (shared by `SPPR_new`, `SPPR_symbolic`, the symbolic helpers, and both Monte-Carlo methods)

Detritus is a recycling pool: dead biomass and faeces flow into it, and detritus-feeders eat from it, which
sends energy back up the web. This loop can amplify SPPR without bound if recycling is too strong, so these
knobs control **how the loop is closed and stabilized**.

**`det_collapse_mode`** — what to do when the recycling system is numerically unstable:

| Value | Behaviour | When to use |
|-------|-----------|-------------|
| `'never'` (default) | Always solve the full coupled detritus system directly. May return **negative** SPPR for unstable models, but **never raises**. The Monte-Carlo samplers rely on this so they can detect and reject unstable draws. | Default; keeps results exact and lets you see instability. |
| `'auto'` | Solve directly **unless** the system is unstable — spectral radius `ρ(B) ≥ 1` or the matrix is ill-conditioned — in which case fall back to a single pooled detritus scaling. | When you want robust, always-finite results. |
| `'always'` | Always use the pooled single-scalar detritus scaling. | Quick, very robust approximation. |

**`det_open_mode`** — how "open" the recycling loop is (how much recycled detritus is actually re-used):

| Value | Effect | Ecological meaning |
|-------|--------|--------------------|
| `'none'` (default) | Closed recycling: all detritus is recycled. | Idealized closed system. |
| `'recycling_loss'` | Damp the recycling matrix `B` by `diag(theta)`. | A fraction `1 − θ` of detritus is buried/exported and lost from the loop. |
| `'source_dilution'` | Damp `B` by `diag(theta)` **and** dilute the source term toward an external SPPR. | As above, but the lost fraction is replaced by material of a known external SPPR value. |

**`det_theta`** — the availability/retention fraction θ (how much detritus stays in the loop):
a single `float` applied to all detritus groups, **or** a `dict` keyed by detritus-group seq (`int`) or
name (`str`). `θ = 1.0` (default) reproduces the fully closed system.

**`det_external_sppr`** — the SPPR assigned to the diluted (replacement) material under
`'source_dilution'`. Same scalar-or-`dict` form as `det_theta`; default `0.0`.

The openness transform applied to the recycling system `(I − B)x = c` is, with θ and `ext` aligned to the
detritus groups:

```
none            : B_open = B,             c_open = c
recycling_loss  : B_open = diag(θ)·B,     c_open = c
source_dilution : B_open = diag(θ)·B,     c_open = θ·c + ext·(1 − θ)
```

With `θ = 1` and `ext = 0` all three reduce to the original closed system.

### 6.6 What the coupled detritus solve actually computes (for the curious)

For `'GE'` / `'With Egestion'`, each detritus group `l` gets an unknown scaling `x_l`. Recycled detritus
production depends on the SPPR already attributed elsewhere, giving:

```
x_l     = (production from non-detritus sources) + Σ_j B[l,j] · x_j
(I − B)·x = c
```

where `c[l]` is the inflow from non-detritus sources and `B[l,j]` is detritus `l`'s dependence on detritus
`j`. The solver tests the **spectral radius** `ρ(B)`: if `ρ(B) < 1` the loop is sub-critical and a finite,
non-negative solution exists; `ρ(B) ≥ 1` means runaway recycling, which is when `'auto'`/`'always'` collapse
to the pooled scaling. The diagnostics of which path was taken are stored on
`model.detritus_resolution_info` after each call.

### 6.7 `SPPR_symbolic` — keeping imported diet explicit

```python
sppr_symbolic, sppr_mat, equations, variables = model.SPPR_symbolic(
    TE_option='GE', diet_import_option='as_DC')
```

Solves the same steady state symbolically (with `sympy`), which lets imported diet be tracked exactly.

**Full signature:**

```python
SPPR_symbolic(TE=None, TE_option='GE', diet_import_option='as_DC', DET_TE_vals=1,
              sppr_det_value=None, collapse_det=None,
              det_collapse_mode='never', det_open_mode='none',
              det_theta=1.0, det_external_sppr=0.0)
```

The new parameter is **`diet_import_option`**:

| Value | Meaning |
|-------|---------|
| `'as_DC'` (default) | Imported diet is kept as a **separate production source** with its own `DIET_SPPR_*`, solved from a second linear system. Each imported group then carries the diet-composition-weighted production it actually requires. This is the more faithful treatment, and the only one for which `is_sppr_balanced` can check the diet-import balance. |
| `'as_PP'` | Imported diet is treated like an **extra primary producer** (its own free SPPR symbol fixed to 1). Simpler, but does not propagate the upstream cost of the imported food. |

`sppr_det_value` (`float | None`): if set, every detritus column is scaled by this fixed value instead of
being solved. All the detritus knobs from §6.5 apply here too.

**Returns** `(sppr_symbolic, sppr_mat, equations, variables)`:
- `sppr_symbolic` (`DataFrame`): the symbolic per-group solution.
- `sppr_mat` (`DataFrame`): the numeric basis matrix (groups × basal sources).
- `equations`, `variables` (`list`): the symbolic system — pass these to `is_sppr_balanced(...,
  diet_import_equations=(equations, variables))` to verify diet-import-aware balance.

### 6.8 Monte-Carlo uncertainty

Transfer efficiency is uncertain, so these methods propagate that uncertainty into SPPR.

**`monte_carlo_SPPR(n_samples=1000, TE_error_percent=10, TE_error_cut_percent=20, TE_option='GE',
DET_TE_vals=1, kind='new', diet_import_option='as_DC', silent=True, <det knobs>)`**

Repeatedly resamples the TE matrix, recomputes SPPR, **rejects any sample with a negative SPPR** (an
unstable/non-physical draw), and averages the survivors.

| Parameter | Type | Default | Meaning |
|-----------|------|---------|---------|
| `n_samples` | `int` | `1000` | Number of TE draws. |
| `TE_error_percent` | `float` | `10` | TE standard deviation as a % of its mean (the gamma CV). |
| `TE_error_cut_percent` | `float` | `20` | Clip each TE sample to ±this % of its mean (drop extreme tails). |
| `TE_option` | `str` | `'GE'` | TE mode (see §5.1). |
| `DET_TE_vals` | `float` | `1` | TE of detritus rows. |
| `kind` | `str` | `'new'` | `'new'` (uses `SPPR_new`) or `'symbolic'` (uses `SPPR_symbolic`). |
| `diet_import_option` | `str` | `'as_DC'` | Passed to `SPPR_symbolic` when `kind='symbolic'`. |
| `silent` | `bool` | `True` | Suppress progress output. |
| `det_collapse_mode`, `det_open_mode`, `det_theta`, `det_external_sppr` | | | Forwarded to every SPPR call (see §6.5). |

**Why gamma sampling?** TE must stay strictly positive. The gamma distribution does this while letting you
specify mean and coefficient of variation directly: with `shape = 1/CV²` and `scale = TE_mean/shape`, the
draws have mean `TE_mean` and `std/mean = CV`. Basal rows are pinned to 1 each draw.

**Returns** `(mean_sppr, accepted_samples_array, rejection_fraction, equations, variables)`. The last two
are `None` when `kind='new'`. `rejection_fraction` tells you how often the model was unstable — a high value
is a red flag about the model or the chosen TE error.
**Raises** `Exception` if `kind` is not `'new'` or `'symbolic'`.

**`monte_carlo_SPPR_2(...)`** is the same loop restricted to `kind='new'`; it pre-allocates from the model
shape and returns the 3-tuple `(mean_sppr, accepted_samples_array, rejection_fraction)`.

---

## 7. Validation helpers

| Method | Returns | Meaning |
|--------|---------|---------|
| `is_model_balanced()` | `(bool, production, consumption)` | Checks the two mass-balance identities hold; returns the recomputed flow Series so you can inspect residuals. |
| `is_sppr_balanced(sppr, diet_import_equations=None)` | `(bool, inflow, outflow)` | Checks an SPPR result is globally self-consistent: primary-production **inflow** equals the export-weighted **outflow**. Pass `(equations, variables)` from `SPPR_symbolic` to use the diet-import-aware check. |
| `balance_model(change_production=False)` | `PPRCalculator` | Returns a deep-copied, exactly mass-balanced version of the model. `change_production=False` absorbs residuals into growth/net_migration (keeping `p` fixed); `True` absorbs them into production. |

The global SPPR balance check verifies:

```
inflow  = total PP production entering the system
outflow = Σ (catch + growth + net_migration) · SPPR
is_balanced = isclose(inflow, outflow)
```

A balanced SPPR means the production you pulled out (as catch, growth, migration) is exactly the production
that came in as primary production — a strong sanity check on any method.

---

## 8. Utilities

**`PPRCalculator.rename_results(results, renaming_dict)`** *(classmethod)* →
`list | DataFrame | Series`. Relabels the index (and columns of DataFrames) of one result or a list of
results, typically with `model.seq2name` (seq → name) or `model.name2seq`. Returns the same type it was
given.

```python
named = PPRCalculator.rename_results(sppr, model.seq2name)
```

---

## 9. End-to-end example

```python
from PPRCalculator import PPRCalculator

# 1. Load and complete a model.
model = PPRCalculator(227)

# 2. Inspect it.
print("groups:", model.n_groups)
print("PP groups:", model.get_PP_seq())
print("balanced?", model.is_model_balanced()[0])

# 3. Compute SPPR several ways.
sppr_1995 = model.SPPR_1995(global_TE=0.1)             # classic chain method
sppr_new, A, L = model.SPPR_new(TE_option='GE')        # full network solver

# 4. Detritus openness sensitivity: lose 30% of recycled detritus.
sppr_open, _, _ = model.SPPR_new(
    TE_option='GE', det_open_mode='recycling_loss', det_theta=0.7)

# 5. Turn SPPR into total PPR and the PPR/NPP ratio.
print("total PPR:", model.get_PPR(sppr_new))
print("PPR/NPP:", model.get_PPR2NPP_ratio(sppr_new))

# 6. Uncertainty bounds (robust mode).
mean, samples, reject_frac, _, _ = model.monte_carlo_SPPR(
    n_samples=200, TE_error_percent=10, kind='new',
    det_collapse_mode='auto')
print(f"rejected {reject_frac:.0%} of unstable draws")

# 7. Verify a result.
print("SPPR balanced?", model.is_sppr_balanced(sppr_new)[0])
```

---

## 10. Choosing a method (rules of thumb)

- **Quick, comparable to the literature** → `SPPR_1995` (single global TE, per-group TL).
- **Full food-web accounting with recycling** → `SPPR_new` with `TE_option='GE'` or `'With Egestion'`.
- **Reproduce the 2015 paper** → `SPPR_2015`, or `SPPR_new(TE_option='TE')` (they should agree).
- **Track imported food explicitly** → `SPPR_symbolic(diet_import_option='as_DC')`.
- **Report uncertainty** → `monte_carlo_SPPR(kind='new')`; use `det_collapse_mode='auto'` if many draws are
  rejected.
- **Models with strong detritus recycling that "explode"** → add openness (`det_open_mode='recycling_loss'`,
  `det_theta < 1`) or set `det_collapse_mode='auto'`/`'always'`.
