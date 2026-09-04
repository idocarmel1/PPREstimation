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
| `TL`   | `model.TL` | Trophic level. Models rarely store one, so it is solved on construction via `get_TL(break_cycles=False, DET_as_PP=True)` and any stored value is kept. |

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
| `model_input` | `str` | **New API**: a path to a per-model JSON file, e.g. `"real_models/EwE_jsons/227_227_Iceland_(1950).json"`. The filename stem must follow `{first_number}_{model_number}_{name}_({year})` — **two** leading numeric tokens (the first is a source/grouping id and is discarded; the second is the model id). For single-model files the two numbers are the same (`227_227_...`); for multi-model source files they differ (`13_10013_Humboldt_Current_(1980)`). |

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
model = PPRCalculator("real_models/EwE_jsons/227_227_Iceland_(1950).json")
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
  md = ModelData("real_models/EwE_jsons/13_10013_Humboldt_Current_(1995-2004).json")
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
| `get_NPP(only_inner=True)` | `float` | Net primary production = total production of PP groups (`only_inner=False` raises). |
| `get_PPR(sppr, only_inner=False, only_pp=False)` | `pd.DataFrame` | Total PPR = `catch · SPPR`, **always a 1-row DataFrame** (see §5.3). |
| `get_PPR2NPP_ratio(sppr, only_pp=False)` | `float` | Fraction of NPP appropriated by the catch = PPR / NPP. |

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
sppr = model.SPPR_new()[0]                 # per-group SPPR (DataFrame: groups × basal sources)
ppr_df = model.get_PPR(sppr)               # 1-row DataFrame of PPR per basal source
ppr_total = ppr_df.sum(axis=1).sum()       # scalar total PPR
```

**`get_PPR(sppr, only_inner=False, only_pp=False)` → `pd.DataFrame` (always a 1-row DataFrame).**
The return type no longer depends on the input: a per-source **DataFrame** input yields one column per
basal source (PP / detritus / import); an already-aggregated **Series** input yields a single `'PPR'`
column. Get the scalar total with `.sum(axis=1).sum()` in either case.

- `only_inner=True` drops the `Import` columns so only **within-system** production is counted (food that
  came in pre-made from outside is excluded).
- `only_pp=True` drops **both** Import and Detritus columns, counting only genuine within-system
  **primary** production. `only_pp` is stronger than `only_inner` and overrides it. The `only_*` flags only
  affect a per-source DataFrame input (a Series is already aggregated over sources).

`get_PPR2NPP_ratio(sppr, only_pp=False)` expresses within-system PPR as a fraction of the system's own NPP.
See `SPPR_Methods.md` §5 for the full inflow/outflow balance interpretation of this ratio.

---

## 6. `PPRCalculator` — SPPR methods

All SPPR methods return a `pd.DataFrame` (or `Series`) of SPPR values. The matrix-based methods return a
**tuple** whose first element is the SPPR table; unpack accordingly.

> **Deep reference:** this section is a practical API/usage guide. For the full derivations — the
> nullspace foundation, the detritus recycling fixed point `sppr_det = a/(1−b)` and its multi-pool
> `(I − B)x = c` form, the openness transforms, and per-method equations — see `SPPR_Methods.md`
> (§4 covers the flow-network methods and detritus handling in depth).

### 6.1 Quick reference

| Method | Returns | One-line description |
|--------|---------|----------------------|
| `SPPR_1986()` | `DataFrame` | Single catch-weighted mean TL, `TE=0.1`, `SPPR = TE^(1−TL)`. |
| `SPPR_1995(global_TE=0.1)` | `DataFrame` | Per-group `SPPR = TE^(1−TL)` with continuous TL. |
| `SPPR_1995_TL_fix(global_TE=0.1)` | `DataFrame` | 1995 with linear interpolation between integer TLs. |
| `SPPR_2015()` | `(SPPR, A, L)` | Leontief matrix-inversion formulation. |
| `SPPR_EwE(TE_option, use_EE=True, return_paths=True, silent=True)` | `(SPPR, A, paths)` | Explicit path enumeration over the network. |
| `SPPR_EwE_Ulanowicz(TE_option, global_TE='mean', use_EE=True)` | `(SPPR, A, L)` | Nullspace reformulation of the EwE path sum. |
| `SPPR_new(...)` | `(SPPR, A, L)` | Primary numeric solver with full detritus recycling. |
| `SPPR_symbolic(...)` | `(sppr_symbolic, sppr_mat, equations, variables)` | Symbolic solver (keeps imported diet explicit). |
| `monte_carlo_SPPR(...)` | `(mean, samples, reject_frac, eqs, vars[, diagnostics])` | Uncertainty propagation over TE; solver chosen by `kind`, its own parameters via `method_kwargs`. |
| `diagnose_sppr(...)` | `dict` | *Not an SPPR method* — grades whether a `SPPR_new` result is trustworthy (§7.1). |

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

**`SPPR_2015()` → `tuple[DataFrame, DataFrame, DataFrame]`.** (No parameters.)
A Leontief input–output formulation. Detritus columns are dissolved by reassigning the PP-derived fraction
of each detritus flow back onto the PP groups, leaving only living compartments. With the
production-normalized transaction matrix `A`:

```
L = (I − A)^(-1)
```

The PP columns of `L` give per-group SPPR; a balancing detritus SPPR is added back at the end. Returns
`(SPPR, A, L)`. See `SPPR_Methods.md` §4 (`SPPR_2015`) for the proof that `A = Z/P` is the same matrix as
`A = DC/TE`.

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
SPPR_new(TE=None, TE_option='GE', DET_TE_vals=1,
         det_collapse_mode='never', det_open_mode='none',
         det_theta=1.0, det_external_sppr=0.0, fix_EE_0_cases=True)
```

| Parameter | Type | Default | Meaning |
|-----------|------|---------|---------|
| `TE` | `pd.DataFrame \| None` | `None` | An explicit TE matrix (e.g. a Monte-Carlo sample). If `None`, built from `TE_option`. |
| `TE_option` | `str` | `'GE'` | `'GE'`, `'TE'`, `'With Egestion'`, `'global'` (see §5.1). |
| `DET_TE_vals` | `float` | `1` | TE assigned to detritus rows when building the TE matrix. |
| `det_collapse_mode` | `str` | `'never'` | Detritus stability strategy (see §6.5). |
| `det_open_mode` | `str` | `'none'` | Detritus openness model (see §6.5). |
| `det_theta` | `float \| dict` | `1.0` | Detritus availability/retention fraction (see §6.5). |
| `det_external_sppr` | `float \| dict` | `0.0` | External SPPR for diluted material under `'source_dilution'` (see §6.5). |
| `fix_EE_0_cases` | `bool` | `True` | Mass-balance fix for the `TE_option='TE'` case: re-credits the PP consumed by `EE=0` dead-end groups to detritus (single-detritus models only; no-op otherwise). See §6.5 and `SPPR_Methods.md` §4. |

**Returns** `(SPPR, A, L)`.
**Raises** `ValueError` (empty nullspace) or `Exception` (bad `TE_option`).

### 6.5 The detritus knobs (shared by `SPPR_new`, `SPPR_symbolic`, and `monte_carlo_SPPR`)

Detritus is a recycling pool: dead biomass and faeces flow into it, and detritus-feeders eat from it, which
sends energy back up the web. This loop can amplify SPPR without bound if recycling is too strong, so these
knobs control **how the loop is closed and stabilized**. For the full derivation of the recycling solve
(the single-pool fixed point `sppr_det = a/(1−b)` and the multi-pool `(I − B)x = c` system) and of each
openness transform, see `SPPR_Methods.md` §4 (`SPPR_new`).

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

**`fix_EE_0_cases`** (default `True`) — a mass-balance fix specific to `TE_option='TE'`. Groups with
`EE=0` (all production dies non-predatorily) get `TE=0` and are severed from the nullspace, which leaks the
PP they consumed and breaks the global PP balance. When `True`, that consumed PP is re-credited to
detritus. It is only active for **single-detritus** models under `'TE'` (a no-op otherwise) and emits a
`RuntimeWarning` when `EE=0` groups are present. See `SPPR_Methods.md` §4 for details.

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
              sppr_det_value=None, det_collapse_mode='never', det_open_mode='none',
              det_theta=1.0, det_external_sppr=0.0, fix_EE_0_cases=True)
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
DET_TE_vals=1, kind='new', method_kwargs=None, exclude_diverged=False,
return_diagnostics=False, diet_import_option='as_DC', silent=True, <det knobs>)`**

Repeatedly resamples the TE matrix, recomputes SPPR, **rejects any sample with a negative SPPR** (an
unstable/non-physical draw), and averages the survivors. Optionally also rejects draws whose
divergence grade is `FAIL` (`exclude_diverged=True`) — worth doing, because a diverged draw can
still return an all-positive SPPR and slip past the negative test.

| Parameter | Type | Default | Meaning |
|-----------|------|---------|---------|
| `n_samples` | `int` | `1000` | Number of TE draws. |
| `TE_error_percent` | `float` | `10` | TE standard deviation as a % of its mean (the gamma CV). |
| `TE_error_cut_percent` | `float` | `20` | Clip each TE sample to ±this % of its mean (drop extreme tails). |
| `TE_option` | `str` | `'GE'` | TE mode (see §5.1). |
| `DET_TE_vals` | `float` | `1` | TE of detritus rows. |
| `kind` | `str` | `'new'` | `'new'` (uses `SPPR_new`) or `'symbolic'` (uses `SPPR_symbolic`). |
| `method_kwargs` | `dict \| None` | `None` | Extra parameters forwarded to the selected solver on every draw, e.g. `{'fix_EE_0_cases': True}`. A key the solver does not accept raises `TypeError`; a key this wrapper already controls (`TE`, `TE_option`, `DET_TE_vals`, `diet_import_option`, `det_*`) raises `ValueError`. Both before any sampling. |
| `exclude_diverged` | `bool` | `False` | Also reject draws whose `diagnose_sppr` divergence status is `FAIL` (b > 1). Requires `kind='new'`; costs no extra solve. |
| `return_diagnostics` | `bool` | `False` | Append a 6th return element with the per-draw accept/reject breakdown. Default keeps the 5-element return. |
| `diet_import_option` | `str` | `'as_DC'` | Passed to `SPPR_symbolic` when `kind='symbolic'`. |
| `silent` | `bool` | `True` | Suppress progress output. |
| `det_collapse_mode`, `det_open_mode`, `det_theta`, `det_external_sppr` | | | Forwarded to every SPPR call (see §6.5). |

**Why gamma sampling?** TE must stay strictly positive. The gamma distribution does this while letting you
specify mean and coefficient of variation directly: with `shape = 1/CV²` and `scale = TE_mean/shape`, the
draws have mean `TE_mean` and `std/mean = CV`. Basal rows are pinned to 1 each draw.

Only strictly positive TEs are sampled. `get_TE` can legitimately return `TE = 0` for a group (float
noise may make it a tiny negative), and gamma requires a positive scale, so those entries are pinned to
their model value for every draw and skipped by the clip band rather than sampled. With all TEs positive
this changes nothing. Such groups make `SPPR ~ 1/TE` near-singular, so their draws are typically rejected
as diverged — check `n_rejected_diverged` in the diagnostics and `diagnose_sppr`'s near-zero-TE warning.

**Returns** `(mean_sppr, accepted_samples_array, rejection_fraction, equations, variables)`. The last two
are `None` when `kind='new'`. `rejection_fraction` tells you how often the model was unstable — a high value
is a red flag about the model or the chosen TE error.
**Raises** `Exception` if `kind` is not `'new'` or `'symbolic'`.

---

## 7. Validation helpers

| Method | Returns | Meaning |
|--------|---------|---------|
| `is_model_balanced()` | `(bool, production, consumption)` | Checks the two mass-balance identities hold; returns the recomputed flow Series so you can inspect residuals. |
| `is_sppr_balanced(sppr, diet_import_equations=None)` | `(bool, inflow, outflow)` | Checks an SPPR result is globally self-consistent: primary-production **inflow** equals the export-weighted **outflow**. Pass `(equations, variables)` from `SPPR_symbolic` to use the diet-import-aware check. |
| `balance_model(change_production=False)` | `PPRCalculator` | Returns a deep-copied, exactly mass-balanced version of the model. `change_production=False` absorbs residuals into growth/net_migration (keeping `p` fixed); `True` absorbs them into production. |
| `diagnose_sppr(TE_option='GE', *, short=False, flat=False, thresholds=None, return_sppr=False, **sppr_kwargs)` | `dict` | One-call trust check for `SPPR_new` output: input-data checks, both convergence conditions, PP balance and the catch footprint, graded `OK` / `WARN` / `FAIL`. See §7.1. |

The global SPPR balance check verifies:

```
inflow  = total PP production entering the system
outflow = Σ (catch + growth + net_migration) · SPPR
is_balanced = isclose(inflow, outflow)
```

A balanced SPPR means the production you pulled out (as catch, growth, migration) is exactly the production
that came in as primary production — a strong sanity check on any method.

### 7.1 `diagnose_sppr` — can I rely on this result?

```python
report = model.diagnose_sppr()
```

`diagnose_sppr` runs `SPPR_new` once and grades everything that decides whether its output is
usable — the input data, the two convergence conditions, the global PP balance — plus the catch
footprint, which is reported but never graded. Use it before trusting any SPPR number from a model
you have not vetted, and to screen Monte-Carlo draws.

The verdict describes a **configuration, not a model**: `b` depends on `TE_option`, `det_theta` and
`det_open_mode`, so the same model can be healthy under one configuration and divergent under
another. The configuration evaluated is echoed back under `report['config']`. For the mathematics of
the two convergence conditions see `SPPR_Methods.md` §4 (`diagnose_sppr`).

**Full signature:**

```python
diagnose_sppr(TE_option='GE', *, short=False, flat=False, thresholds=None,
              return_sppr=False, **sppr_kwargs)
```

| Parameter | Type | Default | Meaning |
|-----------|------|---------|---------|
| `TE_option` | `str` | `'GE'` | As in `SPPR_new` (see §5.1). Moves both convergence numbers. Under `'TE'` there is no recycling matrix, so `b` is reported as `0.0` with a note. |
| `short` | `bool` | `False` | Return only `status`, `model_input`, `divergence`, `balance` — dropping `footprint`, `config` and `warnings`. |
| `flat` | `bool` | `False` | Return one-level `'<section>_<field>'` keys (`sppr_det` expanded per pool, lists replaced by counts) so many models concatenate into a DataFrame. |
| `thresholds` | `dict \| None` | `None` | Per-key overrides of `DEFAULT_DIAGNOSTIC_THRESHOLDS` (module level in `PPRCalculator.py`); missing keys keep their defaults. |
| `return_sppr` | `bool` | `False` | Also return `(SPPR, A, L)` from the internal solve, so you don't pay for the nullspace twice. |
| `**sppr_kwargs` | — | — | Forwarded verbatim to `SPPR_new`: `TE`, `DET_TE_vals`, and the four detritus knobs of §6.5. |

**Returns** a `dict` — three graded sections plus, unless `short=True`, the ungraded `footprint`, the
echoed `config` and a `warnings` list. Every field is documented below.
**Raises** `ValueError` only for an unrecognized `TE_option`, `det_open_mode` or
`det_collapse_mode` — caller mistakes. It never raises on a sick model: a singular `b = 1` solve or a
diverging pooled fallback is reported as `'FAIL'`, and a diagnostic re-solve with
`det_collapse_mode='never'` recovers `b` so the report still explains the failure.

**`report['status']`** — `'OK'`, `'WARN'` or `'FAIL'`: the worst of the three graded sections. Each
section also carries its own `status`. `footprint` is excluded from the verdict on purpose — a large
PPR/NPP is a finding about the ecosystem, not a defect in the calculation.

#### `report['model_input']` — the Ecopath data (independent of `TE_option`)

| Key | Type | Meaning |
|-----|------|---------|
| `status` | `str` | Section verdict. |
| `is_model_balanced` | `bool` | `is_model_balanced()`'s strict boolean. Reported for completeness — the grading uses the residuals below, since the boolean is all-or-nothing and most real models fail it. |
| `p_max_rel_residual` | `float` | Worst group's relative deviation between recomputed and stored **production**. Every SPPR method assumes `p`, `q`, `M0`, `predation` are mutually consistent. |
| `q_max_rel_residual` | `float` | The same for **consumption**. Graded via `model_balance_warn` / `model_balance_fail`. |
| `dc_rows_sum_to_1` | `bool` | Whether every consumer diet row (including `diet_import`) sums to 1 within `dc_row_tol`. Diet fractions must partition intake or `A = DC/TE` misweights every path through that consumer. |
| `dc_max_deviation` | `float` | Largest observed deviation from 1. `ModelData.validate_DC` already raises at load time at `1e-3`, so the default `dc_row_tol=1e-6` is what makes this informative. |
| `n_negative_catch` | `int` | `PPR = C·SPPR`, so a negative catch gives a negative footprint from a healthy SPPR. `WARN`. |
| `n_zero_catch` | `int` | Zero-catch groups are normal; counted, not penalised. |
| `total_catch` | `float` | Summed catch. |
| `has_catch` | `bool` | `False` makes every footprint number identically 0 — vacuous rather than wrong, so `WARN`, and the convergence diagnostics stay valid. |
| `has_ee_issues` | `bool` | `True` if any `EE = 0` or marginal-EE group exists. |
| `n_ee0`, `ee0_groups` | `int`, `list` | Groups whose production is entirely non-predatory death. Under `TE_option='TE'` their whole TE row is 0, severing them from the nullspace and leaking the PP they consumed (see `fix_EE_0_cases`, §6.5). |
| `n_ee_marginal`, `ee_marginal_groups` | `int`, `list` | Groups with `0 < EE < ee_marginal` (default `1e-3`), where `SPPR ~ 1/te` is near-singular. |
| `n_ee_gt_1` | `int` | Groups with `EE > 1` — the classic Ecopath over-consumption flag. |

#### `report['divergence']` — the solve itself

| Key | Type | Meaning |
|-----|------|---------|
| `status` | `str` | Section verdict. |
| `b` | `float` | Detritus recycling gain `ρ(diag(θ)·B)`: how much detritus one unit of detritus regenerates through the death→detritus→consumption→death loop (§6.6). Measured **before** the solve-vs-pool decision, so `det_collapse_mode` cannot mask it. `0.0` under `TE_option='TE'`, which has no recycling matrix. |
| `b_converges` | `bool` | `b < 1`. `False` means the returned SPPR is not a convergent sum. |
| `rho_living` | `float` | `ρ(A_LL)`, the same condition for the living block of `A = DC/TE` — predation cycles and cannibalism. |
| `living_converges` | `bool` | `rho_living < 1`. `False` corrupts the basis `B` is built from, so `b` becomes meaningless too. |
| `sppr_det` | `dict` | `{detritus_seq: sppr_det}` — the primary production charged per unit of each detritus pool. |
| `max_sppr_det` | `float` | The largest of those, graded against `sppr_det_warn` (default `10`). |
| `max_sppr_group` | `dict` | `{'seq', 'tl', 'sppr', 'inv_te'}` for the group carrying the largest total SPPR (`sppr.sum(axis=1)`). Trophic levels come from `get_TL(break_cycles=True, DET_as_PP=True)` — the convention the other SPPR-facing methods use — rather than the `model.TL` attribute, which is built with `break_cycles=False`. |
| `max_tl_group` | `dict` | The same record for the group at the top of the trophic ordering. `inv_te` is `1/te` read off the TE matrix actually in use, so it follows `TE_option` — `p/q` under `'GE'`, `(p/q)(1−M0/p)` under `'TE'`, or the draw itself when an explicit `TE` was passed — and is the per-step production amplification that group's edges carry (`None` where `te = 0`). SPPR rises with trophic depth, so the two records naming the same group is the expected picture; when they disagree, something other than trophic depth dominates the solution. Both are `None` if the solve produced no SPPR. |
| `n_negative_sources` | `int` | Basal-source **columns** containing a negative value. Counted per column, not per group, because a negative detritus column is masked in a group's row total by its positive PP columns. `≥ 1` is `FAIL`. |
| `expect_negatives` | `bool` | The prediction `b ≥ 1` — whether negatives *should* be present. |
| `near_singular_te` | `list` | Group seqs whose TE in the matrix **actually in use** is within `1e-3` of zero, so `SPPR ~ 1/te` blows up. Depends on `TE_option` and on any explicit `TE`, unlike the EE fields above. |
| `solve_error` | `str \| None` | The exception from a failed solve; the other fields are then recovered by the diagnostic re-solve. |

#### `report['balance']` — the PP budget of the result

| Key | Type | Meaning |
|-----|------|---------|
| `status` | `str` | Section verdict, from `rel_gap` against `balance_warn` / `balance_fail` (defaults 1% / 5%). |
| `inflow` | `float` | Primary production entering the system. |
| `outflow` | `float` | The SPPR-weighted export leaving it: `(catch + growth + net_migration)·SPPR`. |
| `rel_gap` | `float` | `\|outflow − inflow\| / \|inflow\|`. An SPPR that converged but does not close the PP budget is arithmetically fine and physically wrong. |
| `is_balanced` | `bool` | `is_sppr_balanced`'s strict boolean, alongside the gap it comes from. |

#### `report['footprint']` — reported, never graded

`get_NPP(only_inner=False)` raises `not implemented yet`, so NPP has one value and the second view is
`only_pp` rather than an outer/inner pair (see §5.3).

| Key | Type | Meaning |
|-----|------|---------|
| `ppr_all` | `float` | Total PPR over every basal source (`only_inner=False`). |
| `ppr_inner` | `float` | With `Import` columns dropped: within-system PPR. |
| `ppr_pp_only` | `float` | With `Import` **and** detritus dropped. The gap from `ppr_all` is how much of the footprint is routed through detritus. |
| `npp` | `float` | `get_NPP(only_inner=True)` — summed production of the `PP` groups. |
| `ppr2npp` | `float` | `ppr_inner / npp`: the %PPR appropriated by the catch. |
| `ppr2npp_pp_only` | `float` | `ppr_pp_only / npp`. A ratio above 1 means the catch needs more primary production than the system makes — a strong hint something is wrong even when every convergence test passed. |

#### `report['config']` and `report['warnings']`

`config` echoes `TE_option`, `det_open_mode`, `det_theta`, `det_external_sppr`,
`det_collapse_mode`, `explicit_TE` (whether a `TE` matrix was supplied), `method`
(`'single_detritus'` / `'multi_detritus'` / `'pooled_detritus_scaling'`), `would_pool`, and `model`
(name and year; `None` for toy or `from_dict` instances). The verdict is only meaningful together
with the configuration that produced it.

`warnings` is one plain-language string per tripped threshold, naming the quantity, its value and the
threshold crossed. The two near-divergence warnings also quote the magnitude being inflated — the worst
`sppr_det` for `b`, and the largest group SPPR with its seq and TL for `rho_living` — so the warning
shows the consequence, not just the ratio.

**Two convergence conditions, and why passing both is still not enough:**

| Condition | Meaning | Failure |
|-----------|---------|---------|
| `b < 1` | The detritus loop: how much detritus one unit of detritus regenerates (§6.6). | `b ≥ 1` — recycling runs away; the returned SPPR is not a convergent sum. |
| `rho_living < 1` | The same condition for the living block of `A = DC/TE`. | `rho_living ≥ 1` — corrupts the basis `B` is built from, so `b` becomes meaningless too. |

Neither is sufficient, because `1/(1−b)` is finite for every `b < 1` but not bounded: as `b`
approaches 1 the recycled contribution grows without limit while both tests still pass. That is why
closeness to 1 (`b_warn`, default `0.7`) and the resulting magnitude `max_sppr_det`
(`sppr_det_warn`, default `10`) are graded too. `sppr_det ≈ 1` is the healthy signature — one unit of
detritus costing about one unit of primary production — and values up to ~10 stay plausible for a
pool fed largely by consumer mortality several trophic steps up. Well beyond that the number is not
usable: North Sea 1991 under `'GE'` converges at `b = 0.89` with `sppr_det = 62.6`.

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

# 3b. Before trusting sppr_new, ask whether this configuration is sound.
#     return_sppr=True reuses the one solve, so this costs no extra nullspace.
report, sppr_new, A, L = model.diagnose_sppr(TE_option='GE', return_sppr=True)
print("status:", report['status'])                     # 'OK' | 'WARN' | 'FAIL'
print("b:", round(report['divergence']['b'], 4),        # detritus loop gain, must be < 1
      "| rho_living:", round(report['divergence']['rho_living'], 4))
print("sppr_det:", report['divergence']['sppr_det'])    # PP cost per unit of each pool
for w in report['warnings']:
    print("  !", w)

# 4. Detritus openness sensitivity: lose 30% of recycled detritus.
sppr_open, _, _ = model.SPPR_new(
    TE_option='GE', det_open_mode='recycling_loss', det_theta=0.7)

# 5. Turn SPPR into total PPR and the PPR/NPP ratio.
#    get_PPR always returns a 1-row DataFrame; collapse it for the scalar total.
print("total PPR:", model.get_PPR(sppr_new).sum(axis=1).sum())
print("PPR/NPP:", model.get_PPR2NPP_ratio(sppr_new))

# 6. Uncertainty bounds (robust mode), gating out draws that diverge as well as
#    draws that go negative, and asking for the per-draw breakdown.
mean, samples, reject_frac, _, _, diag = model.monte_carlo_SPPR(
    n_samples=200, TE_error_percent=10, kind='new',
    method_kwargs={'fix_EE_0_cases': True},
    exclude_diverged=True, return_diagnostics=True,
    det_collapse_mode='auto')
print(f"rejected {reject_frac:.0%} of unstable draws "
      f"({diag['n_rejected_diverged']} diverged, {diag['n_rejected_negative']} negative)")

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
- **Unsure whether a result is usable at all** → `diagnose_sppr()` first. It reports `b` and
  `rho_living` (both must be `< 1`), the per-pool `sppr_det`, and the PP balance gap. Do this before
  interpreting any number from a model you have not vetted — a model can converge cleanly and still be
  unusable.
