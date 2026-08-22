# Graph Report - FishEstimationAI  (2026-08-22)

## Corpus Check
- 265 files · ~3,236,662 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 479 nodes · 776 edges · 41 communities (12 shown, 29 thin omitted)
- Extraction: 97% EXTRACTED · 3% INFERRED · 0% AMBIGUOUS · INFERRED: 22 edges (avg confidence: 0.55)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `0e90e927`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- [[_COMMUNITY_SPPR Solver Core|SPPR Solver Core]]
- [[_COMMUNITY_Model Data Loading|Model Data Loading]]
- [[_COMMUNITY_Model Construction & Balancing|Model Construction & Balancing]]
- [[_COMMUNITY_User Guide|User Guide]]
- [[_COMMUNITY_Ecopath Concepts & Detritus|Ecopath Concepts & Detritus]]
- [[_COMMUNITY_Excel Export & Matrix Utils|Excel Export & Matrix Utils]]
- [[_COMMUNITY_PPR Concepts & Classic Methods|PPR Concepts & Classic Methods]]
- [[_COMMUNITY_README  Project Overview|README / Project Overview]]
- [[_COMMUNITY_Legacy Species Group|Legacy Species Group]]
- [[_COMMUNITY_Diet Composition Concept|Diet Composition Concept]]
- [[_COMMUNITY_Ecopath Models|Ecopath Models]]
- [[_COMMUNITY_Gross Efficiency Concept|Gross Efficiency Concept]]
- [[_COMMUNITY_Model Balance Validation|Model Balance Validation]]
- [[_COMMUNITY_PPRs All Sheet|PPRs All Sheet]]
- [[_COMMUNITY_PPRs Excel Dataset|PPRs Excel Dataset]]
- [[_COMMUNITY_PPRs Inner Sheet|PPRs Inner Sheet]]
- [[_COMMUNITY_SPPR 1986 Method|SPPR 1986 Method]]
- [[_COMMUNITY_SPPR 1995 mTL Method|SPPR 1995 mTL Method]]
- [[_COMMUNITY_SPPR 1995 TL2 Method|SPPR 1995 TL2 Method]]
- [[_COMMUNITY_SPPR 2015 Method|SPPR 2015 Method]]
- [[_COMMUNITY_SPPR EwE Method|SPPR EwE Method]]
- [[_COMMUNITY_SPPR Monte Carlo Method|SPPR Monte Carlo Method]]
- [[_COMMUNITY_SPPR New 2015 Method|SPPR New 2015 Method]]
- [[_COMMUNITY_SPPR New Full Method|SPPR New Full Method]]
- [[_COMMUNITY_SPPR New GE Method|SPPR New GE Method]]
- [[_COMMUNITY_SPPR Symbolic Methods|SPPR Symbolic Methods]]
- [[_COMMUNITY_Transfer Efficiency Concept|Transfer Efficiency Concept]]
- [[_COMMUNITY_Trophic Level Concept|Trophic Level Concept]]
- [[_COMMUNITY_Global PPR Estimation|Global PPR Estimation]]
- [[_COMMUNITY_PPR Estimation Project|PPR Estimation Project]]
- [[_COMMUNITY_Community 30|Community 30]]
- [[_COMMUNITY_Community 31|Community 31]]
- [[_COMMUNITY_Community 32|Community 32]]
- [[_COMMUNITY_Community 33|Community 33]]
- [[_COMMUNITY_Community 34|Community 34]]
- [[_COMMUNITY_Community 35|Community 35]]
- [[_COMMUNITY_Community 36|Community 36]]
- [[_COMMUNITY_Community 37|Community 37]]
- [[_COMMUNITY_Community 38|Community 38]]
- [[_COMMUNITY_Community 40|Community 40]]
- [[_COMMUNITY_Community 41|Community 41]]

## God Nodes (most connected - your core abstractions)
1. `PPRCalculator` - 64 edges
2. `ModelData` - 32 edges
3. `DataFrame` - 28 edges
4. `Series` - 16 edges
5. `DataFrame` - 14 edges
6. `build_model_tables()` - 14 edges
7. `PPRCalculator` - 13 edges
8. `4. Flow-network SPPR methods` - 11 edges
9. `User Guide: `ModelData` and `PPRCalculator`` - 11 edges
10. `MethodSpec` - 10 edges

## Surprising Connections (you probably didn't know these)
- `Project Overview` --references--> `ModelData`  [EXTRACTED]
  CLAUDE.md → ModelData.py
- `ModelData class (user guide)` --references--> `ModelData`  [EXTRACTED]
  USER_GUIDE.md → ModelData.py
- `Project Overview` --references--> `PPRCalculator`  [EXTRACTED]
  CLAUDE.md → PPRCalculator.py
- `PPRCalculator class (user guide)` --references--> `PPRCalculator`  [EXTRACTED]
  USER_GUIDE.md → PPRCalculator.py
- `MethodSpec` --uses--> `ModelData`  [INFERRED]
  create_PPRS_excel.py → ModelData.py

## Import Cycles
- 1-file cycle: `PPRCalculator.py -> PPRCalculator.py`
- 1-file cycle: `create_PPRS_excel.py -> create_PPRS_excel.py`
- 1-file cycle: `tests/test_create_pprs_excel.py -> tests/test_create_pprs_excel.py`
- 1-file cycle: `tests/test_diagnose_sppr.py -> tests/test_diagnose_sppr.py`
- 1-file cycle: `tests/test_monte_carlo_sppr.py -> tests/test_monte_carlo_sppr.py`

## Hyperedges (group relationships)
- **SPPR estimation method family** — user_guide_sppr_1986, user_guide_sppr_1995, user_guide_sppr_2015, user_guide_sppr_ewe, user_guide_sppr_new, user_guide_sppr_symbolic, user_guide_monte_carlo_sppr [EXTRACTED 1.00]
- **Detritus recycling stabilization mechanism** — user_guide_detritus_knobs, user_guide_coupled_detritus_solve, user_guide_det_fate [EXTRACTED 0.85]
- **Model completion pipeline (defaults to LIM to balance)** — user_guide_apply_ecopath_defaults, user_guide_apply_lim, user_guide_mass_balance [EXTRACTED 0.85]

## Communities (41 total, 29 thin omitted)

### Community 0 - "SPPR Solver Core"
Cohesion: 0.06
Nodes (48): ndarray, PPRCalculator, DataFrame, ModelData, Series, Core constructor used by __init__: build the calculator from a loaded ModelData., Return the diet-composition (DC) matrix, optionally redefining detritus rows., Return the flow matrix Z = DC * q (consumption-weighted diet), with DET rows red (+40 more)

### Community 1 - "Model Data Loading"
Cohesion: 0.06
Nodes (38): Any, get_DC(), get_model_data(), get_model_diet_data(), get_model_metadata(), get_seq2name(), load_json_dict(), ModelData (+30 more)

### Community 2 - "Model Construction & Balancing"
Cohesion: 0.07
Nodes (29): black_sea(), _calc(), PPRCalculator, Tests for the method_kwargs / exclude_diverged / return_diagnostics additions to, Keys the wrapper already controls must not be settable twice from two places., Omitting method_kwargs must behave exactly like passing an empty dict., With exclude_diverged=True, FAIL-divergence draws must be dropped from the avera, The gate may only ever remove draws, and its rejections must be labelled 'diverg (+21 more)

### Community 3 - "User Guide"
Cohesion: 0.06
Nodes (33): 10. Choosing a method (rules of thumb), 1. Background: what is PPR / SPPR?, 2. Key ecological quantities and group types, 3.1 Constructor, 3.2 Attributes available after construction, 3.3 Useful methods, 3. `ModelData` — loading a model, 4.1 Primary constructor (+25 more)

### Community 4 - "Ecopath Concepts & Detritus"
Cohesion: 0.06
Nodes (40): Documentation, FishEstimationAI — Claude Instructions, Knowledge Graph, Project Overview, DC (Diet Composition matrix), Ecopath mass-balance framework, Ecopath JSON model format, LIM (Linear Inverse Modeling) (+32 more)

### Community 5 - "Excel Export & Matrix Utils"
Cohesion: 0.08
Nodes (51): _autoformat(), build_footprint_table(), build_groups_table(), build_health_table(), build_mc_table(), build_model_tables(), build_notes_table(), build_sppr_tables() (+43 more)

### Community 6 - "PPR Concepts & Classic Methods"
Cohesion: 0.07
Nodes (9): PPRCalculator, Tests for create_PPRS_excel.py -- the per-model SPPR/PPR Excel exporter and its, SPPR_2015 resolves only PP and Import, so dropping Import already leaves PP alon, SPPR_1986 returns one un-attributed column, so PP cannot be separated out of it., all = every source, inner = drop Import, PP = drop Import and DET., test_method_without_detritus_columns_has_equal_inner_and_pp(), test_non_source_resolved_method_is_nan_in_the_pp_sheet_only(), test_sheets_hold_the_matching_source_sums() (+1 more)

### Community 7 - "README / Project Overview"
Cohesion: 0.05
Nodes (36): 1. What SPPR means, 2. Shared building blocks, 3. Trophic-level SPPR methods, 4. Flow-network SPPR methods, 5. From SPPR to ecosystem footprint, 6. Quick chooser, By goal, `DET_as_PP` and `normalize` (in `get_DC` / `get_Z`) (+28 more)

### Community 30 - "Community 30"
Cohesion: 0.25
Nodes (7): Batch export, Documentation, FishEstimationAI, Project structure, Setup, Usage, What it does

### Community 31 - "Community 31"
Cohesion: 0.19
Nodes (15): Matrix (nullspace) reformulation of the EwE path-summation SPPR.          Inst, # NOTE: trophic levels are resolved further down, after the flow vectors exist -, # TODO: change this function so I can decide which subset of parameters stays co, _find_all_cycles(), _get_circuit_probability(), mat_from_np(), move_scattered_identity(), Removes cycles from a flow matrix Z using the Ulanowicz method.     Z[i, j] rep (+7 more)

### Community 32 - "Community 32"
Cohesion: 0.13
Nodes (15): divergent_report(), _flat_te(), DataFrame, An explicit TE matrix overrides TE_option, and inv_te must follow it., Lower TE amplifies every path, so b must rise as TE falls., b is measured on diag(theta) @ B, so retention loss must lower it., det_collapse_mode is a remedy, not a diagnosis: b is measured pre-decision., A constant TE matrix at `value`, with detritus rows left at 1 (as SPPR_new expec (+7 more)

### Community 33 - "Community 33"
Cohesion: 0.47
Nodes (6): black_sea(), _calc(), PPRCalculator, test_negative_catch_is_warned(), test_zero_catch_is_warned_without_invalidating_divergence(), toy()

## Knowledge Gaps
- **87 isolated node(s):** `Documentation`, `What it does`, `Project structure`, `Batch export`, `Setup` (+82 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **29 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `PPRCalculator` connect `SPPR Solver Core` to `Community 32`, `Model Data Loading`, `Community 33`, `Model Construction & Balancing`, `Ecopath Concepts & Detritus`, `Excel Export & Matrix Utils`, `PPR Concepts & Classic Methods`, `Legacy Species Group`, `Community 31`?**
  _High betweenness centrality (0.405) - this node is a cross-community bridge._
- **Why does `ModelData` connect `Model Data Loading` to `SPPR Solver Core`, `Ecopath Concepts & Detritus`, `Excel Export & Matrix Utils`, `PPR Concepts & Classic Methods`, `Community 31`?**
  _High betweenness centrality (0.206) - this node is a cross-community bridge._
- **Why does `ModelData class (user guide)` connect `Ecopath Concepts & Detritus` to `Model Data Loading`?**
  _High betweenness centrality (0.052) - this node is a cross-community bridge._
- **Are the 9 inferred relationships involving `PPRCalculator` (e.g. with `MethodSpec` and `ModelTables`) actually correct?**
  _`PPRCalculator` has 9 INFERRED edges - model-reasoned connections that need verification._
- **Are the 10 inferred relationships involving `ModelData` (e.g. with `MethodSpec` and `ModelTables`) actually correct?**
  _`ModelData` has 10 INFERRED edges - model-reasoned connections that need verification._
- **Are the 2 inferred relationships involving `DataFrame` (e.g. with `ModelData` and `PPRCalculator`) actually correct?**
  _`DataFrame` has 2 INFERRED edges - model-reasoned connections that need verification._
- **What connects `Load a JSON file (relative to this module's directory) and return it as a dict.`, `Build a ``SpeciesGroupLegacy`` from a flat dict of field values.          Factor`, `Return a human-readable one-line summary of this legacy species group.` to the rest of the system?**
  _205 weakly-connected nodes found - possible documentation gaps or missing edges._