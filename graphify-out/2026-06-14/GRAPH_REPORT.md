# Graph Report - FishEstimationAI  (2026-06-14)

## Corpus Check
- 247 files · ~3,165,872 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 184 nodes · 292 edges · 32 communities (8 shown, 24 thin omitted)
- Extraction: 99% EXTRACTED · 1% INFERRED · 0% AMBIGUOUS · INFERRED: 2 edges (avg confidence: 0.5)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `be221734`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- [[_COMMUNITY_PPRCalculator Core Methods|PPRCalculator Core Methods]]
- [[_COMMUNITY_Community 1|Community 1]]
- [[_COMMUNITY_Community 2|Community 2]]
- [[_COMMUNITY_Ecopath Defaults & LIM|Ecopath Defaults & LIM]]
- [[_COMMUNITY_Excel Export & Entry Points|Excel Export & Entry Points]]
- [[_COMMUNITY_ModelData Species Groups|ModelData Species Groups]]
- [[_COMMUNITY_ModelData IO & JSON Loading|ModelData I/O & JSON Loading]]
- [[_COMMUNITY_Community 7|Community 7]]
- [[_COMMUNITY_SpeciesGroup Model|SpeciesGroup Model]]
- [[_COMMUNITY_Diet Composition Concept|Diet Composition Concept]]
- [[_COMMUNITY_Ecopath Marine Models|Ecopath Marine Models]]
- [[_COMMUNITY_Gross Efficiency Concept|Gross Efficiency Concept]]
- [[_COMMUNITY_Model Balance Validation|Model Balance Validation]]
- [[_COMMUNITY_PPRs All Sheet|PPRs All Sheet]]
- [[_COMMUNITY_PPRs Dataset|PPRs Dataset]]
- [[_COMMUNITY_PPRs Inner Sheet|PPRs Inner Sheet]]
- [[_COMMUNITY_SPPR 1986 Method|SPPR 1986 Method]]
- [[_COMMUNITY_SPPR 1995 mTL Method|SPPR 1995 mTL Method]]
- [[_COMMUNITY_SPPR 1995 TL2 Method|SPPR 1995 TL2 Method]]
- [[_COMMUNITY_SPPR 2015 Method|SPPR 2015 Method]]
- [[_COMMUNITY_SPPR EwE Method|SPPR EwE Method]]
- [[_COMMUNITY_SPPR Monte Carlo GE|SPPR Monte Carlo GE]]
- [[_COMMUNITY_SPPR New 2015 Method|SPPR New 2015 Method]]
- [[_COMMUNITY_SPPR New Full Method|SPPR New Full Method]]
- [[_COMMUNITY_SPPR New GE Method|SPPR New GE Method]]
- [[_COMMUNITY_SPPR Symbolic GE|SPPR Symbolic GE]]
- [[_COMMUNITY_Transfer Efficiency Concept|Transfer Efficiency Concept]]
- [[_COMMUNITY_Trophic Level Concept|Trophic Level Concept]]
- [[_COMMUNITY_Global PPR Estimation|Global PPR Estimation]]
- [[_COMMUNITY_PPR Estimation Project|PPR Estimation Project]]
- [[_COMMUNITY_Community 31|Community 31]]
- [[_COMMUNITY_Community 32|Community 32]]

## God Nodes (most connected - your core abstractions)
1. `PPRCalculator` - 46 edges
2. `ModelData` - 15 edges
3. `Multi-DET Openness, Collapse Modes & Diagnostics — Implementation Plan` - 14 edges
4. `FishEstimationAI` - 9 edges
5. `get_DC()` - 7 edges
6. `mat_from_np()` - 7 edges
7. `get_seq2name()` - 6 edges
8. `remove_cycles()` - 6 edges
9. `Tests executed` - 6 edges
10. `move_scattered_identity()` - 5 edges

## Surprising Connections (you probably didn't know these)
- `PPRCalculator` --uses--> `ModelData`  [INFERRED]
  PPRCalculator.py → ModelData.py
- `ModelData` --uses--> `ModelData`  [INFERRED]
  PPRCalculator.py → ModelData.py
- `main()` --calls--> `PPRCalculator`  [EXTRACTED]
  create_PPRS_excel.py → PPRCalculator.py

## Import Cycles
- 1-file cycle: `PPRCalculator.py -> PPRCalculator.py`

## Hyperedges (group relationships)
- **PPR Estimation Methods Comparison** — converted_pprs_59b295f7_sppr_1986, converted_pprs_59b295f7_sppr_1995_mtl_global_te_01, converted_pprs_59b295f7_sppr_1995_tl2_global_te_01, converted_pprs_59b295f7_sppr_1995_mtl_global_mte, converted_pprs_59b295f7_sppr_ewe, converted_pprs_59b295f7_sppr_2015, converted_pprs_59b295f7_sppr_new_2015, converted_pprs_59b295f7_sppr_new_full, converted_pprs_59b295f7_sppr_new_ge, converted_pprs_59b295f7_sppr_symbolic_ge, converted_pprs_59b295f7_sppr_symbolic_te, converted_pprs_59b295f7_sppr_mc_ge [EXTRACTED 1.00]

## Communities (32 total, 24 thin omitted)

### Community 0 - "PPRCalculator Core Methods"
Cohesion: 0.13
Nodes (13): PPRCalculator, Build the detritus recycling system (I - B) x = c for GE / With Egestion., Largest absolute eigenvalue of M (0 for empty). Used to test whether the, Resolve a per-DET parameter into an np.array aligned with DET_seq.          Ac, Fallback DET scaling: treat all DET groups as one pooled pool, solve the 1-D, Apply openness, decide solve-vs-pool by spectral radius / conditioning, scale th, Monte-Carlo uncertainty propagation over transfer efficiency. Repeatedly resampl, Variant of monte_carlo_SPPR supporting only kind='new'. Pre-allocates the sample (+5 more)

### Community 2 - "Community 2"
Cohesion: 0.13
Nodes (14): Multi-DET Openness, Collapse Modes & Diagnostics — Implementation Plan, New public parameter surface (identical across all methods), Self-review notes, Task 0: Capture full golden baseline (regression oracle), Task 10: Change log, cleanup, git, graphify, Task 1: Add `_resolve_det_param` and `_spectral_radius` helpers, Task 2: Add `_build_det_BC` (numeric recycling system builder), Task 3: Upgrade `_collapse_det_scaling` (points 7 + openness, return diagnostics) (+6 more)

### Community 3 - "Ecopath Defaults & LIM"
Cohesion: 0.19
Nodes (4): ModelData, Applies Ecopath defaults and ensures flows are synced with ratios.         Assu, Uses Linear Inverse Modeling (SLSQP) to fill in missing mass-balance          v, rebalance by changing growth, and net_migration (if change_production=False) or

### Community 4 - "Excel Export & Entry Points"
Cohesion: 0.18
Nodes (15): main(), Ido's matrix (nullspace) reformulation of the EwE path-summation SPPR: instead o, # TODO: change this function so I can decide which subset of parameters stays co, _find_all_cycles(), _get_circuit_probability(), mat_from_np(), move_scattered_identity(), Removes cycles from a flow matrix Z using the Ulanowicz method.     Z[i, j] rep (+7 more)

### Community 5 - "ModelData Species Groups"
Cohesion: 0.08
Nodes (21): get_DC(), get_seq2name(), load_json_dict(), ModelData, ModelData class supporting both old API (model_number) and new API (json_filepat, Initialize ModelData from either a model_number (int) or json_filepath (str)., Load a JSON file and return as dictionary., Initialize ModelData from a JSON filepath.                  Filename format: { (+13 more)

### Community 6 - "ModelData I/O & JSON Loading"
Cohesion: 0.17
Nodes (11): 1. `SPPR_new` golden regression (default params == old), 2. `SPPR_symbolic` golden regression (default params == old), all 12 combos, 3. Feature behaviour tests (`scratch_feature_tests.py`) — 11/11 passed, 4. Monte-Carlo smoke (negatives rejected, never raised), 5. Post-documentation re-verification, Files / methods touched (`PPRCalculator.py`), Multi-DET Openness, Collapse Modes & Diagnostics — Change Log, New public parameters (identical across `SPPR_new`, `SPPR_symbolic`, `monte_carlo_SPPR`, `monte_carlo_SPPR_2`) (+3 more)

### Community 31 - "Community 31"
Cohesion: 0.40
Nodes (4): FishEstimationAI — Claude Instructions, graphify, Knowledge Graph, Project Overview

### Community 32 - "Community 32"
Cohesion: 0.18
Nodes (10): Data, Dependencies, FishEstimationAI, Key Concepts, Notebooks, Overview, Programmatic use, Project Structure (+2 more)

## Knowledge Gaps
- **53 isolated node(s):** `Review points addressed`, `Openness math (per-DET θ, ext aligned to the DET columns)`, `Files / methods touched (`PPRCalculator.py`)`, `1. `SPPR_new` golden regression (default params == old)`, `2. `SPPR_symbolic` golden regression (default params == old), all 12 combos` (+48 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **24 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `PPRCalculator` connect `PPRCalculator Core Methods` to `Community 1`, `Ecopath Defaults & LIM`, `Excel Export & Entry Points`, `ModelData Species Groups`, `Community 7`?**
  _High betweenness centrality (0.262) - this node is a cross-community bridge._
- **Why does `ModelData` connect `ModelData Species Groups` to `PPRCalculator Core Methods`, `Ecopath Defaults & LIM`, `Excel Export & Entry Points`?**
  _High betweenness centrality (0.189) - this node is a cross-community bridge._
- **Why does `SpeciesGroupLegacy` connect `SpeciesGroup Model` to `ModelData Species Groups`?**
  _High betweenness centrality (0.028) - this node is a cross-community bridge._
- **Are the 2 inferred relationships involving `ModelData` (e.g. with `PPRCalculator` and `ModelData`) actually correct?**
  _`ModelData` has 2 INFERRED edges - model-reasoned connections that need verification._
- **What connects `Load a JSON file and return as dictionary.`, `Factory method to create an instance from a dictionary.`, `Factory method to create an instance from a dictionary.` to the rest of the system?**
  _91 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `PPRCalculator Core Methods` be split into smaller, more focused modules?**
  _Cohesion score 0.13124274099883856 - nodes in this community are weakly interconnected._
- **Should `Community 2` be split into smaller, more focused modules?**
  _Cohesion score 0.13333333333333333 - nodes in this community are weakly interconnected._