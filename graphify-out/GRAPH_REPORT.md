# Graph Report - FishEstimationAI  (2026-06-13)

## Corpus Check
- 251 files · ~3,183,292 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 203 nodes · 319 edges · 33 communities (10 shown, 23 thin omitted)
- Extraction: 99% EXTRACTED · 1% INFERRED · 0% AMBIGUOUS · INFERRED: 4 edges (avg confidence: 0.5)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `1ea48fde`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- [[_COMMUNITY_PPRCalculator Core Methods|PPRCalculator Core Methods]]
- [[_COMMUNITY_SPPR 2015 Module|SPPR 2015 Module]]
- [[_COMMUNITY_PPR Computation Methods|PPR Computation Methods]]
- [[_COMMUNITY_Ecopath Defaults & LIM|Ecopath Defaults & LIM]]
- [[_COMMUNITY_Excel Export & Entry Points|Excel Export & Entry Points]]
- [[_COMMUNITY_ModelData Species Groups|ModelData Species Groups]]
- [[_COMMUNITY_ModelData IO & JSON Loading|ModelData I/O & JSON Loading]]
- [[_COMMUNITY_Cycle Detection & Removal|Cycle Detection & Removal]]
- [[_COMMUNITY_SpeciesGroup Model|SpeciesGroup Model]]
- [[_COMMUNITY_EwE SPPR & Utilities|EwE SPPR & Utilities]]
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
1. `PPRCalculator` - 48 edges
2. `ModelData` - 17 edges
3. `Multi-DET Group Support Implementation Plan` - 16 edges
4. `TestIcelandRegression` - 11 edges
5. `get_DC()` - 8 edges
6. `TestHumboldtLoads` - 8 edges
7. `get_seq2name()` - 7 edges
8. `mat_from_np()` - 7 edges
9. `remove_cycles()` - 6 edges
10. `remove_cycles()` - 6 edges

## Surprising Connections (you probably didn't know these)
- `PPRCalculator` --uses--> `ModelData`  [INFERRED]
  PPRCalculator.py → ModelData.py
- `ModelData` --uses--> `ModelData`  [INFERRED]
  PPRCalculator.py → ModelData.py
- `main()` --calls--> `PPRCalculator`  [EXTRACTED]
  create_PPRS_excel.py → PPRCalculator.py
- `TestHumboldtLoads` --uses--> `PPRCalculator`  [INFERRED]
  tests/test_multi_det.py → PPRCalculator.py
- `TestIcelandRegression` --uses--> `PPRCalculator`  [INFERRED]
  tests/test_multi_det.py → PPRCalculator.py

## Import Cycles
- 1-file cycle: `PPRCalculator.py -> PPRCalculator.py`

## Hyperedges (group relationships)
- **PPR Estimation Methods Comparison** — converted_pprs_59b295f7_sppr_1986, converted_pprs_59b295f7_sppr_1995_mtl_global_te_01, converted_pprs_59b295f7_sppr_1995_tl2_global_te_01, converted_pprs_59b295f7_sppr_1995_mtl_global_mte, converted_pprs_59b295f7_sppr_ewe, converted_pprs_59b295f7_sppr_2015, converted_pprs_59b295f7_sppr_new_2015, converted_pprs_59b295f7_sppr_new_full, converted_pprs_59b295f7_sppr_new_ge, converted_pprs_59b295f7_sppr_symbolic_ge, converted_pprs_59b295f7_sppr_symbolic_te, converted_pprs_59b295f7_sppr_mc_ge [EXTRACTED 1.00]

## Communities (33 total, 23 thin omitted)

### Community 0 - "PPRCalculator Core Methods"
Cohesion: 0.17
Nodes (7): PPRCalculator, Args:             n_samples (int, optional): number of sppr samples. Defaults t, Args:             n_samples (int, optional): number of sppr samples. Defaults t, Args:             DET_as_PP (bool, optional): if True, DET row is set to 1. oth, get Z matrix. if DET_as_PP is False (default), DET rows are flow_to_det split by, Args:             TE_option (str): should be one of ['GE', 'TE', 'With Egestion, Args:             global_TE (float, optional): 'mean' or float. Defaults to 0.1

### Community 1 - "SPPR 2015 Module"
Cohesion: 0.19
Nodes (8): decoder(), get_DC(), get_model_diet_data(), get_model_groups_data(), get_seq2name(), Factory method to create an instance from a dictionary., SpeciesGroup, SPPR_2015()

### Community 2 - "PPR Computation Methods"
Cohesion: 0.11
Nodes (18): Background: det_fate matrix semantics, Baseline test (must pass before and after every task), Known risks, Multi-DET Group Support Implementation Plan, Self-Review, Spec coverage, Task 0: Create test file and record Iceland baseline, Task 10: Fix `calc_2015/SPPR_2015.py` and `utils.py` (+10 more)

### Community 3 - "Ecopath Defaults & LIM"
Cohesion: 0.19
Nodes (4): ModelData, Applies Ecopath defaults and ensures flows are synced with ratios.         Assu, Uses Linear Inverse Modeling (SLSQP) to fill in missing mass-balance          v, rebalance by changing growth, and net_migration (if change_production=False) or

### Community 4 - "Excel Export & Entry Points"
Cohesion: 0.10
Nodes (14): main(), Debug script to find where the processing differs between APIs, # TODO: change this function so I can decide which subset of parameters stays co, Args:             TE_option (str): should be one of ['GE', 'TE', 'With Egestion, Humboldt (multi-DET) must load without exception once guard is removed., TestHumboldtLoads, mat_from_np(), move_scattered_identity() (+6 more)

### Community 5 - "ModelData Species Groups"
Cohesion: 0.08
Nodes (23): get_DC(), get_seq2name(), load_json_dict(), ModelData, ModelData class supporting both old API (model_number) and new API (json_filepat, Initialize ModelData from either a model_number (int) or json_filepath (str)., Load a JSON file and return as dictionary., Initialize ModelData from a JSON filepath.                  Filename format: { (+15 more)

### Community 7 - "Cycle Detection & Removal"
Cohesion: 0.38
Nodes (6): find_all_cycles(), get_circuit_probability(), Calculates the probability of a quantum completing the cycle[cite: 97].     The, Phase 1: Enumerates all simple cycles using a depth-first      backtracking sea, Phase 2: Matrix decomposition by grouping cycles into nexuses      based on sha, remove_cycles()

### Community 8 - "SpeciesGroup Model"
Cohesion: 0.29
Nodes (3): Factory method to create an instance from a dictionary., SpeciesGroupLegacy, Comprehensive test of ModelData and PPRCalculator integration. Tests both old A

### Community 9 - "EwE SPPR & Utilities"
Cohesion: 0.33
Nodes (5): Sheet: outflow minus inflow, Sheet: PPRs all, Sheet: PPRs inner, Sheet: SPPR did not explode, Sheet: SPPR is balanced

### Community 31 - "Community 31"
Cohesion: 0.40
Nodes (4): FishEstimationAI — Claude Instructions, graphify, Knowledge Graph, Project Overview

## Knowledge Gaps
- **46 isolated node(s):** `Knowledge Graph`, `Project Overview`, `graphify`, `PPREstimation`, `Background: det_fate matrix semantics` (+41 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **23 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `PPRCalculator` connect `PPRCalculator Core Methods` to `Ecopath Defaults & LIM`, `Excel Export & Entry Points`, `ModelData Species Groups`, `ModelData I/O & JSON Loading`, `SpeciesGroup Model`?**
  _High betweenness centrality (0.272) - this node is a cross-community bridge._
- **Why does `ModelData` connect `ModelData Species Groups` to `PPRCalculator Core Methods`, `SpeciesGroup Model`, `Ecopath Defaults & LIM`, `Excel Export & Entry Points`?**
  _High betweenness centrality (0.158) - this node is a cross-community bridge._
- **Why does `TestIcelandRegression` connect `ModelData I/O & JSON Loading` to `PPRCalculator Core Methods`, `Excel Export & Entry Points`?**
  _High betweenness centrality (0.051) - this node is a cross-community bridge._
- **Are the 3 inferred relationships involving `PPRCalculator` (e.g. with `ModelData` and `TestHumboldtLoads`) actually correct?**
  _`PPRCalculator` has 3 INFERRED edges - model-reasoned connections that need verification._
- **Are the 2 inferred relationships involving `ModelData` (e.g. with `PPRCalculator` and `ModelData`) actually correct?**
  _`ModelData` has 2 INFERRED edges - model-reasoned connections that need verification._
- **What connects `Load a JSON file and return as dictionary.`, `Factory method to create an instance from a dictionary.`, `Factory method to create an instance from a dictionary.` to the rest of the system?**
  _84 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `PPR Computation Methods` be split into smaller, more focused modules?**
  _Cohesion score 0.10526315789473684 - nodes in this community are weakly interconnected._