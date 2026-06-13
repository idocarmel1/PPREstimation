# Graph Report - FishEstimationAI  (2026-06-14)

## Corpus Check
- 248 files · ~3,155,215 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 159 nodes · 269 edges · 29 communities (7 shown, 22 thin omitted)
- Extraction: 99% EXTRACTED · 1% INFERRED · 0% AMBIGUOUS · INFERRED: 4 edges (avg confidence: 0.5)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `48629c82`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- [[_COMMUNITY_PPRCalculator Core Methods|PPRCalculator Core Methods]]
- [[_COMMUNITY_Ecopath Defaults & LIM|Ecopath Defaults & LIM]]
- [[_COMMUNITY_Excel Export & Entry Points|Excel Export & Entry Points]]
- [[_COMMUNITY_ModelData Species Groups|ModelData Species Groups]]
- [[_COMMUNITY_ModelData IO & JSON Loading|ModelData I/O & JSON Loading]]
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
1. `PPRCalculator` - 49 edges
2. `ModelData` - 17 edges
3. `TestIcelandRegression` - 11 edges
4. `TestHumboldtLoads` - 8 edges
5. `get_DC()` - 7 edges
6. `mat_from_np()` - 7 edges
7. `get_seq2name()` - 6 edges
8. `remove_cycles()` - 6 edges
9. `SpeciesGroupLegacy` - 5 edges
10. `SpeciesGroup` - 5 edges

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

## Communities (29 total, 22 thin omitted)

### Community 0 - "PPRCalculator Core Methods"
Cohesion: 0.16
Nodes (8): PPRCalculator, Fallback DET scaling when the coupled (I-B) system diverges (spectral radius >=, Args:             n_samples (int, optional): number of sppr samples. Defaults t, Args:             n_samples (int, optional): number of sppr samples. Defaults t, Args:             DET_as_PP (bool, optional): if True, DET row is set to 1. oth, get Z matrix. if DET_as_PP is False (default), DET rows are flow_to_det split by, Args:             TE_option (str): should be one of ['GE', 'TE', 'With Egestion, Args:             global_TE (float, optional): 'mean' or float. Defaults to 0.1

### Community 3 - "Ecopath Defaults & LIM"
Cohesion: 0.19
Nodes (4): ModelData, Applies Ecopath defaults and ensures flows are synced with ratios.         Assu, Uses Linear Inverse Modeling (SLSQP) to fill in missing mass-balance          v, rebalance by changing growth, and net_migration (if change_production=False) or

### Community 4 - "Excel Export & Entry Points"
Cohesion: 0.18
Nodes (15): main(), # TODO: change this function so I can decide which subset of parameters stays co, Args:             TE_option (str): should be one of ['GE', 'TE', 'With Egestion, _find_all_cycles(), _get_circuit_probability(), mat_from_np(), move_scattered_identity(), Removes cycles from a flow matrix Z using the Ulanowicz method.     Z[i, j] rep (+7 more)

### Community 5 - "ModelData Species Groups"
Cohesion: 0.08
Nodes (22): get_DC(), get_seq2name(), load_json_dict(), ModelData, ModelData class supporting both old API (model_number) and new API (json_filepat, Initialize ModelData from either a model_number (int) or json_filepath (str)., Load a JSON file and return as dictionary., Initialize ModelData from a JSON filepath.                  Filename format: { (+14 more)

### Community 6 - "ModelData I/O & JSON Loading"
Cohesion: 0.11
Nodes (4): Iceland (single-DET) must produce identical results throughout all code changes., Humboldt (multi-DET) must load without exception once guard is removed., TestHumboldtLoads, TestIcelandRegression

### Community 8 - "SpeciesGroup Model"
Cohesion: 0.29
Nodes (3): Factory method to create an instance from a dictionary., SpeciesGroupLegacy, Comprehensive test of ModelData and PPRCalculator integration. Tests both old A

### Community 31 - "Community 31"
Cohesion: 0.40
Nodes (4): FishEstimationAI — Claude Instructions, graphify, Knowledge Graph, Project Overview

## Knowledge Gaps
- **25 isolated node(s):** `Knowledge Graph`, `Project Overview`, `graphify`, `PPREstimation`, `PPR Estimation Project` (+20 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **22 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `PPRCalculator` connect `PPRCalculator Core Methods` to `Ecopath Defaults & LIM`, `Excel Export & Entry Points`, `ModelData Species Groups`, `ModelData I/O & JSON Loading`, `SpeciesGroup Model`?**
  _High betweenness centrality (0.444) - this node is a cross-community bridge._
- **Why does `ModelData` connect `ModelData Species Groups` to `PPRCalculator Core Methods`, `SpeciesGroup Model`, `Ecopath Defaults & LIM`, `Excel Export & Entry Points`?**
  _High betweenness centrality (0.250) - this node is a cross-community bridge._
- **Why does `TestIcelandRegression` connect `ModelData I/O & JSON Loading` to `PPRCalculator Core Methods`?**
  _High betweenness centrality (0.081) - this node is a cross-community bridge._
- **Are the 3 inferred relationships involving `PPRCalculator` (e.g. with `ModelData` and `TestHumboldtLoads`) actually correct?**
  _`PPRCalculator` has 3 INFERRED edges - model-reasoned connections that need verification._
- **Are the 2 inferred relationships involving `ModelData` (e.g. with `PPRCalculator` and `ModelData`) actually correct?**
  _`ModelData` has 2 INFERRED edges - model-reasoned connections that need verification._
- **What connects `Load a JSON file and return as dictionary.`, `Factory method to create an instance from a dictionary.`, `Factory method to create an instance from a dictionary.` to the rest of the system?**
  _61 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `ModelData Species Groups` be split into smaller, more focused modules?**
  _Cohesion score 0.07539118065433854 - nodes in this community are weakly interconnected._