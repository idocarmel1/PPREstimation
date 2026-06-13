# Graph Report - .  (2026-06-13)

## Corpus Check
- 15 files · ~1,867,863 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 140 nodes · 257 edges · 31 communities (9 shown, 22 thin omitted)
- Extraction: 99% EXTRACTED · 1% INFERRED · 0% AMBIGUOUS · INFERRED: 2 edges (avg confidence: 0.5)
- Token cost: 0 input · 0 output

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

## God Nodes (most connected - your core abstractions)
1. `PPRCalculator` - 41 edges
2. `array` - 11 edges
3. `DataFrame` - 10 edges
4. `mat_from_np()` - 7 edges
5. `get_DC()` - 6 edges
6. `ModelData` - 6 edges
7. `remove_cycles()` - 6 edges
8. `remove_cycles()` - 6 edges
9. `SpeciesGroup` - 5 edges
10. `get_seq2name()` - 5 edges

## Surprising Connections (you probably didn't know these)
- `PPRCalculator` --uses--> `ModelData`  [INFERRED]
  PPRCalculator.py → ModelData.py
- `ModelData` --uses--> `ModelData`  [INFERRED]
  PPRCalculator.py → ModelData.py
- `main()` --calls--> `PPRCalculator`  [EXTRACTED]
  create_PPRS_excel.py → PPRCalculator.py
- `remove_cycles_new()` --calls--> `remove_cycles()`  [EXTRACTED]
  utils.py → remove_cycles_fix.py
- `SPPR_2015()` --calls--> `get_seq2name()`  [EXTRACTED]
  utils.py → ModelData.py

## Import Cycles
- 1-file cycle: `PPRCalculator.py -> PPRCalculator.py`

## Hyperedges (group relationships)
- **PPR Estimation Methods Comparison** — converted_pprs_59b295f7_sppr_1986, converted_pprs_59b295f7_sppr_1995_mtl_global_te_01, converted_pprs_59b295f7_sppr_1995_tl2_global_te_01, converted_pprs_59b295f7_sppr_1995_mtl_global_mte, converted_pprs_59b295f7_sppr_ewe, converted_pprs_59b295f7_sppr_2015, converted_pprs_59b295f7_sppr_new_2015, converted_pprs_59b295f7_sppr_new_full, converted_pprs_59b295f7_sppr_new_ge, converted_pprs_59b295f7_sppr_symbolic_ge, converted_pprs_59b295f7_sppr_symbolic_te, converted_pprs_59b295f7_sppr_mc_ge [EXTRACTED 1.00]

## Communities (31 total, 22 thin omitted)

### Community 0 - "PPRCalculator Core Methods"
Cohesion: 0.17
Nodes (6): PPRCalculator, Args:             n_samples (int, optional): number of sppr samples. Defaults t, rebalance by changing growth, and net_migration (if change_production=False) or, Args:             DET_as_PP (bool, optional): if True, DET row is set to 1. oth, get Z matrix. if DET_as_PP is False (default), DET row is flow_to_det. otherwise, Args:             global_TE (float, optional): 'mean' or float. Defaults to 0.1

### Community 1 - "SPPR 2015 Module"
Cohesion: 0.19
Nodes (8): decoder(), get_DC(), get_model_diet_data(), get_model_groups_data(), get_seq2name(), Factory method to create an instance from a dictionary., SpeciesGroup, SPPR_2015()

### Community 2 - "PPR Computation Methods"
Cohesion: 0.31
Nodes (14): array, DataFrame, method_1995_groups(), method_1995_species(), method_1995_TL_fix(), method_2015_full_structure_differential_TE_DET_from_pp(), method_2015_full_structure_TE_mc_DET_from_pp(), method_full_structure() (+6 more)

### Community 3 - "Ecopath Defaults & LIM"
Cohesion: 0.23
Nodes (3): ModelData, Applies Ecopath defaults and ensures flows are synced with ratios.         Assu, Uses Linear Inverse Modeling (SLSQP) to fill in missing mass-balance          v

### Community 4 - "Excel Export & Entry Points"
Cohesion: 0.31
Nodes (8): main(), # TODO: change this function so I can decide which subset of parameters stays co, mat_from_np(), Removes cycles from a flow matrix Z using the Ulanowicz method.     Z[i, j] rep, Convert numpy array to sympy Matrix with optional rational resolution., remove_cycles(), remove_cycles_new(), remove_cycles_old()

### Community 5 - "ModelData Species Groups"
Cohesion: 0.29
Nodes (7): get_DC(), get_seq2name(), ModelData, Get mapping from group sequence number to group name., Get diet composition (DC) and detritus fate matrices., get_model_groups_data(), SPPR_2015()

### Community 6 - "ModelData I/O & JSON Loading"
Cohesion: 0.25
Nodes (4): load_json_dict(), Load a JSON file and return as dictionary., Load JSON file and deserialize to list of SpeciesGroup objects., read_json_SpeciesGroup_list()

### Community 7 - "Cycle Detection & Removal"
Cohesion: 0.38
Nodes (6): find_all_cycles(), get_circuit_probability(), Calculates the probability of a quantum completing the cycle[cite: 97].     The, Phase 1: Enumerates all simple cycles using a depth-first      backtracking sea, Phase 2: Matrix decomposition by grouping cycles into nexuses      based on sha, remove_cycles()

### Community 9 - "EwE SPPR & Utilities"
Cohesion: 0.50
Nodes (3): Args:             TE_option (str): should be one of ['GE', 'TE', 'With Egestion, move_scattered_identity(), Identifies rows and columns that form an identity matrix,     even if they are

## Knowledge Gaps
- **21 isolated node(s):** `PPR Estimation Project`, `Global PPR Estimation`, `PPRs Dataset (Excel Spreadsheet)`, `SPPR 1986 Estimation Method`, `SPPR 1995 TL2 Global TE 0.1 Method` (+16 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **22 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `PPRCalculator` connect `PPRCalculator Core Methods` to `EwE SPPR & Utilities`, `Ecopath Defaults & LIM`, `Excel Export & Entry Points`, `ModelData Species Groups`?**
  _High betweenness centrality (0.274) - this node is a cross-community bridge._
- **Why does `ModelData` connect `ModelData Species Groups` to `PPRCalculator Core Methods`, `Ecopath Defaults & LIM`, `Excel Export & Entry Points`, `ModelData I/O & JSON Loading`?**
  _High betweenness centrality (0.108) - this node is a cross-community bridge._
- **Why does `remove_cycles()` connect `Cycle Detection & Removal` to `Excel Export & Entry Points`?**
  _High betweenness centrality (0.061) - this node is a cross-community bridge._
- **What connects `Load a JSON file and return as dictionary.`, `Factory method to create an instance from a dictionary.`, `Load JSON file and deserialize to list of SpeciesGroup objects.` to the rest of the system?**
  _41 weakly-connected nodes found - possible documentation gaps or missing edges._