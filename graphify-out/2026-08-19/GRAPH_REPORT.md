# Graph Report - FishEstimationAI  (2026-08-19)

## Corpus Check
- 263 files · ~3,225,339 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 344 nodes · 560 edges · 31 communities (10 shown, 21 thin omitted)
- Extraction: 98% EXTRACTED · 2% INFERRED · 0% AMBIGUOUS · INFERRED: 11 edges (avg confidence: 0.6)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `74d4bf12`
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

## God Nodes (most connected - your core abstractions)
1. `PPRCalculator` - 55 edges
2. `DataFrame` - 28 edges
3. `ModelData` - 23 edges
4. `Series` - 16 edges
5. `4. Flow-network SPPR methods` - 11 edges
6. `User Guide: `ModelData` and `PPRCalculator`` - 11 edges
7. `_flat_te()` - 9 edges
8. `6. `PPRCalculator` — SPPR methods` - 9 edges
9. `DataFrame` - 8 edges
10. `get_DC()` - 8 edges

## Surprising Connections (you probably didn't know these)
- `Project Overview` --references--> `ModelData`  [EXTRACTED]
  CLAUDE.md → ModelData.py
- `ModelData class (user guide)` --references--> `ModelData`  [EXTRACTED]
  USER_GUIDE.md → ModelData.py
- `Project Overview` --references--> `PPRCalculator`  [EXTRACTED]
  CLAUDE.md → PPRCalculator.py
- `PPRCalculator class (user guide)` --references--> `PPRCalculator`  [EXTRACTED]
  USER_GUIDE.md → PPRCalculator.py
- `ndarray` --uses--> `ModelData`  [INFERRED]
  PPRCalculator.py → ModelData.py

## Import Cycles
- 1-file cycle: `PPRCalculator.py -> PPRCalculator.py`
- 1-file cycle: `tests/test_diagnose_sppr.py -> tests/test_diagnose_sppr.py`

## Hyperedges (group relationships)
- **SPPR estimation method family** — user_guide_sppr_1986, user_guide_sppr_1995, user_guide_sppr_2015, user_guide_sppr_ewe, user_guide_sppr_new, user_guide_sppr_symbolic, user_guide_monte_carlo_sppr [EXTRACTED 1.00]
- **Detritus recycling stabilization mechanism** — user_guide_detritus_knobs, user_guide_coupled_detritus_solve, user_guide_det_fate [EXTRACTED 0.85]
- **Model completion pipeline (defaults to LIM to balance)** — user_guide_apply_ecopath_defaults, user_guide_apply_lim, user_guide_mass_balance [EXTRACTED 0.85]

## Communities (31 total, 21 thin omitted)

### Community 0 - "SPPR Solver Core"
Cohesion: 0.07
Nodes (40): ndarray, PPRCalculator, DataFrame, Series, Return the diet-composition (DC) matrix, optionally redefining detritus rows., Return the flow matrix Z = DC * q (consumption-weighted diet), with DET rows red, Return the sorted seq IDs of all detritus (DET) groups.          Returns:, Return the sorted seq IDs of all primary-producer (PP) groups.          Returns: (+32 more)

### Community 1 - "Model Data Loading"
Cohesion: 0.06
Nodes (38): Any, get_DC(), get_model_data(), get_model_diet_data(), get_model_metadata(), get_seq2name(), load_json_dict(), ModelData (+30 more)

### Community 2 - "Model Construction & Balancing"
Cohesion: 0.17
Nodes (8): ModelData, Core constructor used by __init__: build the calculator from a loaded ModelData., Normalize the ordering of every Series/DataFrame attribute on the instance., Apply standard Ecopath defaults and sync the mass-balance flows with the ratios., Primary constructor: build the calculator directly from a model identifier., Fill missing mass-balance variables via a per-group Linear Inverse Model (SLSQP), Return the underlying ModelData, or None for toy / from_dict instances., Alternative constructor: rebuild an instance from a dict of pre-existing attribu

### Community 3 - "User Guide"
Cohesion: 0.07
Nodes (28): 10. Choosing a method (rules of thumb), 1. Background: what is PPR / SPPR?, 2. Key ecological quantities and group types, 3.1 Constructor, 3.2 Attributes available after construction, 3.3 Useful methods, 3. `ModelData` — loading a model, 4.1 Primary constructor (+20 more)

### Community 4 - "Ecopath Concepts & Detritus"
Cohesion: 0.11
Nodes (25): DC (Diet Composition matrix), Ecopath mass-balance framework, Ecopath JSON model format, LIM (Linear Inverse Modeling), TL (Trophic Level), SPPR_EwE_Ulanowicz nullspace method, apply_ecopath_defaults (model completion), apply_lim (underdetermined solver) (+17 more)

### Community 5 - "Excel Export & Matrix Utils"
Cohesion: 0.18
Nodes (15): main(), Matrix (nullspace) reformulation of the EwE path-summation SPPR.          Instea, # TODO: change this function so I can decide which subset of parameters stays co, _find_all_cycles(), _get_circuit_probability(), mat_from_np(), move_scattered_identity(), Removes cycles from a flow matrix Z using the Ulanowicz method.     Z[i, j] rep (+7 more)

### Community 6 - "PPR Concepts & Classic Methods"
Cohesion: 0.16
Nodes (15): Documentation, FishEstimationAI — Claude Instructions, Knowledge Graph, Project Overview, NPP (Net Primary Production), Pauly (1986) citation, Pauly & Christensen (1995) citation, PPR (Primary Production Required) (+7 more)

### Community 7 - "README / Project Overview"
Cohesion: 0.06
Nodes (31): 1. What SPPR means, 2. Shared building blocks, 3. Trophic-level SPPR methods, 4. Flow-network SPPR methods, 5. From SPPR to ecosystem footprint, 6. Quick chooser, By goal, `DET_as_PP` and `normalize` (in `get_DC` / `get_Z`) (+23 more)

### Community 8 - "Legacy Species Group"
Cohesion: 0.06
Nodes (28): black_sea(), _calc(), divergent_report(), _flat_te(), DataFrame, PPRCalculator, Tests for PPRCalculator.diagnose_sppr.  Run with either:     pytest tests/test_d, The toy is mass-balanced, acyclic among living groups and far from divergence. (+20 more)

### Community 30 - "Community 30"
Cohesion: 0.29
Nodes (6): Documentation, FishEstimationAI, Project structure, Setup, Usage, What it does

## Knowledge Gaps
- **79 isolated node(s):** `Documentation`, `What it does`, `Project structure`, `Usage`, `Setup` (+74 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **21 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `PPRCalculator` connect `SPPR Solver Core` to `Model Data Loading`, `Model Construction & Balancing`, `Ecopath Concepts & Detritus`, `Excel Export & Matrix Utils`, `PPR Concepts & Classic Methods`, `Legacy Species Group`?**
  _High betweenness centrality (0.332) - this node is a cross-community bridge._
- **Why does `ModelData` connect `Model Data Loading` to `SPPR Solver Core`, `Model Construction & Balancing`, `Ecopath Concepts & Detritus`, `Excel Export & Matrix Utils`, `PPR Concepts & Classic Methods`?**
  _High betweenness centrality (0.228) - this node is a cross-community bridge._
- **Why does `ModelData class (user guide)` connect `Ecopath Concepts & Detritus` to `Model Data Loading`?**
  _High betweenness centrality (0.066) - this node is a cross-community bridge._
- **Are the 3 inferred relationships involving `PPRCalculator` (e.g. with `ModelData` and `DataFrame`) actually correct?**
  _`PPRCalculator` has 3 INFERRED edges - model-reasoned connections that need verification._
- **Are the 5 inferred relationships involving `ModelData` (e.g. with `ndarray` and `PPRCalculator`) actually correct?**
  _`ModelData` has 5 INFERRED edges - model-reasoned connections that need verification._
- **What connects `Load a JSON file (relative to this module's directory) and return it as a dict.`, `Build a ``SpeciesGroupLegacy`` from a flat dict of field values.          Factor`, `Return a human-readable one-line summary of this legacy species group.` to the rest of the system?**
  _159 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `SPPR Solver Core` be split into smaller, more focused modules?**
  _Cohesion score 0.0683526999316473 - nodes in this community are weakly interconnected._