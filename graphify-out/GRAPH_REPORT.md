# Graph Report - .  (2026-06-27)

## Corpus Check
- 251 files · ~3,207,559 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 262 nodes · 455 edges · 30 communities (9 shown, 21 thin omitted)
- Extraction: 98% EXTRACTED · 2% INFERRED · 0% AMBIGUOUS · INFERRED: 9 edges (avg confidence: 0.62)
- Token cost: 49,559 input · 0 output

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

## God Nodes (most connected - your core abstractions)
1. `PPRCalculator` - 49 edges
2. `DataFrame` - 28 edges
3. `ModelData` - 23 edges
4. `Series` - 14 edges
5. `User Guide: `ModelData` and `PPRCalculator`` - 11 edges
6. `FishEstimationAI` - 9 edges
7. `6. `PPRCalculator` — SPPR methods` - 9 edges
8. `DataFrame` - 8 edges
9. `get_DC()` - 8 edges
10. `SPPR (Specific Primary Production Required)` - 8 edges

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

## Hyperedges (group relationships)
- **SPPR estimation method family** — user_guide_sppr_1986, user_guide_sppr_1995, user_guide_sppr_2015, user_guide_sppr_ewe, user_guide_sppr_new, user_guide_sppr_symbolic, user_guide_monte_carlo_sppr [EXTRACTED 1.00]
- **Detritus recycling stabilization mechanism** — user_guide_detritus_knobs, user_guide_coupled_detritus_solve, user_guide_det_fate [EXTRACTED 0.85]
- **Model completion pipeline (defaults to LIM to balance)** — user_guide_apply_ecopath_defaults, user_guide_apply_lim, user_guide_mass_balance [EXTRACTED 0.85]

## Communities (30 total, 21 thin omitted)

### Community 0 - "SPPR Solver Core"
Cohesion: 0.09
Nodes (26): ndarray, DataFrame, Series, Return the net primary production (NPP) of the system.          NPP is the tot, Return the fraction of available NPP appropriated by the catch (PPR / NPP)., 2015-method SPPR: a matrix-inversion (Leontief-style) formulation.          De, Build the detritus recycling system (I - B) x = c for the GE / With Egestion mod, Return the spectral radius (largest absolute eigenvalue) of M.          Used t (+18 more)

### Community 1 - "Model Data Loading"
Cohesion: 0.07
Nodes (33): Any, get_DC(), get_model_data(), get_model_diet_data(), get_model_metadata(), get_seq2name(), load_json_dict(), ModelData (+25 more)

### Community 2 - "Model Construction & Balancing"
Cohesion: 0.09
Nodes (19): PPRCalculator, ModelData, Pauly & Christensen (1986)-style SPPR using a single catch-weighted trophic leve, Core constructor used by __init__: build the calculator from a loaded ModelData., Christensen & Pauly (1995)-style per-group SPPR = TE^(1-TL).          Uses eac, SPPR_1995 variant that linearly interpolates between bracketing integer trophic, Path-enumeration SPPR in the style of Ecopath with Ecosim (EwE) flow-network ana, Normalize the ordering of every Series/DataFrame attribute on the instance. (+11 more)

### Community 3 - "User Guide"
Cohesion: 0.07
Nodes (27): 10. Choosing a method (rules of thumb), 1. Background: what is PPR / SPPR?, 2. Key ecological quantities and group types, 3.1 Constructor, 3.2 Attributes available after construction, 3.3 Useful methods, 3. `ModelData` — loading a model, 4.1 Primary constructor (+19 more)

### Community 4 - "Ecopath Concepts & Detritus"
Cohesion: 0.10
Nodes (27): DC (Diet Composition matrix), Ecopath mass-balance framework, Ecopath JSON model format, LIM (Linear Inverse Modeling), TE / GE (Transfer / Gross Efficiency), TL (Trophic Level), apply_ecopath_defaults (model completion), apply_lim (underdetermined solver) (+19 more)

### Community 5 - "Excel Export & Matrix Utils"
Cohesion: 0.18
Nodes (15): main(), Matrix (nullspace) reformulation of the EwE path-summation SPPR.          Inst, # TODO: change this function so I can decide which subset of parameters stays co, _find_all_cycles(), _get_circuit_probability(), mat_from_np(), move_scattered_identity(), Removes cycles from a flow matrix Z using the Ulanowicz method.     Z[i, j] rep (+7 more)

### Community 6 - "PPR Concepts & Classic Methods"
Cohesion: 0.19
Nodes (13): FishEstimationAI — Claude Instructions, graphify, Knowledge Graph, Project Overview, NPP (Net Primary Production), Pauly (1986) citation, Pauly & Christensen (1995) citation, PPR (Primary Production Required) (+5 more)

### Community 7 - "README / Project Overview"
Cohesion: 0.18
Nodes (10): Data, Dependencies, FishEstimationAI, Key Concepts, Notebooks, Overview, Programmatic use, Project Structure (+2 more)

### Community 8 - "Legacy Species Group"
Cohesion: 0.29
Nodes (4): Build a ``SpeciesGroupLegacy`` from a flat dict of field values.          Factor, Return a human-readable one-line summary of this legacy species group., Convert this legacy group into a single-row DataFrame of its scalar fields., SpeciesGroupLegacy

## Knowledge Gaps
- **58 isolated node(s):** `graphify`, `Overview`, `Project Structure`, `SPPR Methods`, `Data` (+53 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **21 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `ModelData` connect `Model Data Loading` to `SPPR Solver Core`, `Model Construction & Balancing`, `Ecopath Concepts & Detritus`, `Excel Export & Matrix Utils`, `PPR Concepts & Classic Methods`?**
  _High betweenness centrality (0.291) - this node is a cross-community bridge._
- **Why does `PPRCalculator` connect `Model Construction & Balancing` to `SPPR Solver Core`, `Model Data Loading`, `Ecopath Concepts & Detritus`, `Excel Export & Matrix Utils`, `PPR Concepts & Classic Methods`?**
  _High betweenness centrality (0.265) - this node is a cross-community bridge._
- **Why does `ModelData class (user guide)` connect `Ecopath Concepts & Detritus` to `Model Data Loading`?**
  _High betweenness centrality (0.090) - this node is a cross-community bridge._
- **Are the 5 inferred relationships involving `ModelData` (e.g. with `ndarray` and `PPRCalculator`) actually correct?**
  _`ModelData` has 5 INFERRED edges - model-reasoned connections that need verification._
- **What connects `Load a JSON file (relative to this module's directory) and return it as a dict.`, `Build a ``SpeciesGroupLegacy`` from a flat dict of field values.          Factor`, `Return a human-readable one-line summary of this legacy species group.` to the rest of the system?**
  _123 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `SPPR Solver Core` be split into smaller, more focused modules?**
  _Cohesion score 0.09142857142857143 - nodes in this community are weakly interconnected._
- **Should `Model Data Loading` be split into smaller, more focused modules?**
  _Cohesion score 0.07092198581560284 - nodes in this community are weakly interconnected._