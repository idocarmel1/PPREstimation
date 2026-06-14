# Graph Report - .  (2026-06-14)

## Corpus Check
- cluster-only mode — file stats not available

## Summary
- 188 nodes · 350 edges · 34 communities (11 shown, 23 thin omitted)
- Extraction: 99% EXTRACTED · 1% INFERRED · 0% AMBIGUOUS · INFERRED: 5 edges (avg confidence: 0.5)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `7a765b93`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- [[_COMMUNITY_PPRCalculator Core Methods|PPRCalculator Core Methods]]
- [[_COMMUNITY_PPR  NPP Outputs & Balance Checks|PPR / NPP Outputs & Balance Checks]]
- [[_COMMUNITY_Legacy Model Loading & Accessors|Legacy Model Loading & Accessors]]
- [[_COMMUNITY_Ecopath Defaults & LIM|Ecopath Defaults & LIM]]
- [[_COMMUNITY_Excel Export & Entry Points|Excel Export & Entry Points]]
- [[_COMMUNITY_ModelData Species Groups|ModelData Species Groups]]
- [[_COMMUNITY_ModelData IO & JSON Loading|ModelData I/O & JSON Loading]]
- [[_COMMUNITY_Monte-Carlo SPPR_2 & Import Groups|Monte-Carlo SPPR_2 & Import Groups]]
- [[_COMMUNITY_SpeciesGroup Model|SpeciesGroup Model]]
- [[_COMMUNITY_Symbolic SPPR & Detritus System Builder|Symbolic SPPR & Detritus System Builder]]
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
- [[_COMMUNITY_Forced-Balance SPPR & PP Groups|Forced-Balance SPPR & PP Groups]]
- [[_COMMUNITY_README  Project Overview|README / Project Overview]]
- [[_COMMUNITY_Transfer Efficiency & SPPR Dispatch|Transfer Efficiency & SPPR Dispatch]]

## God Nodes (most connected - your core abstractions)
1. `PPRCalculator` - 46 edges
2. `DataFrame` - 28 edges
3. `ModelData` - 18 edges
4. `Series` - 12 edges
5. `FishEstimationAI` - 9 edges
6. `get_DC()` - 8 edges
7. `ndarray` - 7 edges
8. `mat_from_np()` - 7 edges
9. `get_seq2name()` - 6 edges
10. `remove_cycles()` - 6 edges

## Surprising Connections (you probably didn't know these)
- `ndarray` --uses--> `ModelData`  [INFERRED]
  PPRCalculator.py → ModelData.py
- `PPRCalculator` --uses--> `ModelData`  [INFERRED]
  PPRCalculator.py → ModelData.py
- `DataFrame` --uses--> `ModelData`  [INFERRED]
  PPRCalculator.py → ModelData.py
- `ModelData` --uses--> `ModelData`  [INFERRED]
  PPRCalculator.py → ModelData.py
- `Series` --uses--> `ModelData`  [INFERRED]
  PPRCalculator.py → ModelData.py

## Import Cycles
- 1-file cycle: `PPRCalculator.py -> PPRCalculator.py`

## Hyperedges (group relationships)
- **PPR Estimation Methods Comparison** — converted_pprs_59b295f7_sppr_1986, converted_pprs_59b295f7_sppr_1995_mtl_global_te_01, converted_pprs_59b295f7_sppr_1995_tl2_global_te_01, converted_pprs_59b295f7_sppr_1995_mtl_global_mte, converted_pprs_59b295f7_sppr_ewe, converted_pprs_59b295f7_sppr_2015, converted_pprs_59b295f7_sppr_new_2015, converted_pprs_59b295f7_sppr_new_full, converted_pprs_59b295f7_sppr_new_ge, converted_pprs_59b295f7_sppr_symbolic_ge, converted_pprs_59b295f7_sppr_symbolic_te, converted_pprs_59b295f7_sppr_mc_ge [EXTRACTED 1.00]

## Communities (34 total, 23 thin omitted)

### Community 0 - "PPRCalculator Core Methods"
Cohesion: 0.27
Nodes (6): ndarray, Return the spectral radius (largest absolute eigenvalue) of M.          Used t, Resolve a per-DET parameter into a float array aligned with DET_seq., Fallback DET scaling: pool all detritus into one compartment and solve a 1-D pro, Resolve the detritus recycling system: apply openness, choose solve-vs-pool, sca, Primary numeric SPPR solver via the nullspace of L = A - I.          Builds th

### Community 1 - "PPR / NPP Outputs & Balance Checks"
Cohesion: 0.26
Nodes (7): PPRCalculator, Series, Relabel the index (and columns) of one or more SPPR-style results.          Ty, Check an SPPR result is globally self-consistent (inflow == outflow)., Convert a per-group SPPR into total primary production required (PPR) by the cat, Return the net primary production (NPP) of the system.          NPP is the tot, Return the fraction of available NPP appropriated by the catch (PPR / NPP).

### Community 2 - "Legacy Model Loading & Accessors"
Cohesion: 0.14
Nodes (14): Any, get_model_data(), get_model_diet_data(), get_model_metadata(), load_json_dict(), Load a JSON file (relative to this module's directory) and return it as a dict., Load a JSON file and deserialize every object into a ``SpeciesGroupLegacy``., Build a ``SpeciesGroupLegacy`` from a flat dict of field values.          Factor (+6 more)

### Community 3 - "Ecopath Defaults & LIM"
Cohesion: 0.12
Nodes (11): ModelData, Core constructor used by __init__: build the calculator from a loaded ModelData., Normalize the ordering of every Series/DataFrame attribute on the instance., Unpack the fully-filled groups table into the individual named vectors., Primary constructor: build the calculator directly from a model identifier., Apply standard Ecopath defaults and sync the mass-balance flows with the ratios., Fill missing mass-balance variables via a per-group Linear Inverse Model (SLSQP), Check the two Ecopath mass-balance identities hold (within tolerance). (+3 more)

### Community 4 - "Excel Export & Entry Points"
Cohesion: 0.18
Nodes (15): main(), Matrix (nullspace) reformulation of the EwE path-summation SPPR.          Inst, # TODO: change this function so I can decide which subset of parameters stays co, _find_all_cycles(), _get_circuit_probability(), mat_from_np(), move_scattered_identity(), Removes cycles from a flow matrix Z using the Ulanowicz method.     Z[i, j] rep (+7 more)

### Community 5 - "ModelData Species Groups"
Cohesion: 0.09
Nodes (21): get_DC(), get_seq2name(), ModelData, DataFrame, Series, Build a ``SpeciesGroup`` from a flat dict of field values.          Factory used, Convert this group into a single-row DataFrame of all its fields.          Unlik, Container for one Ecopath model's group table, diet composition and detritus fat (+13 more)

### Community 6 - "ModelData I/O & JSON Loading"
Cohesion: 0.25
Nodes (4): Pauly & Christensen (1986)-style SPPR using a single catch-weighted trophic leve, Christensen & Pauly (1995)-style per-group SPPR = TE^(1-TL).          Uses eac, SPPR_1995 variant that linearly interpolates between bracketing integer trophic, Compute the trophic level of every group via the standard linear-algebra definit

### Community 8 - "SpeciesGroup Model"
Cohesion: 0.32
Nodes (4): 2015-method SPPR: a matrix-inversion (Leontief-style) formulation.          De, Return the diet-composition (DC) matrix, optionally redefining detritus rows., Return the flow matrix Z = DC * q (consumption-weighted diet), with DET rows red, Return the sorted seq IDs of all detritus (DET) groups.          Returns:

### Community 9 - "Symbolic SPPR & Detritus System Builder"
Cohesion: 0.29
Nodes (4): Build the detritus recycling system (I - B) x = c for the GE / With Egestion mod, Symbolic SPPR helper, "diet import as PP" variant.          Imported diet is t, Symbolic SPPR helper, "diet import as DC" variant.          Imported diet is k, Return the sorted seq IDs of all regular (consumer) groups.          Returns:

### Community 32 - "README / Project Overview"
Cohesion: 0.18
Nodes (10): Data, Dependencies, FishEstimationAI, Key Concepts, Notebooks, Overview, Programmatic use, Project Structure (+2 more)

### Community 33 - "Transfer Efficiency & SPPR Dispatch"
Cohesion: 0.24
Nodes (6): DataFrame, Path-enumeration SPPR in the style of Ecopath with Ecosim (EwE) flow-network ana, Symbolic SPPR solver: dispatch to the selected diet-import helper.          Ap, Monte-Carlo uncertainty propagation over transfer efficiency.          Repeate, Return a defensive (descending-seq sorted) copy of the per-group parameter table, Build the per-group transfer-efficiency (TE) vector or matrix.          Args:

## Knowledge Gaps
- **30 isolated node(s):** `Series`, `Overview`, `Project Structure`, `SPPR Methods`, `Data` (+25 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **23 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `ModelData` connect `ModelData Species Groups` to `PPRCalculator Core Methods`, `PPR / NPP Outputs & Balance Checks`, `Legacy Model Loading & Accessors`, `Transfer Efficiency & SPPR Dispatch`, `Excel Export & Entry Points`, `Ecopath Defaults & LIM`?**
  _High betweenness centrality (0.307) - this node is a cross-community bridge._
- **Why does `PPRCalculator` connect `PPR / NPP Outputs & Balance Checks` to `PPRCalculator Core Methods`, `Transfer Efficiency & SPPR Dispatch`, `Ecopath Defaults & LIM`, `Excel Export & Entry Points`, `ModelData Species Groups`, `ModelData I/O & JSON Loading`, `Monte-Carlo SPPR_2 & Import Groups`, `SpeciesGroup Model`, `Symbolic SPPR & Detritus System Builder`, `Forced-Balance SPPR & PP Groups`?**
  _High betweenness centrality (0.302) - this node is a cross-community bridge._
- **Why does `DataFrame` connect `Transfer Efficiency & SPPR Dispatch` to `PPRCalculator Core Methods`, `PPR / NPP Outputs & Balance Checks`, `Ecopath Defaults & LIM`, `Excel Export & Entry Points`, `ModelData Species Groups`, `ModelData I/O & JSON Loading`, `Monte-Carlo SPPR_2 & Import Groups`, `SpeciesGroup Model`, `Symbolic SPPR & Detritus System Builder`, `Forced-Balance SPPR & PP Groups`?**
  _High betweenness centrality (0.101) - this node is a cross-community bridge._
- **Are the 5 inferred relationships involving `ModelData` (e.g. with `ndarray` and `PPRCalculator`) actually correct?**
  _`ModelData` has 5 INFERRED edges - model-reasoned connections that need verification._
- **What connects `Series`, `Load a JSON file (relative to this module's directory) and return it as a dict.`, `Build a ``SpeciesGroupLegacy`` from a flat dict of field values.          Factor` to the rest of the system?**
  _99 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Legacy Model Loading & Accessors` be split into smaller, more focused modules?**
  _Cohesion score 0.13970588235294118 - nodes in this community are weakly interconnected._
- **Should `Ecopath Defaults & LIM` be split into smaller, more focused modules?**
  _Cohesion score 0.11904761904761904 - nodes in this community are weakly interconnected._