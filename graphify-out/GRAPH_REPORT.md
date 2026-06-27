# Graph Report - FishEstimationAI  (2026-06-27)

## Corpus Check
- 258 files · ~3,209,104 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 245 nodes · 416 edges · 35 communities (14 shown, 21 thin omitted)
- Extraction: 99% EXTRACTED · 1% INFERRED · 0% AMBIGUOUS · INFERRED: 5 edges (avg confidence: 0.5)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `0a3154ea`
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
- [[_COMMUNITY_Community 7|Community 7]]
- [[_COMMUNITY_Community 8|Community 8]]
- [[_COMMUNITY_Community 9|Community 9]]
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
- [[_COMMUNITY_README  Project Overview|README / Project Overview]]
- [[_COMMUNITY_Transfer Efficiency & SPPR Dispatch|Transfer Efficiency & SPPR Dispatch]]
- [[_COMMUNITY_Community 34|Community 34]]

## God Nodes (most connected - your core abstractions)
1. `PPRCalculator` - 47 edges
2. `DataFrame` - 28 edges
3. `ModelData` - 21 edges
4. `Series` - 14 edges
5. `User Guide: `ModelData` and `PPRCalculator`` - 11 edges
6. `FishEstimationAI` - 9 edges
7. `6. `PPRCalculator` — SPPR methods` - 9 edges
8. `DataFrame` - 8 edges
9. `get_DC()` - 8 edges
10. `Multi-DET TE discrepancy: `SPPR_new` vs `SPPR_symbolic(as_PP)`` - 8 edges

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

## Communities (35 total, 21 thin omitted)

### Community 0 - "PPRCalculator Core Methods"
Cohesion: 0.29
Nodes (4): Build a ``SpeciesGroupLegacy`` from a flat dict of field values.          Factor, Return a human-readable one-line summary of this legacy species group., Convert this legacy group into a single-row DataFrame of its scalar fields., SpeciesGroupLegacy

### Community 1 - "PPR / NPP Outputs & Balance Checks"
Cohesion: 0.21
Nodes (7): Series, Return the fraction of available NPP appropriated by the catch (PPR / NPP)., Unpack the fully-filled groups table into the individual named vectors., Relabel the index (and columns) of one or more SPPR-style results.          Ty, Check the two Ecopath mass-balance identities hold (within tolerance)., Compute the trophic level of every group via the standard linear-algebra definit, Convert a per-group SPPR into total primary production required (PPR) by the cat

### Community 2 - "Legacy Model Loading & Accessors"
Cohesion: 0.07
Nodes (27): 10. Choosing a method (rules of thumb), 1. Background: what is PPR / SPPR?, 2. Key ecological quantities and group types, 3.1 Constructor, 3.2 Attributes available after construction, 3.3 Useful methods, 3. `ModelData` — loading a model, 4.1 Primary constructor (+19 more)

### Community 3 - "Ecopath Defaults & LIM"
Cohesion: 0.17
Nodes (8): ModelData, Core constructor used by __init__: build the calculator from a loaded ModelData., Normalize the ordering of every Series/DataFrame attribute on the instance., Primary constructor: build the calculator directly from a model identifier., Apply standard Ecopath defaults and sync the mass-balance flows with the ratios., Fill missing mass-balance variables via a per-group Linear Inverse Model (SLSQP), Alternative constructor: rebuild an instance from a dict of pre-existing attribu, Return the underlying ModelData, or None for toy / from_dict instances.

### Community 4 - "Excel Export & Entry Points"
Cohesion: 0.18
Nodes (15): main(), Matrix (nullspace) reformulation of the EwE path-summation SPPR.          Inst, # TODO: change this function so I can decide which subset of parameters stays co, _find_all_cycles(), _get_circuit_probability(), mat_from_np(), move_scattered_identity(), Removes cycles from a flow matrix Z using the Ulanowicz method.     Z[i, j] rep (+7 more)

### Community 5 - "ModelData Species Groups"
Cohesion: 0.07
Nodes (34): Any, get_DC(), get_model_data(), get_model_diet_data(), get_model_metadata(), get_seq2name(), load_json_dict(), ModelData (+26 more)

### Community 6 - "ModelData I/O & JSON Loading"
Cohesion: 0.40
Nodes (4): FishEstimationAI — Claude Instructions, graphify, Knowledge Graph, Project Overview

### Community 7 - "Community 7"
Cohesion: 0.22
Nodes (8): 1. Symptom, 2. Root cause, 3. The design question (needs a decision before coding), 4. Suggested approach (if we pursue (b) — recommended), 5. Reproduction, 6. Relevant code (function names, not line numbers — lines drift), 7. Related context, Multi-DET TE discrepancy: `SPPR_new` vs `SPPR_symbolic(as_PP)`

### Community 8 - "Community 8"
Cohesion: 0.19
Nodes (8): PPRCalculator, 2015-method SPPR: a matrix-inversion (Leontief-style) formulation.          De, Run SPPR_new and force exact global balance by solving the single detritus SPPR., Iteratively force the model onto exact mass balance and return a balanced copy., Return a defensive (descending-seq sorted) copy of the per-group parameter table, Return the diet-composition (DC) matrix, optionally redefining detritus rows., Return the flow matrix Z = DC * q (consumption-weighted diet), with DET rows red, Return the sorted seq IDs of all detritus (DET) groups.          Returns:

### Community 9 - "Community 9"
Cohesion: 0.15
Nodes (8): Return the net primary production (NPP) of the system.          NPP is the tot, Symbolic SPPR helper, "diet import as DC" variant.          Imported diet is k, Symbolic SPPR solver: dispatch to the selected diet-import helper.          Ap, Monte-Carlo uncertainty propagation over transfer efficiency.          Repeate, Variant of monte_carlo_SPPR supporting only kind='new'.          Pre-allocates, Return the sorted seq IDs of all primary-producer (PP) groups.          Return, Return the sorted seq IDs of all regular (consumer) groups.          Returns:, Return the sorted seq IDs of all imported-diet (Import) groups.          Impor

### Community 31 - "Community 31"
Cohesion: 0.18
Nodes (7): DataFrame, Pauly & Christensen (1986)-style SPPR using a single catch-weighted trophic leve, Christensen & Pauly (1995)-style per-group SPPR = TE^(1-TL).          Uses eac, SPPR_1995 variant that linearly interpolates between bracketing integer trophic, Path-enumeration SPPR in the style of Ecopath with Ecosim (EwE) flow-network ana, Check an SPPR result is globally self-consistent (inflow == outflow)., Build the per-group transfer-efficiency (TE) vector or matrix.          Args:

### Community 32 - "README / Project Overview"
Cohesion: 0.18
Nodes (10): Data, Dependencies, FishEstimationAI, Key Concepts, Notebooks, Overview, Programmatic use, Project Structure (+2 more)

### Community 33 - "Transfer Efficiency & SPPR Dispatch"
Cohesion: 0.21
Nodes (8): ndarray, Build the detritus recycling system (I - B) x = c for the GE / With Egestion mod, Return the spectral radius (largest absolute eigenvalue) of M.          Used t, Resolve a per-DET parameter into a float array aligned with DET_seq., Fallback DET scaling: pool all detritus into one compartment and solve a 1-D pro, Resolve the detritus recycling system: apply openness, choose solve-vs-pool, sca, Primary numeric SPPR solver via the nullspace of L = A - I.          Builds th, Symbolic SPPR helper, "diet import as PP" variant.          Imported diet is t

### Community 34 - "Community 34"
Cohesion: 0.25
Nodes (7): 1. Symptom, 2. Root cause, 3. Suggested fix, 4. Reproduction, 5. Relevant code, 6. Related context, `SPPR_2015` bug: dropped-group block assignment leaks SPPR=1 into a zero-production basal column

## Knowledge Gaps
- **67 isolated node(s):** `Knowledge Graph`, `Project Overview`, `graphify`, `1. Symptom`, `2. Root cause` (+62 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **21 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `ModelData` connect `ModelData Species Groups` to `Transfer Efficiency & SPPR Dispatch`, `PPR / NPP Outputs & Balance Checks`, `Ecopath Defaults & LIM`, `Excel Export & Entry Points`, `Community 8`, `Community 31`?**
  _High betweenness centrality (0.208) - this node is a cross-community bridge._
- **Why does `PPRCalculator` connect `Community 8` to `Transfer Efficiency & SPPR Dispatch`, `PPR / NPP Outputs & Balance Checks`, `Ecopath Defaults & LIM`, `Excel Export & Entry Points`, `ModelData Species Groups`, `Community 9`, `Community 31`?**
  _High betweenness centrality (0.190) - this node is a cross-community bridge._
- **Why does `DataFrame` connect `Community 31` to `Transfer Efficiency & SPPR Dispatch`, `PPR / NPP Outputs & Balance Checks`, `Ecopath Defaults & LIM`, `Excel Export & Entry Points`, `ModelData Species Groups`, `Community 8`, `Community 9`?**
  _High betweenness centrality (0.062) - this node is a cross-community bridge._
- **Are the 5 inferred relationships involving `ModelData` (e.g. with `ndarray` and `PPRCalculator`) actually correct?**
  _`ModelData` has 5 INFERRED edges - model-reasoned connections that need verification._
- **What connects `Load a JSON file (relative to this module's directory) and return it as a dict.`, `Build a ``SpeciesGroupLegacy`` from a flat dict of field values.          Factor`, `Return a human-readable one-line summary of this legacy species group.` to the rest of the system?**
  _136 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Legacy Model Loading & Accessors` be split into smaller, more focused modules?**
  _Cohesion score 0.07142857142857142 - nodes in this community are weakly interconnected._
- **Should `ModelData Species Groups` be split into smaller, more focused modules?**
  _Cohesion score 0.06802721088435375 - nodes in this community are weakly interconnected._