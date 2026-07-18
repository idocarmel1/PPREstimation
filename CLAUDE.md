# FishEstimationAI — Claude Instructions

## Knowledge Graph

This project has a pre-built knowledge graph in `graphify-out/graph.json` (144 nodes, 287 edges).

**Before answering any question about the codebase or making any code change, you MUST:**

1. Run `graphify query "<your question>"` to orient yourself using the graph.
2. Read only the specific files the graph points to — do not scan the whole project.

Use the saved Python interpreter:
```powershell
& (Get-Content graphify-out\.graphify_python) -m graphify query "<question>"
```

This applies to: architecture questions, "how does X work", "what calls Y", code reviews, refactoring requests, bug fixes — anything that requires understanding the codebase.

**Skip the graph only for:** pure file edits where the target file is already known and specified by the user.

## Project Overview

Fish biomass / PPR (Primary Production Required) estimation using Ecopath marine ecosystem models.
Given a mass-balanced food web, the code computes **SPPR** (Specific Primary Production Required)
per functional group and rolls it up into ecosystem-level footprints (total PPR, %PPR of NPP).

Key files (all SPPR logic lives in `PPRCalculator.py` — there are no separate per-method modules):
- `PPRCalculator.py` — the central computation class. Holds every SPPR method (`SPPR_1986`,
  `SPPR_1995`, `SPPR_1995_TL_fix`, `SPPR_EwE`, `SPPR_EwE_Ulanowicz`, `SPPR_2015`, `SPPR_new`,
  `SPPR_symbolic`), the `get_TE`/`get_TL`/`get_DC`/`get_Z` builders, the detritus recycling solver
  (`_build_det_BC`/`_solve_det_scaling`), Monte-Carlo uncertainty (`monte_carlo_SPPR`), and the
  footprint functions (`get_PPR`, `get_NPP`, `get_PPR2NPP_ratio`).
- `ModelData.py` — loads one Ecopath model from disk; exposes group parameters, the diet-composition
  (`DC`) and detritus-fate matrices, and name/seq lookups.
- `utils.py` — matrix utilities and the Ulanowicz `remove_cycles` cycle-removal algorithm.
- `create_PPRS_excel.py` — command-line entry point (`main()`); runs the methods and exports an
  Excel report to `output/`.
- `real_models/` — JSON model data: `EwE_jsons/`, `new_EwE_jsons/`, `ToyModels/`, `data/`.
- `notebooks/` — analysis and comparison notebooks (e.g. `comp0426`, `comp_paper`,
  `compare_spprs_pipeline`, `GS_distribution`, `tests`).
- `output/` — generated Excel / result files.

## Documentation

- `README.md` — short project overview and quickstart.
- `SPPR_Methods.md` — **the deep reference / source of truth** for every SPPR method, its flags, and
  the underlying math (TE options, flow-network nullspace, detritus recycling `(I−B)x=c`, Monte-Carlo
  Jensen correction, and the `get_PPR` / `get_PPR2NPP_ratio` footprint functions). Consult it before
  changing SPPR behavior.
- `USER_GUIDE.md` — practical guide to `ModelData` and `PPRCalculator` (constructors, attributes,
  every public method) with an end-to-end example; it references `SPPR_Methods.md` for the deep math.

When changing SPPR behavior, keep the code, `SPPR_Methods.md`, and `USER_GUIDE.md` in sync.

## graphify

This project has a knowledge graph at graphify-out/ with god nodes, community structure, and cross-file relationships.

Rules:
- For codebase questions, first run `graphify query "<question>"` when graphify-out/graph.json exists. Use `graphify path "<A>" "<B>"` for relationships and `graphify explain "<concept>"` for focused concepts. These return a scoped subgraph, usually much smaller than GRAPH_REPORT.md or raw grep output.
- If graphify-out/wiki/index.md exists, use it for broad navigation instead of raw source browsing.
- Read graphify-out/GRAPH_REPORT.md only for broad architecture review or when query/path/explain do not surface enough context.
- After modifying code, run `graphify update .` to keep the graph current (AST-only, no API cost).
