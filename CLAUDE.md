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

Key files:
- `PPRCalculator.py` — central hub (41 edges), the main computation class
- `ModelData.py` — data loading and species group management
- `PPR_methods.py` — SPPR estimation method implementations (1986, 1995, 2015, EwE, symbolic)
- `remove_cycles_fix.py` — Ulanowicz cycle removal algorithm
- `utils.py` — matrix utilities
- `calc_2015/SPPR_2015.py` — 2015 method module
- `create_PPRS_excel.py` — entry point, Excel export
- `real_models/` — JSON model data (EwE, species groups, diet matrices)

## graphify

This project has a knowledge graph at graphify-out/ with god nodes, community structure, and cross-file relationships.

Rules:
- For codebase questions, first run `graphify query "<question>"` when graphify-out/graph.json exists. Use `graphify path "<A>" "<B>"` for relationships and `graphify explain "<concept>"` for focused concepts. These return a scoped subgraph, usually much smaller than GRAPH_REPORT.md or raw grep output.
- If graphify-out/wiki/index.md exists, use it for broad navigation instead of raw source browsing.
- Read graphify-out/GRAPH_REPORT.md only for broad architecture review or when query/path/explain do not surface enough context.
- After modifying code, run `graphify update .` to keep the graph current (AST-only, no API cost).
