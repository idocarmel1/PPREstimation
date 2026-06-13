# File Reorganization — Option A (Light Touch)

**Date:** 2026-06-13  
**Approach:** Move only clear misfits; keep all core `.py` modules at root.

## Goal

Reduce root-level clutter from 15 loose files to 5 core Python modules + README/CLAUDE.md, without touching any import paths in the core library.

## Moves

| Source | Destination | Notes |
|--------|-------------|-------|
| `test_integration.py` | `tests/test_integration.py` | Add `sys.path.insert(0, project_root)` after `os.chdir()` so standalone execution still resolves imports |
| `debug_processing.py` | `scripts/debug_processing.py` | Same sys.path fix; create `scripts/` dir |
| `compare_spprs_pipeline.ipynb` | `notebooks/compare_spprs_pipeline.ipynb` | Create `notebooks/` dir |
| `GS_distribution.ipynb` | `notebooks/GS_distribution.ipynb` | |
| `LIM_solver.ipynb` | `notebooks/LIM_solver.ipynb` | |
| `models_sanity_checks.ipynb` | `notebooks/models_sanity_checks.ipynb` | |
| `tests.ipynb` | `notebooks/tests.ipynb` | Exploratory test notebook, not a pytest file |
| `PPRs.xlsx` | `output/PPRs.xlsx` | Create `output/` dir; add to `.gitignore` |

## Files That Stay at Root

`PPRCalculator.py`, `ModelData.py`, `utils.py`, `remove_cycles_fix.py`, `create_PPRS_excel.py`, `README.md`, `CLAUDE.md`

## Import Safety

- Core modules (`PPRCalculator`, `ModelData`, `utils`, `remove_cycles_fix`) stay at root — zero import changes needed anywhere.
- `test_integration.py` and `debug_processing.py` already use hardcoded `os.chdir(project_root)`. After the move, add `sys.path.insert(0, os.getcwd())` immediately after `os.chdir()` so module imports resolve when run as standalone scripts. Pytest runs from project root, so pytest-executed tests already work.
- Notebooks moved to `notebooks/` need a working-directory cell at the top (`os.chdir(project_root)`) so relative paths to `real_models/` etc. still resolve. Most already have this.

## Verification

Run `pytest tests/` from project root after all moves. All existing tests must pass.

## What This Does NOT Touch

- `ToyModels/` — already a directory, left as-is
- `PPR_ratios/` — already a directory, left as-is  
- `real_models/` — data directory, left as-is
- `calc_2015/` — if present, left as-is
- Any imports inside core `.py` files
