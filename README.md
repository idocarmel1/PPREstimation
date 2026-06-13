# FishEstimationAI

Estimates global **Primary Production Required (PPR)** from fish catches using multiple SPPR methods applied to Ecopath marine ecosystem models.

## Overview

PPR quantifies how much oceanic primary production (phytoplankton) is required to sustain observed fish catches, accounting for food-web structure and trophic efficiency losses. This project implements and compares several published SPPR estimation methods — from Pauly & Christensen (1995) to modern symbolic and Monte Carlo approaches — across hundreds of real Ecopath models.

## Project Structure

```
FishEstimationAI/
├── PPRCalculator.py          # Core computation class — the central hub
├── ModelData.py              # Data loading and species group management
├── PPR_methods.py            # SPPR method implementations
├── utils.py                  # Matrix utilities (cycle removal, etc.)
├── remove_cycles_fix.py      # Ulanowicz cycle-removal algorithm
├── create_PPRS_excel.py      # Entry point — runs all methods, exports to Excel
├── calc_2015/
│   └── SPPR_2015.py          # 2015 method module
├── real_models/              # JSON model data (EwE ecosystem models)
└── notebooks/                # Analysis and comparison notebooks
```

## SPPR Methods

| Method | Key | Description |
|--------|-----|-------------|
| SPPR 1986 | `SPPR_1986` | Pauly 1986 original method |
| SPPR 1995 mTL | `SPPR_1995` | Pauly & Christensen 1995 — mean trophic level, global TE |
| SPPR 1995 TL-fix | `SPPR_1995_TL_fix` | 1995 variant using TL2 |
| SPPR 2015 | `SPPR_2015` | Extended 2015 formulation |
| SPPR EwE | `SPPR_EwE` | Ecopath with Ecosim network-based estimate |
| SPPR New | `SPPR_new` | Revised full method with detritus handling |
| SPPR Symbolic | `SPPR_symbolic` | Symbolic/analytical derivation |
| Monte Carlo | `monte_carlo_SPPR` | Uncertainty sampling over gross efficiency (n=100) |

## Data

Models are stored as JSON files in `real_models/` following the Ecopath JSON format. Each model is identified by a numeric model number. `ModelData` supports two APIs:

- **New**: `ModelData(json_filepath="real_models/[227] Gulf of Mexico 1990.json")`
- **Legacy**: `ModelData(model_number=227)`

Key matrices extracted per model: diet composition (DC), detritus fate, trophic levels (TL), transfer efficiencies (TE), and the Z (mortality) matrix.

## Usage

Run all SPPR methods across all available models and export results to Excel:

```bash
python create_PPRS_excel.py
```

Results are saved per model as Excel files. Models with more than one detritus row or zero catch are skipped automatically.

### Programmatic use

```python
from PPRCalculator import PPRCalculator

model = PPRCalculator(227)            # load by model number
print(model.get_PPR())                # total PPR
print(model.SPPR_1995())             # 1995 method estimate
print(model.SPPR_new())              # revised method estimate
print(model.monte_carlo_SPPR(n_samples=100))  # uncertainty bounds
```

## Notebooks

| Notebook | Purpose |
|----------|---------|
| `compare_spprs_pipeline.ipynb` | Full pipeline comparison across methods |
| `comp_paper.ipynb` | Reproduce paper figures |
| `ratios.ipynb` | PPR/NPP ratio analysis |
| `LIM_solver.ipynb` | Linear Inverse Modeling (LIM) for mass balance |
| `models_sanity_checks.ipynb` | Model validation and quality checks |
| `GS_distribution.ipynb` | Gross efficiency distribution analysis |

## Key Concepts

- **PPR** — Primary Production Required: total phytoplankton production needed to support catch
- **SPPR** — Species-specific PPR
- **NPP** — Net Primary Production (ocean baseline)
- **TL** — Trophic Level
- **TE / GE** — Transfer Efficiency / Gross Efficiency between trophic levels
- **DC** — Diet Composition matrix
- **Ecopath** — Mass-balance ecosystem modeling framework (Christensen & Pauly)
- **LIM** — Linear Inverse Modeling, used to fill missing mass-balance values

## Dependencies

- `numpy`, `pandas` — numerical and data operations
- `sympy` — symbolic computation (SPPR_symbolic method)
- `scipy` — optimization (LIM solver via SLSQP)
- `tqdm` — progress bars
