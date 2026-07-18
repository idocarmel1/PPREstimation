# FishEstimationAI

Fish biomass / **PPR** (Primary Production Required) estimation from Ecopath marine
ecosystem models. Given a mass-balanced food web, the tool computes **SPPR** (Specific
Primary Production Required) for every functional group and rolls it up into
ecosystem-level footprints — total PPR appropriated by the catch and %PPR relative to
net primary production.

## What it does

SPPR is the amount of primary production required to sustain one unit of a group's
production. FishEstimationAI estimates it with several methods, from classic
trophic-level approximations (Pauly & Christensen 1986; Christensen & Pauly 1995) to
modern flow-network solvers that trace production through the full diet matrix, handle
detritus recycling, and support uncertainty propagation. Weighting each group's SPPR by
its catch yields the ecosystem PPR footprint.

See [`SPPR_Methods.md`](SPPR_Methods.md) for the full method catalog and math.

## Project structure

- `PPRCalculator.py` — central computation class; SPPR methods, PPR/NPP aggregation, detritus solving.
- `ModelData.py` — model loading, species-group management, diet/detritus-fate matrices.
- `utils.py` — matrix utilities and Ulanowicz cycle removal (`remove_cycles`).
- `create_PPRS_excel.py` — command-line entry point; runs the methods and exports results to Excel.
- `real_models/` — JSON model data (EwE parameters, species groups, diet matrices).
- `graphify-out/` — pre-built knowledge graph of the codebase.

## Usage

```python
from PPRCalculator import PPRCalculator

# Load and mass-balance a model by its number (or JSON filepath).
model = PPRCalculator(227)

print("groups:", model.n_groups)
print("balanced?", model.is_model_balanced()[0])

# Classic trophic-chain method (per-group SPPR).
sppr_1995 = model.SPPR_1995(global_TE=0.1)

# Full flow-network solver; returns (SPPR, A, L).
sppr_new, A, L = model.SPPR_new(TE_option='GE')

# Ecosystem footprint: total PPR appropriated by the catch (1-row DataFrame).
ppr = model.get_PPR(sppr_new)
print("total PPR:", ppr.sum(axis=1).sum())
```

To run the batch pipeline and produce an Excel report:

```bash
python create_PPRS_excel.py
```

## Setup

Requires Python 3.10+ and the scientific Python stack: `numpy`, `pandas`, `scipy`,
`sympy`, `igraph`, and `tqdm`. There is no packaged dependency file yet; install these
manually, e.g.:

```bash
pip install numpy pandas scipy sympy igraph tqdm
```

## Documentation

- [`SPPR_Methods.md`](SPPR_Methods.md) — deep reference on every SPPR method and the underlying mathematics.
- [`USER_GUIDE.md`](USER_GUIDE.md) — practical guide to `ModelData` and `PPRCalculator`, with an end-to-end example.
