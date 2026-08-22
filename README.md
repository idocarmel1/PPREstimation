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
- `create_PPRS_excel.py` — batch exporter; runs every SPPR method over a directory of models and writes one multi-sheet Excel workbook per model, plus a reader that loads them back.
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

# Is that result trustworthy? One call grades the input data, both convergence
# conditions, and the PP budget. return_sppr=True reuses the same solve.
report, sppr_new, A, L = model.diagnose_sppr(return_sppr=True)
print("status:", report['status'])          # 'OK' | 'WARN' | 'FAIL'
print(report['warnings'])

# Ecosystem footprint: total PPR appropriated by the catch (1-row DataFrame).
ppr = model.get_PPR(sppr_new)
print("total PPR:", ppr.sum(axis=1).sum())

# Uncertainty propagation over transfer efficiency. Pick the solver with `kind`, give it
# its own parameters with `method_kwargs`, and drop draws whose TE sample pushed the
# recycling gain past divergence (which a negative-SPPR test alone can miss).
mean, samples, reject_frac, _, _, diag = model.monte_carlo_SPPR(
    n_samples=200, kind='new', TE_option='TE',
    method_kwargs={'fix_EE_0_cases': True},
    exclude_diverged=True, return_diagnostics=True)
print(f"rejected {reject_frac:.0%}: "
      f"{diag['n_rejected_diverged']} diverged, {diag['n_rejected_negative']} negative")
```

### Batch export

`create_PPRS_excel.py` runs all 19 SPPR method configurations over a directory of models
and writes **one workbook per model** into the output directory, alongside a plain-text
run report listing every warning, skipped method and Monte-Carlo rejection:

```bash
# defaults to real_models/EwE_jsons -> output/
python create_PPRS_excel.py [json_dir] [out_dir] [--resume] [--timeout SECONDS]
```

Each workbook holds eight sheets: `groups_df` (per-group parameters); `sppr_PP`,
`sppr_inner` and `sppr_all` — the same flat groups × methods table under three source
scopes (primary producers only; within-system, i.e. adding detritus; and every source,
adding import); `model_health` (`diagnose_sppr` per `TE_option`), `footprint` (PPR and
%NPP per method), `mc_diagnostics` (Monte-Carlo accept/reject breakdown) and `run_notes`
(the conventions needed to read the numbers correctly).

### Per-method time budget

Every method — and every `diagnose_sppr` row — runs in a worker process with a wall-clock
budget, **180 s by default**. The heavy methods (`SPPR_EwE`'s path enumeration, the
Monte-Carlo drivers) can run for hours on a large food web, and once they are down in
numpy/sympy there is no way to interrupt them from inside the interpreter — so a method
that overruns has its worker killed, is left NaN and is recorded as a *timeout*, and the
export moves straight on to the next method. One worker is shared by the whole directory,
so the budget costs a single process spawn, not one per model.

A timeout is not a failure of the method, and the workbook keeps the two apart: the
`status` column of the `run_notes` sheet says `ok` / `failed` / `timeout` / `crashed` per
method, and the run report gives the reason in full. Rerun with a larger budget before
reading anything into a timeout.

```bash
python create_PPRS_excel.py real_models/EwE_jsons output --timeout 600  # 10 minutes
python create_PPRS_excel.py real_models/EwE_jsons output --timeout none # no budget at all
```

```python
cpe.run_directory('real_models/EwE_jsons', 'output', method_timeout=600)
cpe.run_directory('real_models/EwE_jsons', 'output', method_timeout=None)  # in-process
```

From Python you can run a subset and read the result back:

```python
import create_PPRS_excel as cpe

# Skip the heavy path-enumeration method while iterating.
summary = cpe.run_directory('real_models/EwE_jsons', 'output',
                            method_keys=cpe.FAST_METHOD_KEYS)
print(summary['n_written'], 'workbooks;  report:', summary['report'])

tables = cpe.read_pprs_excel('output/435_435_Black_Sea_(1990).xlsx')
tables['sppr_all']['new_GE']    # total SPPR per group for one method
tables['sppr_PP']['new_GE']     # the primary-producer-only share
```

Throughout the workbook **NaN means "not available", never zero** — a method that raised,
one that ran out of time, and a value that cannot be attributed to the sheet's source
scope all stay empty.
`SPPR_1986` and `SPPR_1995` return a single un-attributed SPPR, so they are blank in
`sppr_PP` while carrying their total in `sppr_inner` and `sppr_all`.

## Setup

Requires Python 3.10+ and the scientific Python stack: `numpy`, `pandas`, `scipy`,
`sympy`, `igraph`, `tqdm`, and `openpyxl` (the Excel engine used by the batch exporter).
There is no packaged dependency file yet; install these manually, e.g.:

```bash
pip install numpy pandas scipy sympy igraph tqdm openpyxl
```

Tests are plain `pytest`:

```bash
python -m pytest tests/ -q
```

Note that `PPRCalculator` imports `tqdm.notebook`, so calling `SPPR_EwE` or
`monte_carlo_SPPR` with `silent=False` outside Jupyter raises unless `ipywidgets` is
installed. Pass `silent=True` in scripts — `create_PPRS_excel.py` already does.

## Documentation

- [`SPPR_Methods.md`](SPPR_Methods.md) — deep reference on every SPPR method and the underlying mathematics.
- [`USER_GUIDE.md`](USER_GUIDE.md) — practical guide to `ModelData` and `PPRCalculator`, with an end-to-end example.
