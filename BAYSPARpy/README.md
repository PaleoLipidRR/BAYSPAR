# BAYSPARpy

A Python recode of the MATLAB BAYSPAR package in the repository above — a Bayesian, spatially-varying
calibration for the TEX<sub>86</sub> paleothermometer (Tierney & Tingley, 2014, 2015).

**Prototype.** It lives inside the MATLAB repository for now so its tests can read `ModelOutput/`
directly and be compared against the reference line by line. It will move to
`PaleoLipidRR/BAYSPARpy` once the parameter store is converted (SPEC §8) and the MATLAB golden files
exist. Nothing is published to PyPI pending agreement with J. Tierney and S. B. Malevich.

```bash
cd BAYSPARpy
pip install -e .              # or: PYTHONPATH=src python -m pytest
python -m pytest              # 65 passing, 2 skipped (Phase 0 items)
```

The import package is `baysparpy`, not `bayspar` — `brews/baysparpy` already imports as `bayspar`,
and both need to coexist while this one is validated against it.

## Using it

```python
import baysparpy as bp

# Standard mode: T from TEX86 at a known location (bayspar_tex.m)
out = bp.bayspar_tex(tex86, lon=-17.66, lat=9.17, prior_std=6.0, runname="subT", seed=0)
out.preds            # (N, 3) — the 5th, 50th and 95th percentiles, degrees C
out.prior_mean       # searched from the 1-degree climatology within 500 km
out.grid_loc         # the 20x20 degree cell the parameters came from

# Analogue mode: deep time, where today's location is not informative
an = bp.bayspar_tex_analog(tex86, prior_mean=30.0, prior_std=20.0,
                           search_tol=2 * tex86.std(ddof=1), runname="SST",
                           save_ensemble=True, seed=0)
an.analog_locs       # (n_analogues, 2) centroids, aligned with axis 1 of the ensemble
an.ensemble.shape    # (N, n_analogues, n_draws) — the shape MATLAB documents

# Forward: TEX86 from temperature (TEX_forward.m; note lat before lon, as in MATLAB)
fwd = bp.tex_forward(lat=9.17, lon=-17.66, t=[25.0, 26.0], runname="SST", seed=0)
fwd.metadata["clipped_fraction"]   # share of the ensemble pinned at 0 or 1
```

`predict_seatemp`, `predict_seatemp_analog` and `predict_tex` are the same three functions under
`brews/baysparpy`'s names and argument order, so a notebook ports by changing the import. The
default draw count is MATLAB's 1000, not that package's 5000.

## Where the data comes from

The stores read the MATLAB `ModelOutput/` directory, found by walking up from the package or named
by `BAYSPAR_MODELOUTPUT`. `params_analog.mat` is never read: it is bit-identical to the coretop rows
of `params_standard.mat` (PORTING.md STO-01), so analogue mode indexes the standard store.

## Reviewing it

`../docs/BAYSPARpy/PORTING.md` walks every line of the four MATLAB functions under a stable ID and
says what this code does instead. Ten entries change behaviour and sit in a register at the top;
`tests/test_port.py` has one test per ID and `tests/test_porting_doc.py` fails the build if an entry
loses its test or a difference escapes the register.

Not yet written: plotting, the demo notebooks, the CmdStan and Stan layers (SPEC §9–§10), and the
three test files that need a MATLAB session or `brews/baysparpy` installed. Until those golden tests
exist, this port is checked for internal consistency and against traced intermediates — not yet
against MATLAB itself.
