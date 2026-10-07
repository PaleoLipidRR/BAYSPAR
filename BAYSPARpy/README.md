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

## Compared with `brews/baysparpy`

S. B. Malevich's [baysparpy](https://github.com/brews/baysparpy) is the existing Python port, and a
careful one. `notebooks/compare_implementations.ipynb` runs both against the MATLAB reference on
Jess Tierney's two demo examples and shows the answers side by side. In short:

| | `brews/baysparpy` 0.0.3 | this port |
|---|---|---|
| **Agreement with MATLAB** | prior means exact, analogue cell set exact, medians inside sampling noise | the same, plus 18 tests that assert it against the reference executing |
| **Which posterior draws** | the **first** `nens` of 20,000 | MATLAB's `round(linspace(1, 20000, n))`, spread across the chain |
| **Percentiles** | `np.percentile(interpolation='nearest')` | MATLAB's `prctile` convention (Hazen), verified exactly against the reference's own ensemble |
| **Default ensemble** | 5000 | 1000, MATLAB's |
| **Speed** | ~1.9 s for a 193-point record at 1000 draws, 3.8 of 4 cores busy | ~0.013 s, one core — **~140×** |
| **Analogue ensemble shape** | `(N, analogues, draws)` | the same |
| **τ² pairing in analogue mode** | paired correctly (MATLAB does not) | the same, with `mode="matlab"` to reproduce the reference |
| **NumPy ≥ 2.0** | `Prediction.percentile()` raises — `interpolation=` was removed in NumPy 2.0 | works |
| **Provenance** | — | every translated line documented in `PORTING.md` under a stable ID, with a test per entry |

The measured effect of those choices on Jess's standard-mode demo: draw selection moves the median
0.08 °C on average, the percentile convention 0.004 °C, and reseeding the same implementation moves
it 0.15 °C. The differences are smaller than the sampling noise. In *analogue* mode the draw choice
does show as a systematic ~0.1 °C offset, because every analogue location reuses the same first
1000 draws rather than averaging the difference out.

The speed difference is structural, not tuning: the prior covariance is diagonal, so the posterior
is too and the update is elementwise. `brews/baysparpy` builds an N × N matrix and takes a `solve`
plus a Cholesky once per draw — and once per draw *per analogue location*. That is where both the
wall-clock and the all-cores behaviour come from.

```bash
uv sync                              # .venv with the package, dev tools and brews/baysparpy
uv run jupyter lab notebooks/compare_implementations.ipynb
```

or with pip:

```bash
pip install -e ".[dev,compare]"      # `compare` pulls brews/baysparpy (and cartopy)
jupyter lab notebooks/compare_implementations.ipynb
```

The notebook runs without brews/baysparpy, and reports that it is missing rather than failing.

## Reviewing it

`../docs/BAYSPARpy/PORTING.md` walks every line of the four MATLAB functions under a stable ID and
says what this code does instead. Ten entries change behaviour and sit in a register at the top;
`tests/test_port.py` has one test per ID and `tests/test_porting_doc.py` fails the build if an entry
loses its test or a difference escapes the register.

Not yet written: plotting, the demo notebooks, the CmdStan and Stan layers (SPEC §9–§10), and the
three test files that need a MATLAB session or `brews/baysparpy` installed. Until those golden tests
exist, this port is checked for internal consistency and against traced intermediates — not yet
against MATLAB itself.
