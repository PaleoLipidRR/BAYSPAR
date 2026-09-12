# BAYSPARpy — specification for a Python recode of BAYSPAR

**Status:** draft for review (Ronnakrit Rattanasriampaipong) — not yet agreed with J. Tierney or S. B. Malevich.
**Target repo:** `PaleoLipidRR/BAYSPARpy` (new repository; this document lives in the MATLAB repo until that repo exists).
**Reference implementation:** the MATLAB code in this repository (Tierney & Tingley, 2014, 2015).
**Sibling project:** `PaleoLipidRR/TEXAS` (`texas-psm`), whose conventions this port deliberately mirrors.

Every factual claim about the MATLAB package in §3–§5 was verified against the files in
`ModelOutput/` by the scripts in `audit/`; re-run them with `python audit/audit_modeloutput.py`
and `python audit/bench_closed_form.py`.

**The line-by-line translation record is [`PORTING.md`](PORTING.md)** — 45 numbered entries covering
every line of the four MATLAB functions, with a register of the nine places the Python does not
behave like MATLAB. That is the document to review the recode against; this one states what is being
built and why. `audit/reference_trace.txt` holds the deterministic values to compare against a
MATLAB session.

---

## 1. Goals and non-goals

### Goals

1. **A faithful Python recode of the MATLAB BAYSPAR package**, keeping the original project structure
   and function names recognisable (§6, §7). A user who knows `bayspar_tex.m` should find
   `bayspar_tex()` doing the same thing with the same arguments.
2. **Numerical equivalence with the MATLAB reference**, demonstrated by golden-file tests rather than
   asserted (§11).
3. **A reviewable recode.** Every translated line is accounted for in `PORTING.md` under a stable ID,
   every intentional difference is in its register with a named decision-maker, and a CI check fails
   if an ID there has no corresponding test. A reviewer can check the port without reading the Python,
   and can object to a numbered decision rather than to the port as a whole.
4. **Predictable, modest resource use.** The closed-form prediction must not saturate every core on a
   laptop. This is the direct motivation for the rewrite (§5).
5. **A CmdStan setup and solving path equivalent to TEXAS's** (§9, §10), so that:
   - someone who only wants BAYSPAR never has to install `texas-psm`; and
   - once CmdStan is installed for either package, the *same* installation serves both.
6. **The ability to refit the calibration on new coretop data** (§10.3), which the MATLAB release does
   not provide — it ships fitted posteriors only.
7. **A clean migration path into `texas-psm`** once the TEXAS paper is out (§12).

### Non-goals (v1)

- Reproducing MATLAB's random number *stream*. Draws are stochastic; equivalence is distributional
  (§11.2).
- New science. The default calibration stays TT14/TT15 with the shipped posteriors. Refitting is a
  capability, not a change of the shipped defaults.
- Replacing `brews/baysparpy` on PyPI. That is a conversation with Jess and Brewster, not a technical
  decision (§2).

---

## 2. Naming and distribution

`baysparpy` on PyPI is taken by `brews/baysparpy` 0.0.3, whose import package is `bayspar`.
`texas-psm`'s `dev` extra installs it, so both packages will be importable in the same environment
during validation.

| Item | Value | Rationale |
|---|---|---|
| GitHub repo | `PaleoLipidRR/BAYSPARpy` | Distinct from `brews/baysparpy`; matches the import name. |
| Import package | `baysparpy` | **Must not be `bayspar`** — that is Brewster's import name, and `BAYSPAR` would collide with it on case-insensitive filesystems (macOS, Windows). `import baysparpy` coexists with `import bayspar` in one environment, which the validation suite needs. |
| Distribution name | `bayspar` (reserved, **not published**) | Free on PyPI as of this writing. Nothing is uploaded until Jess and Brewster agree on how this relates to `brews/baysparpy`. |
| Version | `0.1.0.dev0` | Prototype. |

**Decision deferred, deliberately:** if Brewster agrees to hand over the `baysparpy` name, the
distribution name changes to `baysparpy` and nothing else moves — the import name is already
`baysparpy`, so both outcomes converge. Until then the package is installed from a git checkout
(`pip install -e .`), and the README says so in the first paragraph.

---

## 3. What the MATLAB package actually contains (audited)

### 3.1 File inventory

| File | Contents (verified) |
|---|---|
| `bayspar_tex.m` | Standard-mode prediction T ← TEX₈₆. |
| `bayspar_tex_analog.m` | Analogue-mode prediction. |
| `TEX_forward.m` | Forward model TEX₈₆ ← T, standard or analogue. |
| `EarthChordDistances_2.m` | Chordal distance, R = 6378.137 km. |
| `Demo_StandardPrediction.m`, `Demo_AnalogPrediction.m` | Worked examples. |
| `ModelOutput/Output_SpatAg_{SST,subT}/params_standard.mat` | `alpha_samples_comp` (162 × 20000), `beta_samples_comp` (162 × 20000), `tau2_samples` (1 × 20000), `Locs_Comp` (162 × 2, int16). |
| `ModelOutput/Output_SpatAg_{SST,subT}/params_analog.mat` | `alpha_samples` (80 × 20000), `beta_samples` (80 × 20000), `tau2_samples` (1 × 20000). |
| `ModelOutput/obs{SST,subT}.mat` | `locs_st_obs` (37686 × 2 SST / 37105 × 2 subT), `st_obs_ave_vec` — 1° WOA climatology used for the prior mean. |
| `ModelOutput/Data_Input_SpatAg_{SST,subT}.mat` | The calibration inputs (§3.3). |
| `ModelOutput/tex_testdata.mat`, `wilsonlake.mat` | Demo series. |

### 3.2 Three findings that shape the design

1. **The ensemble is 20,000 draws, not 15,000.** Both `ReadMe.md` and the function headers say the
   number of draws "cannot exceed 15000". Every shipped `tau2_samples` is 1 × 20000. The port
   validates `n_draws` against the actual array length and the docs are corrected.

2. **`params_analog.mat` is redundant.** For both SST and subT, `alpha_samples`, `beta_samples` and
   `tau2_samples` in `params_analog.mat` are **bit-identical** to the rows of `params_standard.mat`
   selected by `Data_Input.Locs`, in that order (`np.array_equal` → `True`, all six arrays). The
   analogue files are a 80-row subset of the standard files, not a separate fit.
   *Consequence:* the port stores one parameter set per run (162 cells) plus an 80-element index of
   the data-bearing cells. That removes 97 MB of the 145 MB store and one whole class of "which file
   was this fitted from?" confusion.

3. **`Data_Input` is the coretop calibration dataset**, in a form that supports refitting (§3.3).

### 3.3 `Data_Input` structure (needed for the refit, §10.3)

Fields, verified for SST (subT identical in shape, 906 rather than 903 stacked observations):

| Field | Shape | Meaning |
|---|---|---|
| `Locs` | 80 × 2 int16 | Centroids `[lon, lat]` of the 20° × 20° cells that contain coretop data. All 80 appear in `Locs_Comp`. |
| `ObsAll` | 35 × 80 cell | TEX₈₆ observations, indexed `[site slot, cell]`. |
| `ObsLocs` | 35 × 80 cell | The 1° location of each site slot — exactly one `[lon, lat]` per non-empty entry (570 non-empty). |
| `Target` | 1 × 35 struct | Per site slot: `Field` (80 × 1) target temperature at that site's 1° location, `Err_Stds` (80 × 1) its uncertainty, plus stacking indices. |
| `Inds_Stack`, `Obs_Stack`, `Target_Stack`, `Target_Err_Stds_Stack` | 903 × 1 | The same data flattened; `Inds_Stack` ∈ 1..80 is the cell index. |

So the second index of `ObsAll`/`ObsLocs` is the 20° cell and the first is the *site slot within that
cell* — slot 0 is populated for all 80 cells, slot 34 for only one. Observation counts per slot are
`[135, 104, 95, 86, 56, …, 1]`, summing to 903. Each site has its own target temperature *with a
stated error standard deviation*, which is why the calibration is an errors-in-variables regression.

### 3.4 Grid geometry

`Locs_Comp` is the complete 18 × 9 tiling of the globe: longitudes at odd multiples of 10° from −170
to 170, latitudes at even multiples of 20° from −80 to 80, each cell ±10°. Cell lookup in the MATLAB
code is `abs(Locs_Comp(:,1)-lon) <= 10 & abs(Locs_Comp(:,2)-lat) <= 10`, i.e. **closed on both
edges** — see §4 item 5.

---

## 4. MATLAB → Python behavioural gotchas

Each of these is a place where a naive transcription silently changes results. The port handles each
one explicitly and tests it.

1. **Percentiles.** `prctile(X, [5 50 95], 2)` uses Hazen plotting positions `(i − 0.5)/n` with
   linear interpolation and clamping at the extremes. NumPy's default (`method='linear'`, type 7)
   does **not** match. Use `np.percentile(..., method='hazen')`. A golden-file test pins this against
   MATLAB output on a fixed ensemble.
   (The `sort(Preds,2)` inside the MATLAB call is a no-op — `prctile` sorts anyway.)

2. **Thinning.** `ind_s = round(linspace(1, 20000, Nsamps))`. MATLAB's `round` is half-away-from-zero;
   `np.round` is half-to-even. Use `np.floor(x + 0.5)` and subtract 1 for 0-based indexing. The
   existing Colab notebook in `BAYSPAR/BAYSPAR_recode_colab.ipynb` uses `np.arange(0, N, N//Nsamps)`
   instead and flags it in a comment as a deviation — that deviation does not carry over.
   `brews/baysparpy` takes the **first** `nens` draws (`alpha_samples_comp[jj]` for `jj in
   range(nens)`), which is a third behaviour again; MATLAB deliberately spreads the thinning across
   the full chain.

3. **Reshape order.** MATLAB `reshape` is column-major. Any port of the analogue-mode flattening must
   use `order='F'` or index arithmetic that reproduces it.

4. **`tau2` is mis-paired with `alpha`/`beta` in analogue mode.** In `bayspar_tex_analog.m` (and the
   analogue branch of `TEX_forward.m`):

   ```matlab
   tau2_samples = repmat(tau2_samples, 1, size(alpha_samples,1));   % [d1..dM, d1..dM, ...]
   alpha_samples = reshape(alpha_samples, 1, n_an*M);               % column-major: (loc1,d1),(loc2,d1),...
   ```

   At flat position *k*, `alpha`/`beta` come from draw `floor((k−1)/n_an)+1` while `tau2` comes from
   draw `mod(k−1, M)+1`. For more than one analogue location these are different draws, so the joint
   posterior correlation between (α, β) and τ² is broken. `alpha` and `beta` stay correctly paired
   with each other.
   *Expected impact:* small — the draws are exchangeable within the chain, so the pooled ensemble
   contains the same α, β and τ² values, merely re-paired. *Treatment:* `mode="matlab"` reproduces
   the mis-pairing exactly; the default path pairs by draw index, and a test reports the difference
   in the 5/50/95 percentiles on the Wilson Lake demo so the magnitude is on record rather than
   assumed. (TEXAS's inverse-T model states the same requirement — "ALL PARAMETERS MUST USE THE SAME
   DRAW INDEX m" — for the same reason.)

5. **Sites on a cell boundary match more than one cell.** The lookup is inclusive at ±10°, so a site
   at lon = 0°, or lat = 10°, matches two cells (a corner matches four). In MATLAB the resulting
   2 × M parameter block makes `dats - repmat(...)` a dimension error when N > 1, and *silently
   returns one row per matched cell* when N = 1. `brews/baysparpy` has the same defect via
   `.squeeze()`. The port resolves ties deterministically — nearest centroid by chordal distance,
   with a `UserWarning` naming the cells — and `strict_grid=True` raises instead.

6. **Prior-mean search.** `find(vals < max_dist, 1, 'last')` on the sorted distances is a count of
   observations strictly within 500 km; if none qualify the empty comparison is false and the single
   closest observation is used. Reproduce with `np.searchsorted(sorted_d, 500.0, side='left')` and
   the same fallback, with a stable sort so ties order identically.

7. **Random numbers.** MATLAB `randn`/`randsample` streams cannot be reproduced in NumPy. All public
   functions take `seed`/`rng`; equivalence with MATLAB is tested distributionally (§11.2).

8. **`TEX_forward` with multiple locations.** The vectorised branch requires `len(t) == len(lat) ==
   len(lon)` — one temperature per site, not a time series at each of several sites. The Python
   signature makes this explicit and validates it.

9. **TEX₈₆ clipping.** `TEX_forward` clamps simulated values into [0, 1] *after* sampling. Keep it,
   and record the fraction clipped in the result metadata — for cold subT draws it is not negligible.

10. **`Locs_Comp` is int16.** Cast to float before differencing against float coordinates.

---

## 5. Why a recode rather than `brews/baysparpy` (measured)

`brews/baysparpy` 0.0.3 is a careful port, but its solver is written for a *general* prior covariance
matrix: `bayspar/utils.py::target_timeseries_pred` builds an N × N inverse posterior covariance,
calls `np.linalg.solve` against an N × N identity, and takes a Cholesky factor — **once per posterior
draw**, and in analogue mode once per draw *per analogue location*.

As actually called, the prior covariance is always `np.eye(nd) * prior_std**-2` — diagonal. The
posterior covariance is therefore diagonal too, and the entire N × N linear algebra is avoidable.
The MATLAB original never does it: `bayspar_tex.m` solves elementwise.

Measured on the container this spec was written in (4 cores, Python 3.11, NumPy 2.4.6, OpenBLAS), on
the `lopes_santos2010` subT demo series, `n_draws = 1000`:

| Record length N | Dense per-draw solve | Vectorised closed form | Speed-up |
|---|---|---|---|
| 50 | 0.18 s | ~1 ms | ×170 |
| 200 | 1.6 – 2.2 s | 7 – 24 ms | ×70 – ×320 |
| 500 | 13.6 – 16.5 s | 31 – 91 ms | ×150 – ×530 |

(Three repeats. The dense column is stable; the closed-form column is milliseconds on a shared
container and therefore noisy, so read the ratio as "two to three orders of magnitude", not as a
precise figure.)

And the CPU load, on the same dense loop (N = 300, 400 draws):

| | wall | CPU time | cores busy |
|---|---|---|---|
| default BLAS | 1.9 – 2.0 s | 7.4 – 7.9 s | **4.0 of 4** |
| `OMP_NUM_THREADS=1` | 3.6 s | 3.6 s | 1.0 |

That is the "it uses all my CPUs" behaviour: LAPACK inside the per-draw loop hands the work to a
threaded BLAS. Note the default is faster in wall-clock — the threads are not the bug, the O(N³) work
per draw is. At the package default of `nens=5000`, and in analogue mode across ~20 analogue cells, a
500-point record costs on the order of `15 s × 5 × 20 ≈ 25 minutes` of fully-saturated CPU, for
arithmetic that takes a few seconds.

**Design rule:** the closed-form path is pure elementwise NumPy — no BLAS call, no thread pool, one
core. Any future parallelism is opt-in and explicit (`n_jobs`), never a side effect of a library.

Two further deviations from MATLAB worth recording, both verified in the 0.0.3 source:
- it uses the first `nens` draws rather than MATLAB's spread thinning (§4 item 2);
- it ships only `*_comp` parameter files, and analogue mode reads α, β from those at the analogue
  cells. Because `params_analog.mat` is bit-identical to that subset (§3.2 finding 2), this happens
  to be **numerically correct** — worth stating explicitly, since it looks like a substantive
  deviation and is not.

---

## 6. Architecture and package layout

The MATLAB structure is preserved where it can be: one module per `.m` file, same names, plus a
modern API layered on top.

```
BAYSPARpy/
├── pyproject.toml
├── README.md
├── CHANGELOG.md
├── CITATION.cff
├── CLAUDE.md                        # working conventions, mirroring TEXAS's
├── LICENSE                          # inherit the MATLAB repo's licence
├── environment.yml                  # conda path, ships CmdStan (see §9)
├── src/baysparpy/
│   ├── __init__.py                  # public API re-exports
│   ├── constants.py                 # GRID_HALF_SPACE=10, EARTH_RADIUS_KM=6378.137,
│   │                                # MAX_DIST_KM=500, MIN_NUM=1, RECOMMENDED_CMDSTAN_VERSION
│   ├── bayspar_tex.py               # <- bayspar_tex.m
│   ├── bayspar_tex_analog.py        # <- bayspar_tex_analog.m
│   ├── tex_forward.py               # <- TEX_forward.m
│   ├── distance.py                  # <- EarthChordDistances_2.m
│   ├── predict.py                   # modern API (§7.2), wraps the three above
│   ├── results.py                   # Prediction container: ensemble, percentiles, metadata
│   ├── modelparams/                 # <- ModelOutput/Output_SpatAg_*/
│   │   ├── core.py                  # Draws: lazy load, thinning, cell lookup
│   │   └── registry.py              # store locations + SHA-256 checksums (§8)
│   ├── observations/                # <- ModelOutput/obs*.mat, Data_Input_*.mat
│   │   ├── seatemp.py               # prior-mean search
│   │   └── coretops.py              # Data_Input -> tidy frame, analogue search
│   ├── stan/
│   │   ├── compiler.py              # model compilation + cache
│   │   ├── invT.py                  # §10.2 inverse-T sampling
│   │   ├── calibrate.py             # §10.3 coretop refit
│   │   └── io.py                    # posterior <-> NetCDF
│   ├── stan_models/
│   │   ├── invT_bayspar_marginal.stan
│   │   ├── invT_bayspar_marginal_ar1.stan
│   │   └── bayspar_spatial_calibration.stan
│   ├── _cmdstan/                    # vendored from TEXAS, parity-tested (§9)
│   │   ├── paths.py  install.py  doctor.py  system_info.py
│   ├── data/                        # bundled 1000-draw default + example series
│   └── plot.py                      # predictplot, analogmap, densityplot
├── notebooks/
│   ├── demo_standard_prediction.ipynb   # <- Demo_StandardPrediction.m
│   ├── demo_analog_prediction.ipynb     # <- Demo_AnalogPrediction.m
│   └── demo_refit_coretops.ipynb        # new (§10.3)
├── tests/
├── docs/
└── tools/
    └── convert_modeloutput.py       # .mat -> NetCDF, one-time, reproducible
```

### Port map

| MATLAB | Python | Notes |
|---|---|---|
| `bayspar_tex(dats, lon, lat, prior_std, runname, ...)` | `baysparpy.bayspar_tex(...)` | Same argument order and defaults. |
| `bayspar_tex_analog(dats, prior_mean, prior_std, search_tol, runname, ...)` | `baysparpy.bayspar_tex_analog(...)` | Same. |
| `TEX_forward(lat, lon, t, ...)` | `baysparpy.tex_forward(...)` | Same (note lat before lon, as in MATLAB). |
| `EarthChordDistances_2(a, b)` | `baysparpy.distance.earth_chord_distances(a, b)` | Vectorised, no N·M intermediate. |
| `Output_Struct` | `Prediction` dataclass | `.preds`, `.site_loc`, `.grid_loc`, `.prior_mean`, `.prior_std`, `.ensemble`, `.metadata`; `.to_dataframe()`, `.to_xarray()`. |
| `ModelOutput/` | `baysparpy.modelparams` + cache dir | §8. |
| `Demo_*.m` | `notebooks/demo_*.ipynb` | Same figures. |

---

## 7. Public API

### 7.1 MATLAB-faithful layer

```python
bayspar_tex(dats, lon, lat, prior_std, runname="SST", *,
            n_draws=1000, save_ensemble=False, mode="modern",
            seed=None, strict_grid=False) -> Prediction

bayspar_tex_analog(dats, prior_mean, prior_std, search_tol, runname="SST", *,
                   n_draws=1000, save_ensemble=False, mode="modern",
                   seed=None) -> Prediction

tex_forward(lat, lon, t, runname="SST", type="standard", search_tol=None, *,
            n_draws=1000, seed=None) -> ForwardPrediction
```

`mode="matlab"` reproduces the reference exactly, including the τ² mis-pairing of §4 item 4 and the
legacy percentile convention; `mode="modern"` (default) pairs draws correctly. Both are tested; the
difference between them is reported, not hidden.

### 7.2 Modern layer

Deliberately aligned with TEXAS's `predict.py` so that calls port with a rename:

```python
predict_seatemp(tex, lat, lon, prior_std, temptype="sst", prior_mean=None, *,
                n_draws=1000, seed=None, backend="closed_form") -> Prediction
predict_seatemp_analog(tex, prior_std, search_tol, temptype="sst", prior_mean=None, ...) -> Prediction
predict_tex(seatemp, lat, lon, temptype="sst", ...) -> ForwardPrediction
```

`backend="closed_form"` (default, §5) or `backend="stan"` (§10.2). The two are the same model; the
Stan backend differs in *how* calibration uncertainty is integrated (§10.1) and is the entry point
for priors the closed form cannot express.

Name compatibility with `brews/baysparpy` is intentional: the modern names and argument names match
his, so a notebook switches by changing `from bayspar import …` to `from baysparpy import …`. Where
his defaults differ from MATLAB (`nens=5000` vs 1000) the MATLAB default wins and the difference is
documented.

### 7.3 Environment / toolchain

```python
baysparpy.doctor()                    # environment report, same layout as TEXAS.doctor()
baysparpy.install_cmdstan()           # opt-in CmdStan install
baysparpy.download_params(...)        # fetch the full 20,000-draw store (§8)
baysparpy.set_cache_dir(path)
```
Console scripts: `bayspar-doctor`, `bayspar-install-cmdstan`.

---

## 8. Data packaging

The parameter store is 145 MB of `.mat`. After dropping the redundant analogue files (§3.2) it is
97 MB as float64, 44 MB as float32 (measured, compressed `.npz`; the float32 round-trip error on α is
4.3 × 10⁻⁸ relative).

**Decision: no Zenodo deposit. Everything lives in the repository, with a small tier in the wheel.**

| Tier | Contents | Size | Where |
|---|---|---|---|
| Bundled | Exactly the MATLAB 1000-draw thinning of α, β, τ² for SST and subT, float64, plus `Locs_Comp`, the 80-cell index, the coretop `Data_Input` tables, and the WOA prior-mean vectors | ~6 MB | in the wheel |
| Full store | All 20,000 draws, float64, NetCDF | ~97 MB | in the repository, via Git LFS (`git lfs pull`), as TEXAS tracks its data |

The bundled tier *is* the MATLAB default ensemble, so `pip install` plus the demo works offline and
nothing about the standard result depends on fetching anything. Asking for `n_draws > 1000`, or for a
custom thinning, needs the full store and raises a message naming the `git lfs pull` (or the
`baysparpy.download_params()` call, if tier 3 below ever exists).

*Contingency, not built now:* if a PyPI release later makes a 97 MB checkout unacceptable and a
network tier is needed, host it under a PaleoLipidRR-owned record the same way TEXAS does, with the
same registry shape, SHA-256 checksums and `download_*()` API as `TEXAS/utils/download.py`. No
third-party deposit is created in the meantime.

Storage format: NetCDF via xarray, dimensions `(cell, draw)`, coordinates `lon`, `lat`, attributes
recording provenance (source `.mat` filename, its SHA-256, conversion date, `baysparpy` version).
`tools/convert_modeloutput.py` performs the conversion and is committed, so the store is reproducible
from the MATLAB files rather than a binary blob of unknown origin.

---

## 9. CmdStan integration (shared with TEXAS)

**Requirement (user):** someone using only BAYSPARpy must be able to set up CmdStan from BAYSPARpy
alone; and a CmdStan set up through either package must work for both.

**Decision: vendor TEXAS's toolchain module into `baysparpy/_cmdstan/`, with a parity test.**
Neither package depends on the other. Because both resolve CmdStan the same way — same env var, same
search order, same `set_cmdstan_path` call — one installation serves both automatically.

Ported from `TEXAS/utils/`:

- `find_cmdstan(min_version="2.23.0")` — search order preserved exactly: `CMDSTAN` env var →
  `CONDA_PREFIX/bin/cmdstan` → `sys.prefix/bin/cmdstan` → highest `cmdstan-*` under `/opt/cmdstan`,
  `~/.cmdstan`, `/usr/local/cmdstan` → whatever cmdstanpy already holds. Always calls
  `set_cmdstan_path()`.
- `install_cmdstan(version=RECOMMENDED_CMDSTAN_VERSION, ...)` — opt-in, never runs on import; no-ops
  when a working install resolves; steps aside for conda-managed installs; auto-enables `overwrite`
  for a half-built directory; prints `doctor()` afterwards.
- `doctor()` — CmdStan presence and version, compiler, toolchain conflicts, package versions.
- `RECOMMENDED_CMDSTAN_VERSION` pinned to the same value as TEXAS (2.36.0 at the time of writing).

**Parity test** (`tests/test_cmdstan_parity.py`): skipped when `texas-psm` is absent; when present it
asserts that `baysparpy._cmdstan.find_cmdstan()` and `TEXAS.utils.paths.find_cmdstan()` resolve to
the same path, that the recommended versions match, and that the documented search order in both
docstrings is identical. Drift becomes a failing test rather than a support ticket.

**Interop test** (`tests/test_cmdstan_interop.py`, marked `slow`): install CmdStan via
`baysparpy.install_cmdstan()`, then compile a TEXAS model with it, and vice versa.

Environment paths, in the TEXAS order of preference: conda-lock/`environment.yml` (bundles CmdStan),
pip/uv + `bayspar-install-cmdstan`, and a Docker image for the fully reproducible case.

---

## 10. Stan models

### 10.1 What Stan buys, stated honestly

The shipped BAYSPAR posteriors are already fitted, and the standard-mode inversion is conjugate, so
**Stan is not needed to reproduce MATLAB**. It is needed for three things:

1. **Correct marginalisation.** MATLAB draws, for each posterior sample *m*, one value from the
   conditional posterior `p(T | y, θ_m)` and pools with equal weight. That is
   `∫ p(T | y, θ) p(θ) dθ`. The marginal posterior is
   `p(T | y) ∝ p(T) · (1/M) Σ_m N(y | α_m + β_m T, τ²_m)`, which is the same mixture weighted by each
   draw's evidence `Z_m = ∫ p(T) N(y | α_m + β_m T, τ²_m) dT`. The two agree when the evidence is
   flat across draws and differ otherwise. The difference is quantified on the demo series and
   reported, not assumed small.
2. **Priors the closed form cannot express** — AR(1) or random-walk structure down a core, a
   truncated prior (no reconstruction below the freezing point), per-sample prior means, TEX₈₆
   measurement error.
3. **Refitting** (§10.3).

The closed form stays the default. Stan is opt-in.

### 10.2 Inverse temperature: `invT_bayspar_marginal.stan`

Directly analogous to `TEXAS/stan_models/invT_gen_logi_fixed_univ_marginal_unconstrained.stan`, with
BAYSPAR's linear forward model in place of the generalised logistic:

```stan
data {
  int<lower=1> N;                 // samples in the record
  vector[N] tex;                  // TEX86 observations
  vector[N] prior_mu_t;           // per-sample prior mean (degC)
  real<lower=0> prior_sigma_t;    // prior SD (degC)
  int<lower=1> M;                 // calibration draws
  vector[M] alpha;                // draw m of alpha  (paired by m)
  vector[M] beta;                 // draw m of beta
  vector<lower=0>[M] tau;         // draw m of sqrt(tau2)
  int<lower=1> grainsize;         // reduce_sum chunking; grainsize=N disables
}
parameters { vector[N] t_est; }
model {
  target += reduce_sum(ll_chunk, linspaced_int_array(N,1,N), grainsize,
                       tex, t_est, prior_mu_t, prior_sigma_t, alpha, beta, tau);
}
```

with `ll_chunk` adding `normal_lpdf(t_seg | mu_seg, prior_sigma_t)` and, per sample,
`log_sum_exp_m normal_lpdf(tex_i | alpha[m] + beta[m]*t_i, tau[m]) − log(M)`. The prior lives inside
the chunk function for the same reason as in TEXAS: `reduce_sum` partitions over *n*.

- **Analogue mode** needs no separate model: the pooled (α, β, τ) draws across analogue cells are
  passed as the M vectors, correctly paired by draw index. The mixture over locations is exactly what
  marginalising over that pooled set does.
- **`invT_bayspar_marginal_ar1.stan`** adds `t_est[n] ~ normal(t_est[n-1], sigma_rw)` for
  age-ordered records, with the age vector supplied so the innovation variance scales with Δt.
- **Threading.** `threads_per_chain` defaults to 1 and `grainsize` to N — the parallelism is opt-in,
  in keeping with §5. `STAN_THREADS` compilation matches TEXAS so the model cache is compatible.

Acceptance: on a record where the closed form and the marginal posterior should agree (flat evidence
across draws), the two must agree in 5/50/95 to within Monte-Carlo error; where they disagree, the
disagreement is characterised in `docs/marginalisation.md`.

### 10.3 Refitting on new coretops: `bayspar_spatial_calibration.stan`

**This capability does not exist in the MATLAB release** — the repository ships fitted posteriors but
not the sampler that produced them. The model must therefore be reconstructed from Tierney & Tingley
(2014) and validated against the shipped posteriors before it is trusted. That reconstruction is the
largest single work item in this spec and is scheduled last (§13, Phase 4).

Model, from §3.3 and TT14:

```
for cell c, site j in c, observation k:
    tex[c,j,k]  ~ Normal(alpha[c] + beta[c] * T_true[c,j], tau)
    T_obs[c,j]  ~ Normal(T_true[c,j], T_err_sd[c,j])        // errors-in-variables
    alpha ~ MVN(mu_alpha * 1, Sigma_alpha(sigma_alpha, phi))  // GP over cell centroids
    beta  ~ MVN(mu_beta  * 1, Sigma_beta(sigma_beta,  phi))
    Sigma(d) = sigma^2 * exp(-d / phi),  d = chordal distance between centroids
```

- **Seven of the 903 SST coretops carry a target error standard deviation of exactly zero**
  (none of the subT ones do; found while testing the port, pinned by `test_STO_04`). The
  errors-in-variables term divides by that, so the model must floor it, drop those observations, or
  treat their targets as known exactly — a decision to make explicitly rather than at a crash.
- Priors on `mu_alpha, mu_beta, sigma_alpha, sigma_beta, phi, tau` to be taken from TT14 §3 —
  **to be confirmed against the paper before implementation**; the kernel family (exponential vs
  Matérn) and whether α and β share one range parameter are the two specific points to check.
- `generated quantities` draws α, β at all 162 `Locs_Comp` centroids conditional on the 80 fitted
  cells — this is precisely the off-line interpolation step the MATLAB ReadMe describes
  ("performed off-line and results simply looked up"), so the refit pipeline emits a drop-in
  replacement for `params_standard`, and the 80-cell "analogue" set falls out as a subset (§3.2).
- 80 cells → dense 80 × 80 covariance; Cholesky per iteration is cheap. Non-centred parameterisation
  for α and β.

**Acceptance criterion.** Refit on the *shipped* `Data_Input` must reproduce the shipped posteriors:
posterior means within Monte-Carlo error at every data-bearing cell (target |Δα| < 0.01,
|Δβ| < 5 × 10⁻⁴, given α ∈ [−0.33, 0.68] and β ∈ [−0.006, 0.039]), posterior SD ratios in
[0.9, 1.1], and R̂ < 1.01. If the reconstruction does not meet this, the difference is documented and
the refit ships as clearly-labelled "reimplemented, not identical to TT14" rather than being quietly
presented as the same model.

API:

```python
baysparpy.stan.calibrate(coretops, target_field, *, runname="SST",
                         chains=4, iter_sampling=1000, threads_per_chain=1,
                         seed=None) -> CalibrationPosterior
CalibrationPosterior.to_netcdf(path)      # usable directly as a parameter store
```

`coretops` is a tidy DataFrame (`lon`, `lat`, `tex86`, `target_t`, `target_t_sd`); a loader turns
the shipped `Data_Input` into one, so "refit the original" and "refit with my new coretops" are the
same call with a different frame.

---

## 11. Validation and test plan

### 11.1 Golden files (deterministic)

Generated once from MATLAB and committed under `tests/golden/`:

| Quantity | Tolerance |
|---|---|
| `earth_chord_distances` on a fixed lon/lat grid | 1e-10 relative |
| Prior mean for each demo site, both runs | exact (float64 equality) |
| Grid-cell index for a list of sites incl. boundary cases | exact |
| Thinning index vector for `n_draws` ∈ {10, 999, 1000, 1001, 20000} | exact |
| Percentiles of a fixed ensemble matrix (`prctile` parity) | 1e-12 |
| Analogue cell selection for the Wilson Lake demo at several tolerances | exact set equality |
| Per-draw posterior mean and SD (pre-noise) for a fixed record | 1e-10 relative |

Deterministic quantities are the bulk of the port; testing them exactly is what makes the stochastic
part credible.

### 11.2 Distributional equivalence

MATLAB and Python draw from different RNG streams, so the ensembles are compared, not matched:
5/50/95 percentiles from a 20,000-draw MATLAB run and a 20,000-draw Python run must agree within the
Monte-Carlo standard error of the percentile (checked as |Δ| < 3 × MCSE for each of the three), on all
four demo series (`castaneda2010`, `lopes_santos2010`, `shevenell2011`, `wilsonlake`) × both runs ×
both modes.

### 11.3 Cross-implementation

Against `brews/baysparpy` (installed in the test env): percentiles must agree within MC error for
standard mode. Analogue mode is expected to agree as well (§5) — if it does not, that is a finding
worth reporting upstream, and the test records the discrepancy rather than asserting agreement.

### 11.4 Regression / performance

A benchmark test asserts the closed-form path stays within a factor of 2 of the recorded timings in
§5 and that it performs no BLAS call (checked with `threadpoolctl` — the thread pool must show zero
usage during a prediction).

### 11.5 The translation record

`tests/test_port.py` holds one test per `PORTING.md` entry, named for its ID
(`test_BT_09`, `test_BTA_05`, …), and `tests/test_porting_doc.py` asserts that every ID in the
document has such a test and that every `DEVIATION` and `DEFECT` entry also appears in the
document's register. The record cannot fall out of step with the code without the build going red.

Until the golden files have been generated in a MATLAB session, every golden test is marked
`xfail(strict=False)` with the reason "golden file not yet generated from MATLAB" — so the suite says
what has and has not been checked against the reference, rather than passing vacuously.

### 11.6 Smoke and API

Mirroring TEXAS's suite: `test_imports.py`, `test_public_api_docs.py` (every public symbol has a
docstring), `test_optional_deps.py` (core import works without cartopy/cmdstanpy), plus the CmdStan
parity and interop tests of §9.

---

## 12. Forward path into `texas-psm`

Once the TEXAS paper is published, BAYSPAR becomes one proxy system model among several. To make that
a move rather than a rewrite:

| BAYSPARpy | TEXAS counterpart | Alignment action |
|---|---|---|
| `predict_seatemp(tex, ...)` | `predict_T_from_proxyObs(proxyObs, prior_mu_t, prior_sigma_t, ...)` | Keep `prior_mean`/`prior_std` as aliases of `prior_mu_t`/`prior_sigma_t`; return the same dict keys (`x_vals`, `p5`, `p50`, `p95`, `ensemble`, `metadata`). |
| `Prediction` | TEXAS returns dicts/`xr.Dataset` | `Prediction.to_dict()` emits the TEXAS shape exactly. |
| `_cmdstan/` | `TEXAS/utils/{paths,install,doctor}.py` | Vendored copy, parity-tested (§9) — deleting it and importing TEXAS's is a one-line change. |
| `modelparams/registry.py` | `TEXAS/utils/download.py` | Same registry shape, same cache-dir env-var convention. |
| Posterior NetCDF | TEXAS posterior `.nc` | Same dimension names and provenance attributes. |
| `stan_models/*.stan` | `TEXAS/stan_models/*.stan` | Same header-comment conventions, same `reduce_sum`/`grainsize` pattern. |

The target end state is `texas-psm` exposing BAYSPAR as a registered PSM, with BAYSPARpy remaining
installable standalone for people who only want TEX₈₆.

---

## 13. Milestones

**Phase 0 — scaffold (≈1 day).** Repo, `pyproject.toml`, CI (lint + tests, three OSes), `CLAUDE.md`,
licence and citation files, `tools/convert_modeloutput.py`, the NetCDF store tracked in Git LFS.

**Phase 1 — the port (≈3 days).** *Substantially done — the package is in `BAYSPARpy/`.*
`distance`, `stores`, `bayspar_tex`, `bayspar_tex_analog`, `tex_forward`, `Prediction`, the modern
API and the error types are written, with 65 passing tests: one per `PORTING.md` entry, the
reference-trace check, the record/code meta-test, and the no-BLAS check.
Still outstanding in this phase: plotting (`predictplot`, `analogmap`, `densityplot`), the demo
notebooks, and the three test files that need a MATLAB session or `brews/baysparpy` installed
(§11.1–11.3). The package currently reads `ModelOutput/` directly, which is why it lives in this
repository for now; Phase 0's conversion is what frees it to move. **This is the deliverable that replaces the current workflow** and
is worth cutting a `0.1.0` prototype tag at.

**Phase 2 — toolchain (≈1 day).** `_cmdstan/` vendored, `doctor`, `install_cmdstan`, console scripts,
parity and interop tests, installation docs including the conda and Docker paths.

**Phase 3 — Stan inversion (≈2 days).** `invT_bayspar_marginal.stan` and the AR(1) variant, the
`backend="stan"` path, the marginalisation comparison write-up.

**Phase 4 — refit (≈1 week, uncertain).** TT14 model reconstruction, `bayspar_spatial_calibration.stan`,
the `calibrate()` API, and the reproduction check of §10.3. Uncertain because it depends on details of
TT14 that are not in this repository; the first task in the phase is to read the paper against §3.3
and write down the model, before any Stan is written.

**Phase 5 — docs and release (≈2 days).** Jupyter Book, the three demo notebooks, migration notes for
`brews/baysparpy` users, and — only if Jess and Brewster agree — a PyPI release.

---

## 14. Open decisions

### Still open

1. **PyPI name.** `bayspar`, or inherit `baysparpy` from Brewster? Needs Jess and Brewster. Nothing
   is published until then; the import name `baysparpy` is safe under either outcome.
2. **Licence.** The MATLAB repo's licence carries over; confirm with Jess that it permits
   redistributing the parameter files under the new repository.
3. **TT14 model details for the refit** (§10.3): kernel family, whether α and β share a range
   parameter, and the exact priors. Requires reading the paper; ideally confirmed with Jess, who has
   the original sampler.
4. **Two entries in the `PORTING.md` register are Jess's call**, not a maintainer's: `BTA-05` (should
   the default pair τ² with the α, β of the same draw, against MATLAB's behaviour?) and `STO-01`
   (storing one parameter set rather than two bit-identical ones).

### Settled

| | Decision |
|---|---|
| Default `n_draws` | **1000**, following Jess — not `brews/baysparpy`'s 5000. It is also the bundled tier (§8), so the default result needs no data fetch. |
| Arctic warning | **None.** The calibration excludes coretops north of 70° N, but a *paleo* site above that latitude is a different situation, and the port will not second-guess it. MATLAB does not warn; neither will we. |
| Data hosting | **No Zenodo deposit.** Bundled tier in the wheel, full store in the repository via Git LFS (§8). A hosted tier, if ever needed, goes under PaleoLipidRR with the same mechanism TEXAS uses. |

---

## References

Tierney, J. E., & Tingley, M. P. (2014). A Bayesian, spatially-varying calibration model for the
TEX₈₆ proxy. *Geochimica et Cosmochimica Acta*, 127, 83–106. https://doi.org/10.1016/j.gca.2013.11.026

Tierney, J. E., & Tingley, M. P. (2015). A TEX₈₆ surface sediment database and extended Bayesian
calibration. *Scientific Data*, 2, 150029. https://doi.org/10.1038/sdata.2015.29
