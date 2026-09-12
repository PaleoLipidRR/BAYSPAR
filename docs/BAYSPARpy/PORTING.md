# MATLAB → Python translation record

This is the document to review. It walks every line of the four MATLAB functions and states what
the Python does instead, so that a reviewer can check the recode **without reading the Python
source**, and can disagree with a specific numbered decision rather than with the port as a whole.

Companion documents: [`SPEC.md`](SPEC.md) (what is being built and why),
[`audit/reference_trace.txt`](audit/reference_trace.txt) (deterministic values to compare against
MATLAB).

**The code these entries describe is in [`../../BAYSPARpy/`](../../BAYSPARpy/)**, developed in this
repository so its tests can read `ModelOutput/` directly, and to be split out to
`PaleoLipidRR/BAYSPARpy` once the port is validated. Run it with
`cd BAYSPARpy && PYTHONPATH=src python -m pytest`.

---

## How to review this

**Every entry has an ID** (`BT-07`, `TXF-04`, …). Quote the ID when you disagree — in a review
comment, an issue, or a note in the margin. IDs are stable: they never get renumbered, and a
withdrawn decision keeps its ID with a `WITHDRAWN` status.

**Every entry has a status.** Only two of the five need your judgement:

| Status | Meaning | Needs review? |
|---|---|---|
| `EXACT` | Same arithmetic in the same order. Results agree bit-for-bit, or to 1e-12 where an expression is re-associated. | No — just check the translation is what you'd write. |
| `EQUIV` | Different expression, same result, pinned by a test. | Skim. |
| `EXTENSION` | New capability. Defaults reproduce MATLAB exactly; the behaviour only changes if a caller opts in. | Skim. |
| **`DEVIATION`** | **Intentionally different behaviour.** | **Yes.** |
| **`DEFECT`** | **MATLAB behaviour we believe is wrong.** `mode="matlab"` reproduces it; the default does not. | **Yes.** |

**Every entry names its test.** Test IDs follow the entry ID: `BT-07` is pinned by
`tests/test_port.py::test_BT_07`. A CI check asserts that every ID in this document has a
corresponding test and that every `DEVIATION` and `DEFECT` also appears in the register below — so
this document cannot drift away from the code without the build going red.

**Line numbers** refer to the MATLAB files as committed in this repository at the time of writing
(commit `22e997f`).

---

## Register: everything that does not behave like MATLAB

Ten entries. If you only read one part of this document, read this table.

| ID | Status | What changes | Why | Decides |
|---|---|---|---|---|
| [`BTA-05`](#bta-05) | DEFECT | τ² is paired with the α, β of the **same** draw in analogue mode | MATLAB's `repmat`/`reshape` pair them by different indices, breaking the joint posterior | Jess |
| [`BTA-08`](#bta-08) | DEFECT | The analogue ensemble is returned as `(N, n_analogues, n_draws)` | The MATLAB docstring promises that shape; the code returns 2-D | Ronnakrit |
| [`BT-08`](#bt-08) | DEFECT | A site on a cell boundary resolves to the nearest centroid, with a warning | MATLAB errors for N > 1 and silently duplicates rows for N = 1 | Ronnakrit |
| [`TXF-07`](#txf-07) | DEVIATION | Forward-model output size follows `n_out`, default 1000 | MATLAB hard-codes 1000 columns regardless of the ensemble used | Ronnakrit |
| [`GEN-03`](#gen-03) | DEVIATION | Explicit seeded RNG per call | MATLAB uses the global stream; results are not reproducible run to run | Ronnakrit |
| [`BT-05`](#bt-05) | EXTENSION | `max_dist` and `min_num` are keyword arguments | Hard-coded in MATLAB at 500 km / 1, which stay the defaults | — |
| [`BT-01`](#bt-01) | EXTENSION | Keyword arguments instead of `varargin` position | Same defaults; `n_draws` validated against the real ensemble length | — |
| [`TXF-01`](#txf-01) | EXTENSION | Same, for the forward model | Keeps MATLAB's `(lat, lon, t)` order and its error messages | — |
| [`STO-01`](#sto-01) | EXTENSION | `params_analog` is derived from the standard store by cell index | The two files are bit-identical (audited); storing one removes 97 MB | Jess |
| [`TXF-05`](#txf-05) | EXTENSION | Forward model accepts `n_draws` thinning | MATLAB always uses all 20,000 draws, then subsamples 1000 — which stays the default | — |

Two further things reviewers often expect to be deviations and are **not**: the `mode="matlab"`
switch reproduces `BTA-05` and `BTA-08` exactly for anyone validating against MATLAB output, and
`brews/baysparpy`'s analogue mode — which reads α, β from the *standard* parameter file — turns out
to be numerically identical to MATLAB's, for the reason in [`STO-01`](#sto-01).

---

## Conventions applied everywhere

<a id="gen-01"></a>
### GEN-01 · EXACT · Index base

MATLAB indices are 1-based and inclusive; Python's are 0-based and half-open. Every index that
crosses the boundary (thinning vectors, cell indices, `Inds_Stack`) is converted at the point it
enters or leaves the store, never silently. Stored index arrays keep MATLAB's 1-based convention on
disk and are converted on load, so a reviewer comparing a `.nc` file against a `.mat` file sees the
same numbers.

<a id="gen-02"></a>
### GEN-02 · EXACT · Column vectors

```matlab
dats = dats(:);          % bayspar_tex.m:82, bayspar_tex_analog.m:65
t = t(:); lat = lat(:); lon = lon(:);   % TEX_forward.m:35-37
```
```python
dats = np.asarray(dats, dtype=float).ravel()
```
`.ravel()` on a 1-D or single-column input is exactly `(:)`. For a genuinely 2-D input MATLAB
flattens column-major and NumPy row-major, so the port **raises** on anything with more than one
non-trivial dimension rather than flattening it differently. MATLAB would accept it silently.

<a id="gen-03"></a>
### GEN-03 · DEVIATION · Random numbers

MATLAB draws from the global stream (`randn`, `normrnd`, `randsample`), so results change run to
run and cannot be reproduced. Every Python entry point takes `seed=None` and builds
`np.random.default_rng(seed)`; the generator is threaded through explicitly rather than taken from
the global NumPy state. `seed=None` keeps MATLAB's "different every time" behaviour as the default.

Equivalence with MATLAB is therefore **distributional**, never elementwise — see the test plan in
SPEC §11. Everything upstream of the first random draw is compared exactly, which is what
`audit/reference_trace.txt` is for.

<a id="gen-04"></a>
### GEN-04 · EXACT · Types

All arithmetic is float64. `Locs_Comp` and `Data_Input.Locs` are int16 on disk and are cast to
float before any comparison with a float coordinate; `Inds_Stack` is uint8 and is cast to int
before indexing. (uint8 arithmetic saturates in both languages, which is exactly the kind of thing
that produces a wrong answer with no error.)

<a id="gen-05"></a>
### GEN-05 · EQUIV · `repmat` → broadcasting

MATLAB materialises `repmat(x, Nd, 1)`; NumPy broadcasts. Same values, no allocation. Shapes are
annotated in the source at each step, because a broadcasting mistake here produces a plausible
wrong answer rather than an error — see [`BT-09`](#bt-09).

<a id="gen-06"></a>
### GEN-06 · EQUIV · `length()`

`length(x)` is `max(size(x))`, which for a 162 × 20000 matrix is 20000, not 162. The port never
translates it literally; each use is replaced by the dimension the surrounding code actually means,
named explicitly (`n_cells`, `n_draws`, `n_obs`).

<a id="gen-07"></a>
### GEN-07 · EQUIV · Errors

MATLAB `error(...)` strings are preserved verbatim as the message of a typed Python exception, so
that a user searching for "Your search tolerance is too narrow" finds the same text.

| MATLAB message | Python exception |
|---|---|
| `Your search tolerance is too narrow` | `SearchToleranceError(ValueError)` |
| `To use analog mode, enter a search tolerance in TEX units` | `ValueError` |
| `please enter "analog" to specify analog mode` | `ValueError` |
| `You entered too many or too few arguments` | replaced by the signature itself |

---

## `EarthChordDistances_2.m` → `baysparpy/distance.py`

<a id="ecd-01"></a>
### ECD-01 · EXACT · Earth radius

```matlab
RR = 6378.137;   % line 10
```
```python
EARTH_RADIUS_KM = 6378.137   # WGS-84 equatorial radius, as in the MATLAB source
```

<a id="ecd-02"></a>
### ECD-02 · EQUIV · The pairing

```matlab
Pts_paired_Vec = [kron(llPoints1, ones(M,1)), kron(ones(N,1), llPoints2)];   % line 17
```
The Kronecker products build an N·M × 4 matrix of every point pair. NumPy broadcasts
`points1[:, None, :]` against `points2[None, :, :]` instead — same pairs, without materialising
N·M rows. This matters: the prior-mean search pairs one site against 37,686 WOA cells, and a
gridded application would pair thousands.

<a id="ecd-03"></a>
### ECD-03 · EXACT · The formula

```matlab
Half_Angles = asin(sqrt( sin((lat1-lat2)*pi/180/2).^2 ...
              + cos(lat1*pi/180).*cos(lat2*pi/180).*sin(abs(lon1-lon2)*pi/180/2).^2 ));
Chords = 2*RR*sin(Half_Angles);                                           % lines 18-19
```
Transcribed term for term, including the operation order, so rounding matches.

Two notes for the reviewer rather than changes:
- `abs()` on the longitude difference is a no-op — `sin(x)²` is even — but it is kept, because
  removing it would change nothing and cost a line of justification.
- This is the **chord** through the sphere, not the great-circle arc: `2R·sin(a/2)` rather than
  `R·a`. At the 500 km prior-mean cutoff the chord is 0.13 km (0.026%) shorter than the arc. The
  cutoff is a modelling choice, not a measurement, so this is immaterial — but the function is
  named for what it computes and the port keeps that name.

<a id="ecd-04"></a>
### ECD-04 · EXACT · Output orientation

```matlab
Chord_Dists_Mat = reshape(Chords_as_vec, M, N)';    % line 22
```
Result is N × M, entry (i, j) the distance from `points1[i]` to `points2[j]`. The Python returns
the same orientation, and the docstring says so, because the transpose at the end of the MATLAB is
easy to lose.

<a id="ecd-05"></a>
### ECD-05 · EXACT · Degenerate inputs

Identical points give exactly 0 (checked: 180°E to 180°W is 0.000000000 km). A single point against
a single point returns a 1 × 1 array, not a scalar, matching MATLAB.

---

## `bayspar_tex.m` → `baysparpy/bayspar_tex.py`

<a id="bt-01"></a>
### BT-01 · EXTENSION · Argument handling

```matlab
ng = nargin;                          % lines 48-59
if ng == 7,  Nsamps = varargin{1}; ens_sel = varargin{2};
elseif ng == 6, Nsamps = varargin{1}; ens_sel = 0;
elseif ng == 5, Nsamps = 1000; ens_sel = 0; end
```
```python
def bayspar_tex(dats, lon, lat, prior_std, runname="SST", *,
                n_draws=1000, save_ensemble=False, mode="modern",
                seed=None, strict_grid=False, max_dist=500.0, min_num=1):
```
Same defaults (1000 draws, no ensemble). Two differences, both additive: `n_draws` is validated
against the actual ensemble length (20,000 — the MATLAB header's "cannot exceed 15000" is stale,
see SPEC §3) and raises `EnsembleSizeError` rather than indexing out of range; and an `ng` outside
5–7 leaves `Nsamps` undefined in MATLAB, producing a confusing downstream error, where Python
raises at the call.

<a id="bt-02"></a>
### BT-02 · EQUIV · Loading parameters

```matlab
load(['ModelOutput/', 'Output_SpatAg_', runname, '/params_standard'], ...
     'alpha_samples_comp','beta_samples_comp','tau2_samples','Locs_Comp');   % lines 62-63
```
Replaced by a lazily-loaded, cached store (SPEC §8). Two behavioural consequences worth naming:
MATLAB resolves `ModelOutput/` relative to the **current working directory**, so the demos only run
from the repository root; the Python store resolves against the package and the cache directory, so
it works from anywhere. And MATLAB re-reads ~49 MB on every call; the Python store reads once per
process.

<a id="bt-03"></a>
### BT-03 · EQUIV · `runname`

```matlab
if runname=="SST" ... elseif runname=="subT" ...      % lines 66-70
```
An unrecognised `runname` falls through both branches in MATLAB, leaving `locs_st_obs` undefined
and failing later with an unrelated message. Python validates against `{"SST", "subT"}` up front and
raises naming both. Lower-case `"sst"`/`"subt"` are accepted as aliases, for compatibility with
`brews/baysparpy`'s `temptype` argument.

<a id="bt-04"></a>
### BT-04 · EXACT · Grid spacing

```matlab
grid_half_space = 10;      % line 73 — "grid spacing is hard-coded here"
```
`constants.GRID_HALF_SPACE = 10.0`. Not a parameter: it is a property of the fitted model, and the
162-cell tiling in the parameter file depends on it.

<a id="bt-05"></a>
### BT-05 · EXTENSION · Prior-mean search radius

```matlab
min_num = 1;      % line 76
max_dist = 500;   % line 79
```
Hard-coded in MATLAB, with a comment block explaining that `min_num=K, max_dist=0` uses the closest
K points and that the paper used `min_num=1, max_dist=500`. The port exposes both as keyword
arguments with exactly those defaults, so the documented alternative is reachable without editing
the source.

<a id="bt-06"></a>
### BT-06 · EXACT · Thinning

```matlab
ind_s = round(linspace(1, length(tau2_samples), Nsamps));   % line 86
alpha_samples_comp = alpha_samples_comp(:, ind_s);          % lines 87-89
```
```python
def _matlab_thinning(n_total: int, n_draws: int) -> np.ndarray:
    # MATLAB round() is half-away-from-zero; np.round is half-to-even.
    return np.floor(np.linspace(1, n_total, n_draws) + 0.5).astype(int) - 1
```
The comment in the MATLAB ("so as to use the full span of the ensemble even if few samples are
used") is the point of the whole line, and two existing ports lose it: the Colab notebook uses a
step-`arange` (and flags the difference itself), `brews/baysparpy` takes the first `nens` draws.
For `Nsamps=1000` the indices are 1, 21, 41, … 20000 — see `reference_trace.txt` §A, which tabulates
them for five values of `Nsamps` including the rounding-sensitive 999 and 1001.

Note the indices are **not** guaranteed distinct for `Nsamps` near `n_total`; MATLAB does not
deduplicate and neither does the port.

<a id="bt-07"></a>
### BT-07 · EXACT · Prior mean

```matlab
dists_prior = EarthChordDistances_2([lon, lat], locs_st_obs);        % line 112
[vals_dist_prior, inds_dist_prior] = sort(dists_prior);              % line 114
num_below_dist = find(vals_dist_prior < max_dist, 1, 'last');        % line 116
if num_below_dist > min_num
    prior_mean = mean(st_obs_ave_vec(inds_dist_prior(1:num_below_dist)));
else
    prior_mean = mean(st_obs_ave_vec(inds_dist_prior(1:min_num)));   % lines 119-124
end
```
```python
d = earth_chord_distances((lon, lat), locs_st_obs)[0]
order = np.argsort(d, kind="stable")                 # stable: ties order as MATLAB's sort does
n_below = int(np.searchsorted(d[order], max_dist, side="left"))   # strict <, on sorted distances
take = n_below if n_below > min_num else min_num
prior_mean = float(st_obs_ave_vec[order[:take]].mean())
```
Three subtleties, all pinned by tests:
- `find(..., 1, 'last')` on a sorted vector is a *count*, not a position — it happens to be both,
  which is why `side="left"` (strict `<`) is the right search.
- When nothing is within `max_dist`, MATLAB's `find` returns `[]`, the comparison `[] > 1` is
  `false`, and the branch falls to the single closest observation. The Python reproduces the
  fallback explicitly rather than relying on empty-comparison semantics.
- MATLAB `sort` is stable, so tied distances keep their original order; `kind="stable"` is required,
  not cosmetic, because averaging a different tied subset changes the prior mean.

Expected values for all three demo series and both runs are in `reference_trace.txt` §C — e.g.
`lopes_santos2010` under subT gives `PriorMean = 19.1673884681` from 56 observations within 500 km.

<a id="bt-08"></a>
### BT-08 · DEFECT · Grid-cell lookup

```matlab
inder_g = find(abs(Locs_Comp(:,1)-lon) <= grid_half_space ...
             & abs(Locs_Comp(:,2)-lat) <= grid_half_space);    % line 133
```
The comparison is inclusive at both edges, and `Locs_Comp` tiles the globe, so a site exactly on a
cell boundary matches **two** cells (a corner, four). Boundaries are at even multiples of 10°
longitude and odd multiples of 10° latitude — lon = 0°, lat = 10° and similar are not exotic
positions.

What MATLAB then does: `alpha_samples_comp` becomes 2 × M, and at line 149
`dats - repmat(alpha_samples_comp, Nd, 1)` is a dimension error for N > 1 — but for **N = 1** it
broadcasts and returns one row per matched cell, silently. `brews/baysparpy` has the same defect via
`.squeeze()`.

The port resolves the tie to the nearest centroid by chordal distance and issues a `UserWarning`
naming both cells; `strict_grid=True` raises `AmbiguousGridCellError` instead. `mode="matlab"` does
not reproduce the silent duplication — there is no result worth reproducing.

<a id="bt-09"></a>
### BT-09 · EXACT · The posterior

```matlab
pmu      = repmat(ones(Nd,1)*prior_mean, 1, Nsamps);                      % line 144
pinv_cov = repmat(prior_std, Nd, Nsamps).^-2;                             % line 145
sigmaS   = sqrt(tau2_samples);                                            % line 146
post_mean_num = pinv_cov .* pmu ...
              + repmat(sigmaS,Nd,1).^-2 .* repmat(beta_samples_comp,Nd,1) ...
                .* (dats - repmat(alpha_samples_comp,Nd,1));              % line 149
post_mean_den = pinv_cov + repmat(beta_samples_comp,Nd,1).^2 .* repmat(sigmaS,Nd,1).^-2;
post_mean = post_mean_num ./ post_mean_den;                               % line 151
post_sig  = sqrt(post_mean_den.^-1);                                      % line 152
```
```python
pinv_cov = prior_std ** -2.0                       # scalar; MATLAB tiles it to (N, M)
den = pinv_cov + beta ** 2 / tau2                  # (M,)   -- beta, tau2 are (M,)
num = pinv_cov * prior_mean + (beta / tau2) * (dats[:, None] - alpha)   # (N, M)
post_mean = num / den                              # (N, M) by broadcasting
post_sd = np.sqrt(1.0 / den)                       # (M,)
```
This is the conjugate normal update: precision-weighted mean of the prior and the likelihood. Three
things to check when reviewing:
- `sigmaS.^-2` is `1/tau2` — MATLAB takes the square root at line 146 and squares it back at 149.
  The port divides by `tau2` directly. Same value to 1e-16; the round trip through `sqrt` is the
  only reason this is `EXACT` rather than bit-identical, and the test allows 1e-12.
- `post_mean_den` does not depend on `dats`, so it is (M,) not (N, M). MATLAB computes the full
  (N, M) tile; the port broadcasts. This is most of the memory saving.
- `dats[:, None]` is the only place N enters. If that indexing is wrong the result is still
  plausible-looking, which is why `reference_trace.txt` §E prints `post_mean` and `post_sig` for
  `dats(1)` across the first five thinned draws — compare those five numbers first.

<a id="bt-10"></a>
### BT-10 · EQUIV · The draw

```matlab
Preds = post_mean + randn(Nd, Nsamps) .* post_sig;    % line 153
```
```python
preds = post_mean + rng.standard_normal(post_mean.shape) * post_sd
```
Different stream ([`GEN-03`](#gen-03)), same distribution. Note the draw is **independent per
sample and per draw** — there is no temporal structure in the prior, which is what makes the closed
form diagonal and the dense solve in `brews/baysparpy` unnecessary (SPEC §5).

<a id="bt-11"></a>
### BT-11 · EXACT · Percentiles

```matlab
Output_Struct.Preds = prctile(sort(Preds,2), [5 50 95], 2);    % line 156
```
```python
preds_pct = np.percentile(preds, [5, 50, 95], axis=1, method="hazen").T
```
MATLAB `prctile` uses Hazen plotting positions `(i − 0.5)/n` with linear interpolation and clamping
at the extremes; NumPy's default (type 7) does not match and would bias the 5th and 95th toward the
middle. `reference_trace.txt` §F shows the three candidate methods against
`prctile(1:10, [5 50 95])` → `1.0, 5.5, 10.0`.

The inner `sort(Preds,2)` is a no-op — `prctile` sorts internally — and is dropped.

<a id="bt-12"></a>
### BT-12 · EQUIV · Output

```matlab
Output_Struct.Preds / .SiteLoc / .GridLoc / .PriorMean / .PriorStd / .PredsEns   % lines 95-103, 159-161
```
A `Prediction` dataclass with `preds`, `site_loc`, `grid_loc`, `prior_mean`, `prior_std`, and
`ensemble` (`None` unless `save_ensemble=True`), plus `metadata` recording `runname`, `n_draws`,
`mode`, `seed`, the store version and the thinning indices — so a result carries enough provenance
to be reproduced. `.to_dict()` emits the MATLAB field names for anyone diffing against a saved
`Output_Struct`.

---

## `bayspar_tex_analog.m` → `baysparpy/bayspar_tex_analog.py`

Lines 46–57 (`nargin`), 65 (`dats(:)`), 68–71 (thinning) and 116–128 (the solve and percentiles)
are identical to their `bayspar_tex.m` counterparts and are covered by
[`BT-01`](#bt-01), [`GEN-02`](#gen-02), [`BT-06`](#bt-06), [`BT-09`](#bt-09), [`BT-11`](#bt-11).
What follows is what differs.

<a id="bta-01"></a>
### BTA-01 · EQUIV · Parameter source

```matlab
load(['ModelOutput/', 'Output_SpatAg_', runname, '/params_analog'], ...
     'alpha_samples','beta_samples','tau2_samples');            % lines 61-62
load(['ModelOutput/Data_Input_SpatAg_', runname], 'Data_Input');  % line 63
```
The port reads the standard store and selects the 80 data-bearing cells by index. This is not an
approximation — see [`STO-01`](#sto-01), where the two files are shown to be bit-identical.

<a id="bta-02"></a>
### BTA-02 · EXACT · Cell-mean TEX₈₆

```matlab
for i = 1:N_bg
    spatialMean(i) = mean(Data_Input.Obs_Stack(Data_Input.Inds_Stack == i));   % lines 94-97
end
```
```python
sums = np.bincount(inds_stack, weights=obs_stack, minlength=n_cells + 1)[1:]
counts = np.bincount(inds_stack, minlength=n_cells + 1)[1:]
spatial_mean = sums / counts          # every cell has >= 1 observation (audited)
```
Equivalent because `Inds_Stack` covers 1..80 with no empty cell — checked by the audit script, so
the division is safe. Note this averages **observations**, not sites: a cell with eight cores at one
1° location and one at another weights the first location eight times. That is MATLAB's behaviour
and the port keeps it.

<a id="bta-03"></a>
### BTA-03 · EXACT · Analogue selection

```matlab
inder_g = spatialMean >= (mean(dats)-search_tol) & spatialMean <= (mean(dats)+search_tol);  % line 99
```
Inclusive at both ends, comparison against the mean of the whole input series. Reproduced exactly.
For the Wilson Lake demo (`search_tol = std(dats)*2`), 24 of the 80 cells are selected — listed with
their centroids and cell means in `reference_trace.txt` §D, which is the fastest way to confirm this
line ported correctly.

Note `std(dats)` in the demo normalises by N−1; NumPy's default is N. The port's demo notebook uses
`ddof=1`.

<a id="bta-04"></a>
### BTA-04 · EXACT · Empty selection

```matlab
if sum(inder_g)==0, error('Your search tolerance is too narrow'), else, end   % lines 100-103
```
Same message, as `SearchToleranceError` ([`GEN-07`](#gen-07)). The empty `else` branch is dropped.

<a id="bta-05"></a>
### BTA-05 · DEFECT · τ² is mis-paired with α and β

```matlab
tau2_samples  = repmat(tau2_samples, 1, size(alpha_samples,1));     % line 107
alpha_samples = reshape(alpha_samples, 1, n_an*M);                  % line 109
beta_samples  = reshape(beta_samples,  1, n_an*M);                  % line 110
```
`reshape` is column-major, so flat position *k* of α and β holds
`(location = (k−1) mod n_an, draw = ⌊(k−1)/n_an⌋)`. `repmat` tiles τ² by draw, so the same position
holds `draw = (k−1) mod M`. For `n_an > 1` those are different draws.

α and β stay correctly paired **with each other** — both are reshaped identically — so the
regression line at each flat position is a real posterior sample. Only its residual variance comes
from somewhere else in the chain.

**Measured impact: none detectable.** On the Wilson Lake demo at 5000 draws across its 24 analogue
cells, the mean absolute difference between the two pairings is 0.024 / 0.011 / 0.033 °C for the
5th / 50th / 95th percentiles — against 0.027 / 0.012 / 0.034 °C between two runs of the *same*
pairing at different seeds. The difference is Monte-Carlo noise. That is what the exchangeability
argument predicts: the pooled ensemble holds the same α, β and τ² values, merely re-paired, and τ²
is close to independent of α and β in this posterior.

It still breaks the joint posterior, and it is the exact constraint TEXAS's inverse model states in
capitals ("ALL PARAMETERS MUST USE THE SAME DRAW INDEX m"), so the default pairs correctly. The point
of measuring was to be able to say the fix costs nothing rather than to assume it.

```python
# modern: one (n_an, M) block per parameter, flattened identically -> draw index preserved
alpha = alpha[sel]                                   # (n_an, M)
beta = beta[sel]                                     # (n_an, M)
tau2_pooled = np.broadcast_to(tau2, (n_sel, n_draws))  # (n_an, M) -- same draw, every location
```
`mode="matlab"` reproduces the mis-pairing exactly. `tests/test_port.py::test_BTA_06` runs both on
the Wilson Lake demo at 20,000 draws and reports the difference in the 5/50/95 percentiles, so the
magnitude is on the record rather than assumed. **Jess should decide** which is the default.

<a id="bta-06"></a>
### BTA-06 · EXACT · Pooling order

The flattening order matters for [`BTA-08`](#bta-08) and for anyone reading the ensemble. MATLAB's
column-major reshape gives *location fastest*: positions 1..n_an are draw 1 at each location, then
draw 2, and so on. The port keeps the parameters as a 2-D `(n_analogues, n_draws)` block and never
flattens them, which makes the order a property of the array shape rather than of a convention.

<a id="bta-07"></a>
### BTA-07 · EXACT · Analogue locations

```matlab
Output_Struct.AnLocs = Data_Input.Locs(inder_g,:);    % line 112
```
`Prediction.analog_locs`, as `(lon, lat)` pairs in `Data_Input.Locs` order — the same order as the
ensemble's location axis, which is what makes [`BTA-08`](#bta-08) interpretable.

<a id="bta-08"></a>
### BTA-08 · DEFECT · Ensemble shape

The MATLAB header (line 41) and this repository's ReadMe both promise:

> `.PredsEns` — Nd by No. analog locations by Nsamps array of predictions … the second dimension
> corresponds to the locations in `.AnLocs`.

The code returns `Preds` of shape `Nd × (n_analogues · Nsamps)` — 2-D, with the location axis
flattened into the draw axis. Any downstream code that trusts the documentation and indexes the
second dimension as a location gets a draw instead.

The port returns `(N, n_analogues, n_draws)`, as documented, with axis 1 aligned to
`analog_locs`. `mode="matlab"` returns the flattened 2-D array for bit-comparison; the docstring
gives the reshape that converts one to the other:

```python
ens_3d = ens_2d.reshape(n_obs, n_draws, n_analogues).transpose(0, 2, 1)   # location fastest
```

---

## `TEX_forward.m` → `baysparpy/tex_forward.py`

<a id="txf-01"></a>
### TXF-01 · EXTENSION · Argument handling

```matlab
if ng==6, runname=varargin{1}; type=varargin{2}; stol=varargin{3};
elseif ng==5, error('To use analog mode, enter a search tolerance in TEX units');
elseif ng==4, runname=varargin{1}; type="standard";
elseif ng==3, runname='SST'; type="standard";
else, error('You entered too many or too few arguments'); end          % lines 12-27
```
```python
def tex_forward(lat, lon, t, runname="SST", type="standard", search_tol=None, *,
                n_draws=None, n_out=1000, seed=None):
```
Note the argument order is `(lat, lon, t)` — latitude first, the opposite of `bayspar_tex`. That is
MATLAB's order and the port keeps it rather than silently "fixing" it; the docstring flags it, and
`brews/baysparpy` made the same choice. `type="analog"` without `search_tol` raises with MATLAB's
message.

`type` shadows the Python builtin. It is kept for signature compatibility; internally the parameter
is bound to `mode_`.

<a id="txf-02"></a>
### TXF-02 · EQUIV · Mode validation

```matlab
if ~strcmp(type,"standard") && ~strcmp(type,"analog")
    error('please enter "analog" to specify analog mode')       % lines 30-32
end
```
Validated against `{"standard", "analog"}` before any file is read.

<a id="txf-03"></a>
### TXF-03 · EQUIV · Standard-mode cell lookup

```matlab
Nloc = length(lon);
for i=1:Nloc
    inder_g(i) = find(abs(Locs_Comp(:,1)-lon(i))<=grid_half_space & ...
                      abs(Locs_Comp(:,2)-lat(i))<=grid_half_space);   % lines 56-59
end
alpha_samples = alpha_samples_comp(inder_g, :);
```
Vectorised, with the same boundary tie-breaking as [`BT-08`](#bt-08) (here MATLAB's `find` returning
two indices would make the assignment `inder_g(i) = [...]` fail outright, which is at least loud).

**The multiple-location semantics are easy to misread**, and the ReadMe's "tex_forward can now
handle multiple lat/lon locations!" does not spell them out: line 89 pairs `t` element-wise with the
per-location parameter rows, so this is *one temperature per site*, not a time series at each of
several sites. The port requires `len(t) == len(lat) == len(lon)` when more than one location is
given, and raises a message saying exactly that.

<a id="txf-04"></a>
### TXF-04 · EXACT · Analogue search is in temperature units

```matlab
spatialMean(i) = mean(Data_Input.Target_Stack(Data_Input.Inds_Stack==i));   % line 70
inder_g = spatialMean >= (mean(t)-stol) & spatialMean <= (mean(t)+stol);    % line 73
```
`Target_Stack` — the **temperature** field — where `bayspar_tex_analog.m` uses `Obs_Stack`, the
TEX₈₆ field ([`BTA-02`](#bta-02)). This is correct in both cases (each searches in the units of its
own input) and is the single easiest line in the package to port wrongly, since the two functions
are otherwise near-identical. Separate tests assert each reads the field it should.

`search_tol` is therefore in °C here and in TEX₈₆ units in `bayspar_tex_analog`. The port names the
argument `search_tol` in both, as MATLAB does, and says so in both docstrings.

<a id="txf-05"></a>
### TXF-05 · EXTENSION · No thinning

Unlike the two prediction functions, the forward model does not thin: it uses all 20,000 draws and
then subsamples the *output*. The port keeps that as the default (`n_draws=None`), and accepts an
explicit `n_draws` to thin the parameters first, using the same rule as [`BT-06`](#bt-06).

The τ² mis-pairing of [`BTA-05`](#bta-05) is present here too (lines 81–84 are the same three
lines), and is handled the same way.

<a id="txf-06"></a>
### TXF-06 · EQUIV · The draw

```matlab
tex86 = normrnd(t .* beta_samples + alpha_samples, repmat(sqrt(tau2_samples), length(t), 1));  % line 89
```
```python
mu = t[:, None] * beta + alpha                      # (N, M)
tex = rng.normal(mu, np.sqrt(tau2))                 # tau2 broadcasts along N
```

<a id="txf-07"></a>
### TXF-07 · DEVIATION · Output size

```matlab
iters = length(alpha_samples);
randind = randsample(iters, 1000);      % without replacement
tex86 = tex86(:, randind);              % lines 91-93
```
MATLAB always returns exactly 1000 columns, whatever the ensemble size, and the ReadMe describes the
output as "N x 1000". The port makes this `n_out=1000` — same default, and a caller who wants the
full ensemble can ask for it. `randsample(n, k)` draws **without** replacement; the port uses
`rng.choice(n_total, size=n_out, replace=False)` and raises if `n_out > n_total` rather than
erroring inside the sampler.

<a id="txf-08"></a>
### TXF-08 · EXACT · Clipping

```matlab
tex86(tex86>1) = 1;
tex86(tex86<0) = 0;         % lines 95-96
```
Kept — TEX₈₆ is a ratio in [0, 1] and the linear forward model does not know that. The port also
records the clipped fraction in `metadata["clipped_fraction"]`: for cold subT draws it is not
negligible, and a user comparing forward-modelled to measured TEX₈₆ should know when a tenth of the
ensemble is pinned at a bound.

---

## Data files → the stores

<a id="sto-01"></a>
### STO-01 · EXTENSION · One parameter store, not two

Audited (`audit/audit_modeloutput.py`, section 2): for both SST and subT,

```
alpha_samples == alpha_samples_comp[Data_Input.Locs]   -> np.array_equal True
beta_samples  == beta_samples_comp[Data_Input.Locs]    -> np.array_equal True
tau2_samples  identical                                -> np.array_equal True
```

`params_analog.mat` is the 80 data-bearing rows of `params_standard.mat`, in `Data_Input.Locs`
order — not a separate fit. The port stores the 162-cell array once plus an 80-element index, which
removes 97 MB of the 145 MB store and the standing question of which file a given α came from.

**This is the one entry where a reviewer should check the audit rather than take the document's
word**, because everything about analogue mode rests on it. Run
`python docs/BAYSPARpy/audit/audit_modeloutput.py` and read section 2.

<a id="sto-02"></a>
### STO-02 · EQUIV · `.mat` → NetCDF

`tools/convert_modeloutput.py` is committed and deterministic: dimensions `(cell, draw)`,
coordinates `lon`/`lat`, and attributes recording the source `.mat` filename, its SHA-256, the
conversion date and the converter version. A reviewer can re-run it and diff the result. float64
throughout — float32 would be 44 MB instead of 97 MB but changes the last digits, and the point of
this package is to match a reference implementation.

<a id="sto-03"></a>
### STO-03 · EQUIV · Observation files

`obsSST.mat` / `obssubT.mat` carry `locs_st_obs` and `st_obs_ave_vec` — 37,686 and 37,105 1° WOA
cells respectively. Note the two runs have **different** location vectors (subT drops cells where
the 0–200 m gamma-weighted average is undefined), so the store keeps them per-run rather than
sharing one grid.

<a id="sto-04"></a>
### STO-04 · EQUIV · `Data_Input` → a tidy table

The nested 35 × 80 cell arrays become one flat table with columns
`cell_index, tex86, target_t, target_t_sd`, one row per observation (903 for SST, 906 for subT).

**Found while testing this:** seven of the 903 SST observations carry a target error standard
deviation of *exactly zero* (none of the subT ones do). An errors-in-variables likelihood divides by
that quantity, so the Stan refit (SPEC §10.3) has to floor or special-case them; the test pins the
count so the fix is driven by the data rather than by a crash three phases from now. The 35-element first axis is the **site slot within a
cell**, not a time or an ensemble member — see SPEC §3.3. This table is the input format for the
refit (SPEC §10.3), so "refit the original" and "refit with my coretops" become the same call.

<a id="sto-05"></a>
### STO-05 · EQUIV · Demo data

`tex_testdata.mat` and `wilsonlake.mat` ship as CSV, matching `brews/baysparpy`'s example files so
that comparison scripts can load one file for both packages.

---

## What the tests actually assert

| Test | Covers | Kind | Status |
|---|---|---|---|
| `test_port.py::test_<ID>` | one entry above, 45 of them | unit | 43 pass, 2 skip (`STO-02`, `STO-05` are Phase 0) |
| `test_reference_trace.py` | the package reproduces every value in `audit/reference_trace.txt` | golden-ish | passing |
| `test_porting_doc.py` | every ID here has a test and every test an ID; every DEVIATION and DEFECT is in the register; IDs have no gaps | meta | passing |
| `test_no_blas.py` | a prediction makes no BLAS call and allocates no N × N intermediate | performance | passing |
| `test_golden.py` | the same quantities against values generated **in MATLAB** | golden | **not written** — needs one MATLAB session |
| `test_distributional.py` | 5/50/95 from 20,000-draw runs, Python vs MATLAB, within 3 × the percentile's Monte-Carlo SE | statistical | **not written** — same dependency |
| `test_vs_baysparpy.py` | same, against `brews/baysparpy` in the same environment | cross-impl | not written |

`test_reference_trace.py` checks the package against values this port itself produced, so it catches
drift but cannot catch a shared misreading of the MATLAB. Only `test_golden.py` can, and it needs a
MATLAB session — until then the record says so rather than implying the port is validated against
the reference.

---

## Changelog of this document

| Date | Change |
|---|---|
| 2026-09-11 | First version: 4 MATLAB files, 45 entries, 10 in the register. No Python written yet — every entry describes intended behaviour. |
| 2026-09-12 | The port landed in `BAYSPARpy/`; 43 of 45 entries now have a passing test and 2 are skipped as Phase 0 work. Two entries changed as a result: `BTA-05` gains the measured size of the τ² mis-pairing (Monte-Carlo noise, on the Wilson Lake demo), and `STO-04` records seven SST coretops whose target error SD is exactly zero, which the refit will have to handle. |
