# Handoff

Where the BAYSPARpy work stands, and what to do next on your own machine.
Last updated 2026-09-12, at commit `9f2ab49` plus the golden-test scaffolding.

---

## Pick it up

```bash
git clone https://github.com/PaleoLipidRR/BAYSPAR
cd BAYSPAR
git checkout claude/cool-hopper-11jy6s     # or master, if this has been merged

cd BAYSPARpy
python -m venv .venv && source .venv/bin/activate     # or conda, as you prefer
pip install -e ".[dev]"
python -m pytest                                       # expect 65 passed, 13 skipped
```

Python ≥ 3.10, numpy and scipy; `[dev]` adds pytest and threadpoolctl. The package finds the
MATLAB `ModelOutput/` directory by walking up from its own location, so a checkout of this
repository is all it needs; `BAYSPAR_MODELOUTPUT` overrides that if you move the data.

Quick check that it works on something real:

```python
import numpy as np, scipy.io as sio, baysparpy as bp
rec = sio.loadmat("../ModelOutput/tex_testdata.mat")["lopes_santos2010"][0, 0]
out = bp.bayspar_tex(np.ravel(rec["tex86"]), lon=-17.6635, lat=9.166,
                     prior_std=6.0, runname="subT", seed=0)
out.prior_mean      # 19.1673884681
out.preds[:3]       # 5/50/95, degrees C
```

---

## The one job that needs your machine: golden files

Everything currently passing checks the port against itself or against intermediates the port
produced. **Nothing has yet been compared against MATLAB.** A shared misreading of the reference
would sail through. Closing that is one MATLAB session:

```matlab
>> cd /path/to/BAYSPAR
>> addpath('docs/BAYSPARpy/audit')
>> generate_golden
```

It writes `docs/BAYSPARpy/audit/golden_matlab.mat` and prints the same values in the layout of
`audit/reference_trace.txt`, so you can eyeball the two side by side before trusting either. Then:

```bash
cd BAYSPARpy && python -m pytest tests/test_golden.py -v
```

Those 11 tests currently skip with the reason "golden file not yet generated from MATLAB". Once the
`.mat` exists they run: thinning indices, chordal distances, prior means and their observation
counts, grid-cell lookup, the Wilson Lake analogue set, per-draw `post_mean`/`post_sig`, and the
`prctile` convention — all to 1e-10 or exact. Commit the `.mat` (it is a few KB) so CI can run them
too.

`generate_golden.m` needs base MATLAB plus the Statistics Toolbox (for `prctile`). It has not been
executed — it was written from the source, not tested — so if it throws, the fix is likely a
one-liner and the values it computes are all mirrored in `audit/reference_trace.py` if you would
rather cross-check by hand.

**If a golden test fails**, the failing quantity localises the problem to a single `PORTING.md`
entry — that is what the ID scheme is for. Quote the ID in the issue.

---

## What is done

| | |
|---|---|
| `docs/BAYSPARpy/SPEC.md` | Scope, naming, package layout, data plan, CmdStan and Stan plan, validation plan, milestones, decisions |
| `docs/BAYSPARpy/PORTING.md` | 45-entry translation record; 10 in the register; every entry has a test |
| `docs/BAYSPARpy/audit/` | Three Python scripts + `reference_trace.txt`, and now `generate_golden.m` |
| `BAYSPARpy/src/baysparpy/` | `distance`, `stores`, `bayspar_tex`, `bayspar_tex_analog`, `tex_forward`, `results`, `predict`, `errors`, `constants`, `utils` |
| `BAYSPARpy/tests/` | 65 passing: one per PORTING.md entry, the reference-trace check, the record↔code meta-test, the no-BLAS check. 13 skipped: 11 golden, 2 Phase 0 |

Two findings the code produced that the documents could not:

- **`BTA-05` is free.** MATLAB pairs τ² with α, β from different draws in analogue mode. On Wilson
  Lake at 5000 draws the difference between that and correct pairing is 0.024/0.011/0.033 °C across
  the three percentiles — against 0.027/0.012/0.034 °C between two runs of the *same* pairing at
  different seeds. Monte-Carlo noise. Jess can decide the default knowing the fix costs nothing.
- **Seven of the 903 SST coretops have a target error SD of exactly zero** (subT has none). The
  errors-in-variables refit divides by that. Pinned by `test_STO_04`; caveat in SPEC §10.3.

## What is not done

In the order I would do it:

1. **Golden files** (above). Everything else is less valuable until the port is validated.
2. **`test_vs_baysparpy.py`** — `pip install baysparpy` and compare percentiles. Needs both packages
   in one environment, which is why the import name is `baysparpy` and not `bayspar`.
3. **Plotting** — `predictplot`, `analogmap`, `densityplot`, matching the two MATLAB demo figures.
4. **The two demo notebooks**, replacing `Demo_StandardPrediction.m` and `Demo_AnalogPrediction.m`.
5. **Phase 0 leftovers** — `tools/convert_modeloutput.py` (`.mat` → NetCDF with provenance) and the
   demo CSV export. These are what free the package to leave this repository; `test_STO_02` and
   `test_STO_05` skip until they exist.
6. **Phase 2, CmdStan** — vendor TEXAS's `paths.py`/`install.py`/`doctor.py` into
   `baysparpy/_cmdstan/`, plus the parity and interop tests (SPEC §9).
7. **Phase 3, the Stan inversion**, then **Phase 4, the refit** (SPEC §10).

## Open questions

- **PyPI name** — `bayspar` (reserved, unpublished) or inherit `baysparpy` from Brewster. Needs him
  and Jess. Nothing is published either way.
- **Licence** — confirm with Jess that redistributing the parameter files under a new repository is
  fine.
- **`BTA-05` and `STO-01`** — Jess's call, per the register. The measurement above should make
  `BTA-05` easy.
- **TT14 kernel and priors** for the refit — a reading task before any Stan is written.

## When you split the repository out

`BAYSPARpy/` becomes the root of `PaleoLipidRR/BAYSPARpy`. Three things move or change with it:

- `docs/BAYSPARpy/*.md` and `audit/` move too — `tests/test_porting_doc.py` and `test_golden.py`
  walk up looking for `docs/BAYSPARpy/`, and skip with an explanation if it is absent.
- `stores.find_modeloutput()` will no longer find `ModelOutput/` by walking up. That is what item 5
  above replaces; until then, `BAYSPAR_MODELOUTPUT` bridges it.
- `git subtree split -P BAYSPARpy -b bayspar-py-only` preserves the history if you want it.
