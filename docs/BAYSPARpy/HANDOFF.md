# Handoff

Where the BAYSPARpy work stands, and what to do next on your own machine.
Last updated 2026-09-12, after validating the port against the original MATLAB code.

---

## Pick it up

```bash
git clone https://github.com/PaleoLipidRR/BAYSPAR
cd BAYSPAR
git checkout claude/cool-hopper-11jy6s     # or master, if this has been merged

cd BAYSPARpy
python -m venv .venv && source .venv/bin/activate     # or conda, as you prefer
pip install -e ".[dev]"
python -m pytest                                       # expect 83 passed, 2 skipped
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

## Validation against the reference: done, with one caveat

This was the open job at handoff. It is closed — **the port is now checked against the original
MATLAB code actually executing**, not only against itself.

GNU Octave 8.4 (with the `statistics` package) runs the `.m` sources, so the reference values were
generated without a MATLAB licence:

```bash
cd ~/BAYSPAR
octave --no-gui --quiet --eval "pkg load statistics; addpath('docs/BAYSPARpy/audit'); generate_golden"
octave --no-gui --quiet --eval "pkg load statistics; addpath('docs/BAYSPARpy/audit'); verify_original"
cd BAYSPARpy && python -m pytest        # 83 passed, 2 skipped
```

Under MATLAB, the same two scripts, from the repository root:

```bash
matlab -batch "addpath('docs/BAYSPARpy/audit'); generate_golden"
matlab -batch "addpath('docs/BAYSPARpy/audit'); verify_original"
```

(The `>>` in an earlier draft of this note was the MATLAB prompt — those lines are typed *inside*
MATLAB, not in a shell. Sorry for the confusion. Both forms are now in each script's header.)

Two reference files are committed, so the tests run anywhere:

| File | Made by | Feeds |
|---|---|---|
| `audit/golden_matlab.mat` | `generate_golden.m` — recomputes the intermediates inline | `tests/test_golden.py`, 11 tests |
| `audit/original_matlab.mat` | `verify_original.m` — calls `bayspar_tex.m`, `bayspar_tex_analog.m` and `TEX_forward.m` unmodified | `tests/test_vs_original.py`, 7 tests |

The second matters more: the first checks the arithmetic against expressions I transcribed, so it
could not catch a misreading of the functions; the second runs the functions themselves. It confirms
the prior mean and grid cell exactly, the analogue set exactly, and the percentiles distributionally
(one reference run against the mean of 25 seeded Python runs, judged at 5σ per cell with the
family-wise rate in mind).

**The sharpest result:** `np.percentile(matlab_ensemble, [5, 50, 95], method="hazen")` returns
MATLAB's own saved `Preds` **exactly — all 579 values, difference 0.0**, where NumPy's default
method does not. `PORTING.md` BT-11 is now proven on real data rather than on `1:10`.

**The caveat: this was Octave, not MATLAB.** Octave is a re-implementation; its `prctile` comes from
the statistics package and its `sort` and `randn` are its own. Almost everything checked is plain
arithmetic where the two agree, and the exact `prctile` reproduction above is strong evidence for
the one convention that could plausibly differ — but a run under MATLAB proper would settle it. If
you have a licence handy, re-run the two commands above; the committed `.mat` files will be
overwritten and the same 18 tests will either still pass or tell you precisely which quantity moved.

**Found while doing this, in the reference rather than the port:** `bayspar_tex.m` line 66 is
`if runname=="SST"`. MATLAB reads `"SST"` as a string and compares it whole; Octave reads it as a
char array and compares elementwise, so a four-character runname (`'subT'`) raises *nonconformant
arguments*. The consequence is only that `verify_original.m` exercises the standard-mode path as
`'SST'` under Octave. `bayspar_tex_analog.m` and `TEX_forward.m` run either way.

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

1. ~~Golden files~~ — done, see above. Optionally re-run under MATLAB proper to remove the Octave
   caveat.
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
