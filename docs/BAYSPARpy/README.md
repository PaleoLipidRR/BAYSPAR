# BAYSPARpy planning materials

Design work for a Python recode of this MATLAB package, to live in a new repository
(`PaleoLipidRR/BAYSPARpy`). Nothing here changes the MATLAB code.

| File | What it is |
|---|---|
| [`SPEC.md`](SPEC.md) | The specification: scope, naming, audit of the MATLAB package, port design, CmdStan/Stan plan, validation plan, milestones, open decisions. |
| [`audit/audit_modeloutput.py`](audit/audit_modeloutput.py) | Verifies every factual claim SPEC §3 and §8 make about `ModelOutput/`. Exits non-zero if one stops holding. |
| [`audit/bench_closed_form.py`](audit/bench_closed_form.py) | Prototype of the closed-form standard prediction, benchmarked against the dense per-draw solve used by `brews/baysparpy`. Backs SPEC §5. |

Both scripts need only `numpy` and `scipy`, and run from the repository root:

```bash
python docs/BAYSPARpy/audit/audit_modeloutput.py
python docs/BAYSPARpy/audit/bench_closed_form.py
```

Open questions needing Jess (and Brewster, on the PyPI name) are collected in SPEC §14.
