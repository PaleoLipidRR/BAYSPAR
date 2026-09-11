# BAYSPARpy planning materials

Design and review material for a Python recode of this MATLAB package, to live in a new repository
(`PaleoLipidRR/BAYSPARpy`). Nothing here changes the MATLAB code.

| File | What it is | Who reads it |
|---|---|---|
| [`SPEC.md`](SPEC.md) | What is being built and why: scope, naming, audit of the MATLAB package, package layout and API, data packaging, the CmdStan and Stan plan, validation plan, milestones, open decisions. | Anyone deciding whether this is the right project. |
| [`PORTING.md`](PORTING.md) | **The translation record.** 45 numbered entries walking every line of the four MATLAB functions and stating what the Python does instead, with a register of the nine places it does not behave like MATLAB. | Anyone reviewing whether the recode is faithful. |
| [`audit/reference_trace.txt`](audit/reference_trace.txt) | Deterministic values — thinning indices, distances, prior means, grid cells, analogue selections, per-draw posteriors, percentile conventions — to compare against a MATLAB session. | Whoever has MATLAB open. |
| [`audit/audit_modeloutput.py`](audit/audit_modeloutput.py) | Verifies every factual claim SPEC §3 and §8 make about `ModelOutput/`; exits non-zero if one stops holding. | — |
| [`audit/bench_closed_form.py`](audit/bench_closed_form.py) | Prototype of the closed-form standard prediction, benchmarked against the dense per-draw solve used by `brews/baysparpy`. Backs SPEC §5. | — |
| [`audit/reference_trace.py`](audit/reference_trace.py) | Regenerates `reference_trace.txt`. | — |

All three scripts need only `numpy` and `scipy`, and run from the repository root:

```bash
python docs/BAYSPARpy/audit/audit_modeloutput.py
python docs/BAYSPARpy/audit/bench_closed_form.py
python docs/BAYSPARpy/audit/reference_trace.py > docs/BAYSPARpy/audit/reference_trace.txt
```

## How to review

Start with the register at the top of [`PORTING.md`](PORTING.md) — nine rows, each naming who
decides. Everything else in that document is either an exact translation or a proven equivalence,
and is there so the port can be checked line by line without reading Python.

Decisions still needing Jess (or Brewster, on the PyPI name) are collected in SPEC §14.
