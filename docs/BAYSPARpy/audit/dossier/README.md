# Comparison dossier

An HTML page that sets the original MATLAB functions, `brews/baysparpy` and BAYSPARpy side by side:
the math each function computes, every step of the code in three columns with real line numbers,
and the results of running all three on the demo inputs.

```bash
cd BAYSPARpy
uv sync                                                   # installs brews/baysparpy too
uv run python ../docs/BAYSPARpy/audit/dossier/gen_data.py # runs all three, writes data.json (~30 s)
uv run python ../docs/BAYSPARpy/audit/dossier/build.py    # writes dossier.html
```

| File | What it does |
|---|---|
| `gen_data.py` | Runs BAYSPARpy and brews/baysparpy on the inputs stored in `../original_matlab.mat`, reads the MATLAB (Octave) results from the same file, and saves every number the page charts. Includes a timing run, so those figures depend on the machine. |
| `template.html` | The page. Each `<!--stage ... -->` block names the file and line range to show for each implementation, its status, and a caption. |
| `build.py` | Pulls those line ranges from the source files, highlights them with Pygments, and inlines `data.json`. Fails if a range runs past the end of a file, so an edited source cannot silently shift the excerpts. |

`data.json` and `dossier.html` are build outputs and are not committed. Line ranges in
`template.html` are tied to the current sources; after editing a ported function, check its stages
still point at the right lines.
