# BAYSPARpy site map

How the modules in `src/baysparpy/` depend on each other, which MATLAB file each
one ports, what data they read, and which tests exercise them.

The line-by-line translation record lives in `docs/BAYSPARpy/PORTING.md` at the
root of the BAYSPAR repository. IDs such as `BT-06` in the docstrings refer to it.

## Modules at a glance

| Module | Layer | Ports | Role |
|---|---|---|---|
| `__init__.py` | package | — | Re-exports the public API |
| `predict.py` | modern API | — | `predict_seatemp`, `predict_seatemp_analog`, `predict_tex`: thin, keyword-friendly wrappers |
| `bayspar_tex.py` | MATLAB-faithful | `bayspar_tex.m` | Standard mode: TEX86 → temperature at a known site. Also holds `solve_posterior` |
| `bayspar_tex_analog.py` | MATLAB-faithful | `bayspar_tex_analog.m` | Analogue mode: TEX86 → temperature, ignoring location. Also holds `find_analogs` |
| `tex_forward.py` | MATLAB-faithful | `TEX_forward.m` | Forward model: temperature → TEX86 |
| `distance.py` | MATLAB-faithful | `EarthChordDistances_2.m` | `earth_chord_distances` |
| `stores.py` | data | the `load` calls | Finds `ModelOutput/` and reads the `.mat` files into `Draws`, `SeaTempObs`, `CoreTops` (cached) |
| `results.py` | output | `Output_Struct` | `Prediction`, `ForwardPrediction` dataclasses |
| `utils.py` | shared | small MATLAB idioms | `as_column`, `matlab_thinning`, `prctile`, `check_runname` |
| `constants.py` | shared | hard-coded numbers | `EARTH_RADIUS_KM`, `GRID_HALF_SPACE`, `MAX_DIST_KM`, `MIN_NUM`, `DEFAULT_N_DRAWS`, `RUNNAMES`, … |
| `errors.py` | shared | MATLAB `error(...)` | `BaysparError` and its three subclasses |

## Import graph

An arrow `A --> B` means *A imports from B*. Modules higher up depend on those
below; nothing points upward, so there are no import cycles.

```mermaid
flowchart TD
    init["__init__.py<br/>(public API)"]

    subgraph modern["Modern API"]
        predict["predict.py"]
    end

    subgraph faithful["MATLAB-faithful layer"]
        tex["bayspar_tex.py"]
        analog["bayspar_tex_analog.py"]
        fwd["tex_forward.py"]
    end

    subgraph data["Data & output"]
        stores["stores.py"]
        results["results.py"]
        distance["distance.py"]
    end

    subgraph shared["Shared"]
        utils["utils.py"]
        constants["constants.py"]
        errors["errors.py"]
    end

    init --> predict
    init -.->|re-exports| tex & analog & fwd & stores & results & distance & utils & constants & errors

    predict --> tex
    predict --> analog
    predict --> fwd
    predict --> results
    predict --> constants

    analog -->|solve_posterior| tex
    fwd -->|find_analogs| analog

    tex --> stores
    tex --> results
    tex --> utils
    tex --> constants

    analog --> stores
    analog --> results
    analog --> utils
    analog --> errors
    analog --> constants

    fwd --> stores
    fwd --> results
    fwd --> utils
    fwd --> constants

    stores --> distance
    stores --> utils
    stores --> errors
    stores --> constants

    results -.->|prctile, lazily| utils
    distance --> constants
    utils --> constants
    utils --> errors
```

Two links cross between the three ported functions, and both reuse a shared step
instead of copying it:

- `bayspar_tex_analog` calls `bayspar_tex.solve_posterior` for the conjugate normal update.
- `tex_forward` calls `bayspar_tex_analog.find_analogs` for its analogue cell search.

## Call flow for each entry point

### Standard mode: `predict_seatemp` / `bayspar_tex`

```mermaid
flowchart LR
    A["predict_seatemp"] --> B["bayspar_tex"]
    B --> C["check_runname<br/>as_column"]
    B --> D["get_draws(runname)"]
    D --> D2["Draws.cell_index(lon, lat)"]
    D2 --> E1["earth_chord_distances"]
    B --> F["matlab_thinning"]
    B --> G["get_seatemp(runname)<br/>.prior_mean(...)"]
    G --> E2["earth_chord_distances"]
    B --> H["solve_posterior"]
    B --> I["prctile 5/50/95"]
    I --> J["Prediction"]
```

### Analogue mode: `predict_seatemp_analog` / `bayspar_tex_analog`

```mermaid
flowchart LR
    A["predict_seatemp_analog"] --> B["bayspar_tex_analog"]
    B --> C["get_draws + get_coretops"]
    B --> F["matlab_thinning"]
    B --> K["find_analogs<br/>(on CoreTops.cell_means('tex86'))"]
    K -->|none found| X["SearchToleranceError"]
    K --> L["analog_cell_rows"]
    L --> H["solve_posterior<br/>(from bayspar_tex.py)"]
    H --> I["prctile 5/50/95"]
    I --> J["Prediction"]
```

### Forward model: `predict_tex` / `tex_forward`

```mermaid
flowchart LR
    A["predict_tex"] --> B["tex_forward"]
    B --> C["get_draws(runname)"]
    B --> F["matlab_thinning"]
    B -->|standard| D["Draws.cell_index per site"]
    B -->|analogue| K["get_coretops<br/>find_analogs (from bayspar_tex_analog.py)<br/>analog_cell_rows"]
    D --> J["ForwardPrediction"]
    K --> J
```

## Data files read by `stores.py`

`find_modeloutput()` uses `$BAYSPAR_MODELOUTPUT` if it is set. Otherwise it
walks up from the package to the first `ModelOutput/` that contains
`obsSST.mat`. Every file is loaded once and cached (`lru_cache`).

```mermaid
flowchart LR
    subgraph MO["ModelOutput/"]
        P["Output_SpatAg_{SST,subT}/<br/>params_standard.mat"]
        O["obsSST.mat / obssubT.mat"]
        DI["Data_Input_SpatAg_{SST,subT}.mat"]
    end
    P --> GD["get_draws → Draws"]
    O --> GS["get_seatemp → SeaTempObs"]
    DI --> GC["get_coretops → CoreTops"]
```

`params_analog.mat` is never read. It is bit-identical to a subset of rows in
`params_standard.mat`, and `analog_cell_rows` selects those rows (PORTING.md STO-01).

## Tests

```mermaid
flowchart LR
    conf["conftest.py<br/>fixtures: modeloutput, lopes, wilsonlake"] --> stores["stores.find_modeloutput"]

    tp["test_port.py<br/>one test per PORTING.md ID"] --> pkg["baysparpy<br/>(all modules)"]
    tg["test_golden.py"] --> pkg
    tg --> gm[("audit/golden_matlab.mat")]
    tvo["test_vs_original.py"] --> pkg
    tvo --> om[("audit/original_matlab.mat")]
    trt["test_reference_trace.py"] --> pkg
    trt --> rt[("audit/reference_trace.txt")]
    tnb["test_no_blas.py"] --> pkg
    tpd["test_porting_doc.py"] --> pd[("docs/BAYSPARpy/PORTING.md")]
    tpd --> tp
```

| Test file | Checks |
|---|---|
| `test_port.py` | Each deliberate MATLAB → Python decision, one test per PORTING.md ID |
| `test_golden.py` | Output matches values generated in MATLAB by `audit/generate_golden.m` |
| `test_vs_original.py` | Output matches the original, unmodified MATLAB functions (`audit/verify_original.m`) |
| `test_reference_trace.py` | Intermediate values match `audit/reference_trace.txt` |
| `test_no_blas.py` | The closed-form path never calls a threaded BLAS |
| `test_porting_doc.py` | Every PORTING.md ID has a test and every test ID exists in PORTING.md |

The `audit/` paths are under `docs/BAYSPARpy/` in the BAYSPAR repository.

## Notebook

`notebooks/compare_implementations.ipynb` imports `baysparpy as ours` and,
when it is installed, `bayspar as brews`. It runs both on Jess Tierney's two
demos and compares them against the MATLAB results in
`audit/original_matlab.mat`.
