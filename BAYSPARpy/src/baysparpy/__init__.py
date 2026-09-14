"""BAYSPAR — a Bayesian, spatially-varying calibration for the TEX86 proxy.

A Python recode of the MATLAB package of Tierney & Tingley (2014, 2015). The
translation is documented line by line in ``docs/BAYSPARpy/PORTING.md``; each
deliberate difference from MATLAB carries an ID quoted in the relevant docstring.

Cite the original methodology:

    Tierney, J. E., & Tingley, M. P. (2014). A Bayesian, spatially-varying
    calibration model for the TEX86 proxy. Geochimica et Cosmochimica Acta, 127,
    83-106. https://doi.org/10.1016/j.gca.2013.11.026
"""
from importlib.metadata import PackageNotFoundError as _NotFound, version as _version

from .bayspar_tex import bayspar_tex, solve_posterior
from .bayspar_tex_analog import bayspar_tex_analog, find_analogs
from .constants import (
    DEFAULT_N_DRAWS, EARTH_RADIUS_KM, GRID_HALF_SPACE, MAX_DIST_KM, MIN_NUM, RUNNAMES,
)
from .distance import earth_chord_distances
from .errors import (
    AmbiguousGridCellError, BaysparError, EnsembleSizeError, SearchToleranceError,
)
from .predict import predict_seatemp, predict_seatemp_analog, predict_tex
from .results import ForwardPrediction, Prediction
from .stores import (
    CoreTops, Draws, SeaTempObs, find_modeloutput, get_coretops, get_draws, get_seatemp,
)
from .tex_forward import tex_forward
from .utils import matlab_thinning, prctile

try:
    __version__ = _version("bayspar")
except _NotFound:      # running from a source checkout that was never installed
    __version__ = "0.1.0.dev0"

__all__ = [
    "__version__",
    # MATLAB-faithful layer
    "bayspar_tex", "bayspar_tex_analog", "tex_forward", "earth_chord_distances",
    # modern layer
    "predict_seatemp", "predict_seatemp_analog", "predict_tex",
    # results
    "Prediction", "ForwardPrediction",
    # stores
    "Draws", "SeaTempObs", "CoreTops",
    "get_draws", "get_seatemp", "get_coretops", "find_modeloutput",
    # pieces worth reusing
    "solve_posterior", "find_analogs", "matlab_thinning", "prctile",
    # errors
    "BaysparError", "SearchToleranceError", "EnsembleSizeError",
    "AmbiguousGridCellError",
    # constants
    "EARTH_RADIUS_KM", "GRID_HALF_SPACE", "MAX_DIST_KM", "MIN_NUM",
    "DEFAULT_N_DRAWS", "RUNNAMES",
]
