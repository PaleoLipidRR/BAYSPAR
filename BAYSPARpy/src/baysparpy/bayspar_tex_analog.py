"""Analogue-mode prediction — a port of bayspar_tex_analog.m.

For deep time, where the site's modern location is not informative: the
regression parameters are pooled across calibration cells whose mean modern
TEX86 is close to the input series' mean.
"""
from __future__ import annotations

import numpy as np

from .bayspar_tex import solve_posterior
from .constants import DEFAULT_N_DRAWS
from .errors import SearchToleranceError
from .results import Prediction
from .stores import analog_cell_rows, get_coretops, get_draws
from .utils import as_column, check_runname, matlab_thinning, prctile


def find_analogs(cell_means: np.ndarray, centre: float, search_tol: float) -> np.ndarray:
    """Cells whose mean falls within a tolerance of a value (PORTING.md BTA-03).

    Args:
        cell_means: (80,) mean per calibration cell.
        centre: Value to search around — the mean of the input series.
        search_tol: Half-width of the window, inclusive at both ends.

    Returns:
        Indices of the selected cells.

    Raises:
        SearchToleranceError: if the window catches nothing.
    """
    sel = np.flatnonzero((cell_means >= centre - search_tol)
                         & (cell_means <= centre + search_tol))
    if sel.size == 0:
        raise SearchToleranceError()
    return sel


def bayspar_tex_analog(dats, prior_mean: float, prior_std: float, search_tol: float,
                       runname: str = "SST", *, n_draws: int = DEFAULT_N_DRAWS,
                       save_ensemble: bool = False, mode: str = "modern",
                       seed=None) -> Prediction:
    """Predict temperature from TEX86 by analogy, without using the site location.

    Args:
        dats: (N,) TEX86 values. Analogue selection uses their mean.
        prior_mean: Prior mean temperature, °C.
        prior_std: Prior standard deviation, °C. Do not make this narrow.
        search_tol: Tolerance in TEX86 units, compared against each calibration
            cell's mean TEX86. (In :func:`~baysparpy.tex_forward.tex_forward` the
            same argument is in °C — see PORTING.md TXF-04.)
        runname: ``"SST"`` or ``"subT"``.
        n_draws: Posterior draws per analogue location.
        save_ensemble: Keep the (N, n_analogues, n_draws) ensemble.
        mode: ``"modern"`` pairs τ² with the α and β of the same draw;
            ``"matlab"`` reproduces the reference implementation's mis-pairing
            and its flattened 2-D ensemble (PORTING.md BTA-05, BTA-08).
        seed: Seed for the draw.

    Returns:
        A :class:`~baysparpy.results.Prediction` carrying `analog_locs`.

    Raises:
        SearchToleranceError: if no cell falls within `search_tol`.
        ValueError: if `mode` is not one of the two.
    """
    if mode not in ("modern", "matlab"):
        raise ValueError(f"mode must be 'modern' or 'matlab', got {mode!r}")
    runname = check_runname(runname)
    tex = as_column(dats, "dats")
    rng = np.random.default_rng(seed)

    draws = get_draws(runname)
    coretops = get_coretops(runname)
    ind = matlab_thinning(draws.n_draws, n_draws)
    sel = find_analogs(coretops.cell_means("tex86"), float(tex.mean()), search_tol)

    rows = analog_cell_rows(draws, coretops)[sel]      # PORTING.md STO-01
    alpha = draws.alpha[np.ix_(rows, ind)]             # (n_an, M)
    beta = draws.beta[np.ix_(rows, ind)]
    tau2 = draws.tau2[ind]                             # (M,)
    n_an = len(sel)

    if mode == "matlab":
        # Column-major flattening, and tau2 tiled by draw -- so at flat position k
        # alpha/beta come from draw k//n_an and tau2 from draw k % M. See BTA-05.
        alpha_f = alpha.ravel(order="F")
        beta_f = beta.ravel(order="F")
        tau2_f = np.tile(tau2, n_an)
        post_mean, post_sd = solve_posterior(tex, alpha_f, beta_f, tau2_f,
                                             prior_mean, prior_std)
        preds = post_mean + rng.standard_normal(post_mean.shape) * post_sd
        pooled = preds
        ensemble = preds
    else:
        # Same draw index for alpha, beta and tau2 at every location.
        pinv_cov = prior_std ** -2.0
        den = pinv_cov + beta ** 2 / tau2                                  # (n_an, M)
        num = pinv_cov * prior_mean + (beta / tau2) * (tex[:, None, None] - alpha)
        post_mean = num / den                                              # (N, n_an, M)
        post_sd = np.sqrt(1.0 / den)
        preds = post_mean + rng.standard_normal(post_mean.shape) * post_sd
        pooled = preds.reshape(len(tex), -1)   # percentiles pool over locations too
        ensemble = preds

    return Prediction(
        preds=prctile(pooled, [5, 50, 95], axis=1),
        prior_mean=prior_mean,
        prior_std=prior_std,
        analog_locs=coretops.locs[sel],
        ensemble=ensemble if save_ensemble else None,
        metadata={"runname": runname, "n_draws": n_draws, "mode": mode, "seed": seed,
                  "thinning_index": ind, "analog_cells": sel, "search_tol": search_tol},
    )
