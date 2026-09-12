"""Standard-mode prediction of temperature from TEX86 — a port of bayspar_tex.m."""
from __future__ import annotations

import numpy as np

from .constants import DEFAULT_N_DRAWS, MAX_DIST_KM, MIN_NUM
from .results import Prediction
from .stores import get_draws, get_seatemp
from .utils import as_column, check_runname, matlab_thinning, prctile


def solve_posterior(tex: np.ndarray, alpha: np.ndarray, beta: np.ndarray,
                    tau2: np.ndarray, prior_mean: float,
                    prior_std: float) -> tuple[np.ndarray, np.ndarray]:
    """Conjugate normal update, one posterior per observation and draw.

    A port of ``bayspar_tex.m`` lines 144-152 (PORTING.md BT-09). MATLAB tiles
    every term to (N, M); the prior precision and the denominator do not depend
    on the observations, so this broadcasts instead — which is what removes the
    O(N³) per-draw linear algebra other ports carry (SPEC §5).

    Args:
        tex: (N,) TEX86 observations.
        alpha: (M,) intercept draws.
        beta: (M,) slope draws.
        tau2: (M,) residual variance draws, paired with `alpha` and `beta`.
        prior_mean: Prior mean temperature, °C.
        prior_std: Prior standard deviation, °C.

    Returns:
        ``(post_mean, post_sd)`` — (N, M) and (M,).
    """
    pinv_cov = prior_std ** -2.0
    den = pinv_cov + beta ** 2 / tau2                              # (M,)
    num = pinv_cov * prior_mean + (beta / tau2) * (tex[:, None] - alpha)   # (N, M)
    return num / den, np.sqrt(1.0 / den)


def bayspar_tex(dats, lon: float, lat: float, prior_std: float,
                runname: str = "SST", *, n_draws: int = DEFAULT_N_DRAWS,
                save_ensemble: bool = False, seed=None, strict_grid: bool = False,
                max_dist: float = MAX_DIST_KM, min_num: int = MIN_NUM,
                prior_mean: float | None = None) -> Prediction:
    """Predict temperature from a TEX86 series at a known location.

    Standard mode: the regression parameters come from the calibration cell that
    contains the site, so it assumes oceanographic conditions like today's. For
    deep time use :func:`~baysparpy.bayspar_tex_analog.bayspar_tex_analog`.

    Args:
        dats: (N,) TEX86 observations from one location.
        lon: Longitude, -180 to 180.
        lat: Latitude, -90 to 90.
        prior_std: Prior standard deviation on temperature, °C. Do not make this
            narrow; 6-10 is a reasonable start.
        runname: ``"SST"`` or ``"subT"``.
        n_draws: Posterior draws to mix over, spread across the chain. The store
            holds 20,000 (not the 15,000 the MATLAB docs claim).
        save_ensemble: Keep the full (N, n_draws) ensemble on the result.
        seed: Seed for the draw. ``None`` matches MATLAB's non-reproducibility.
        strict_grid: Raise instead of choosing when the site is exactly on a cell
            boundary (PORTING.md BT-08).
        max_dist: Prior-mean search radius, km. MATLAB hard-codes 500.
        min_num: Observations to fall back on when nothing is within `max_dist`.
        prior_mean: Override the searched prior mean, °C. MATLAB always searches.

    Returns:
        A :class:`~baysparpy.results.Prediction`.
    """
    runname = check_runname(runname)
    tex = as_column(dats, "dats")
    rng = np.random.default_rng(seed)

    draws = get_draws(runname)
    ind = matlab_thinning(draws.n_draws, n_draws)
    cell = draws.cell_index(lon, lat, strict=strict_grid)
    alpha = draws.alpha[cell, ind]
    beta = draws.beta[cell, ind]
    tau2 = draws.tau2[ind]

    n_below = None
    if prior_mean is None:
        prior_mean, n_below = get_seatemp(runname).prior_mean(
            lon, lat, max_dist=max_dist, min_num=min_num)

    post_mean, post_sd = solve_posterior(tex, alpha, beta, tau2, prior_mean, prior_std)
    preds = post_mean + rng.standard_normal(post_mean.shape) * post_sd

    return Prediction(
        preds=prctile(preds, [5, 50, 95], axis=1),
        prior_mean=prior_mean,
        prior_std=prior_std,
        site_loc=(float(lon), float(lat)),
        grid_loc=tuple(draws.locs[cell]),
        ensemble=preds if save_ensemble else None,
        metadata={"runname": runname, "n_draws": n_draws, "mode": "standard",
                  "seed": seed, "thinning_index": ind, "grid_cell": cell,
                  "n_obs_within_max_dist": n_below},
    )
