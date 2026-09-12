"""Forward model: TEX86 from temperature — a port of TEX_forward.m."""
from __future__ import annotations

import numpy as np

from .bayspar_tex_analog import find_analogs
from .constants import FORWARD_N_OUT
from .results import ForwardPrediction
from .stores import analog_cell_rows, get_coretops, get_draws
from .utils import as_column, check_runname, matlab_thinning


def tex_forward(lat, lon, t, runname: str = "SST", type: str = "standard",
                search_tol: float | None = None, *, n_draws: int | None = None,
                n_out: int = FORWARD_N_OUT, seed=None, mode: str = "modern",
                strict_grid: bool = False) -> ForwardPrediction:
    """Forward-model TEX86 from temperature.

    Note the argument order is ``(lat, lon, t)`` — latitude first, the opposite of
    :func:`~baysparpy.bayspar_tex.bayspar_tex`. That is MATLAB's order, kept.

    With several locations, temperatures pair with them **elementwise**: this is
    one temperature per site, not a series at each of several sites
    (PORTING.md TXF-03).

    Args:
        lat: Latitude, or one per temperature.
        lon: Longitude, or one per temperature.
        t: (N,) temperatures in °C.
        runname: ``"SST"`` or ``"subT"``.
        type: ``"standard"`` or ``"analog"``. Named as in MATLAB, despite
            shadowing the builtin.
        search_tol: Required for analogue mode, **in °C** — the comparison is
            against each cell's mean target temperature, not its TEX86.
        n_draws: Thin the parameters first. MATLAB uses all 20,000 (the default).
        n_out: Columns to return. MATLAB hard-codes 1000 (PORTING.md TXF-07).
        seed: Seed for the draw and the subsample.
        mode: ``"matlab"`` reproduces the τ² mis-pairing in analogue mode.
        strict_grid: Raise on a site exactly on a cell boundary.

    Returns:
        A :class:`~baysparpy.results.ForwardPrediction`, TEX86 clipped to [0, 1].

    Raises:
        ValueError: on an unknown `type`, a missing `search_tol`, or a
            location/temperature length mismatch.
    """
    if type not in ("standard", "analog"):
        raise ValueError('please enter "analog" to specify analog mode')
    if type == "analog" and search_tol is None:
        raise ValueError("To use analog mode, enter a search tolerance in TEX units")
    runname = check_runname(runname)
    t = as_column(t, "t")
    lat = as_column(lat, "lat")
    lon = as_column(lon, "lon")
    rng = np.random.default_rng(seed)

    draws = get_draws(runname)
    ind = (matlab_thinning(draws.n_draws, n_draws) if n_draws is not None
           else np.arange(draws.n_draws))
    tau2 = draws.tau2[ind]

    if type == "standard":
        if lat.size != lon.size:
            raise ValueError(f"lat has {lat.size} values and lon has {lon.size}")
        cells = [draws.cell_index(lo, la, strict=strict_grid) for lo, la in zip(lon, lat)]
        alpha = draws.alpha[np.ix_(cells, ind)]        # (n_loc, M)
        beta = draws.beta[np.ix_(cells, ind)]
        if len(cells) == 1:
            mu = t[:, None] * beta[0] + alpha[0]       # every t at the one site
            tau2_b = tau2
        elif len(cells) == t.size:
            mu = t[:, None] * beta + alpha             # one t per site, elementwise
            tau2_b = tau2
        else:
            raise ValueError(
                f"{len(cells)} locations and {t.size} temperatures: tex_forward pairs "
                "them elementwise, so pass one temperature per location (or a single "
                "location for a whole series)"
            )
    else:
        coretops = get_coretops(runname)
        # Analogue search here is in TEMPERATURE units, against Target_Stack --
        # bayspar_tex_analog searches Obs_Stack, in TEX86 units. PORTING.md TXF-04.
        sel = find_analogs(coretops.cell_means("target_t"), float(t.mean()), search_tol)
        rows = analog_cell_rows(draws, coretops)[sel]
        alpha = draws.alpha[np.ix_(rows, ind)]
        beta = draws.beta[np.ix_(rows, ind)]
        if mode == "matlab":
            alpha_f, beta_f = alpha.ravel(order="F"), beta.ravel(order="F")
            tau2_b = np.tile(tau2, len(sel))
        else:
            alpha_f, beta_f = alpha.ravel(), beta.ravel()
            tau2_b = np.tile(tau2, len(sel))   # ravel() is draw-fastest, so tile matches
        mu = t[:, None] * beta_f + alpha_f

    tex = rng.normal(mu, np.sqrt(tau2_b))
    n_total = tex.shape[1]
    if n_out > n_total:
        raise ValueError(f"n_out={n_out} exceeds the {n_total} available draws")
    tex = tex[:, rng.choice(n_total, size=n_out, replace=False)]

    clipped = float(np.mean((tex > 1) | (tex < 0)))
    np.clip(tex, 0.0, 1.0, out=tex)
    return ForwardPrediction(
        tex86=tex,
        metadata={"runname": runname, "type": type, "mode": mode, "seed": seed,
                  "n_draws": ind.size, "n_out": n_out, "clipped_fraction": clipped},
    )
