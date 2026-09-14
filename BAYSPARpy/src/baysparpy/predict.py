"""The modern API.

Function and argument names follow ``brews/baysparpy`` so a notebook ports by
changing ``from bayspar import ...`` to ``from baysparpy import ...``. Two
defaults differ, deliberately: ``n_draws`` is MATLAB's 1000 rather than 5000
(SPEC §14), and the argument is named ``n_draws`` rather than ``nens``.
"""
from __future__ import annotations

from .bayspar_tex import bayspar_tex
from .bayspar_tex_analog import bayspar_tex_analog
from .constants import DEFAULT_N_DRAWS
from .results import ForwardPrediction, Prediction
from .tex_forward import tex_forward


def predict_seatemp(tex, lat: float, lon: float, prior_std: float,
                    temptype: str = "sst", prior_mean: float | None = None, *,
                    n_draws: int = DEFAULT_N_DRAWS, seed=None,
                    save_ensemble: bool = False, **kwargs) -> Prediction:
    """Predict sea temperature from TEX86 at a known location.

    Args:
        tex: (N,) TEX86 observations from one location.
        lat: Latitude, -90 to 90.
        lon: Longitude, -180 to 180.
        prior_std: Prior standard deviation, °C.
        temptype: ``"sst"`` or ``"subt"`` (``"SST"``/``"subT"`` also accepted).
        prior_mean: Prior mean, °C. Searched from the climatology when omitted.
        n_draws: Posterior draws to mix over.
        seed: Seed for the draw.
        save_ensemble: Keep the full ensemble on the result.
        **kwargs: Passed to :func:`~baysparpy.bayspar_tex.bayspar_tex`.

    Returns:
        A :class:`~baysparpy.results.Prediction`.
    """
    return bayspar_tex(tex, lon, lat, prior_std, temptype, n_draws=n_draws,
                       seed=seed, save_ensemble=save_ensemble,
                       prior_mean=prior_mean, **kwargs)


def predict_seatemp_analog(tex, prior_std: float, search_tol: float,
                           temptype: str = "sst", prior_mean: float | None = None, *,
                           n_draws: int = DEFAULT_N_DRAWS, seed=None,
                           save_ensemble: bool = False, **kwargs) -> Prediction:
    """Predict sea temperature from TEX86 by analogy, ignoring the site location.

    Args:
        tex: (N,) TEX86 values; their mean drives analogue selection.
        prior_std: Prior standard deviation, °C.
        search_tol: Tolerance in TEX86 units.
        temptype: ``"sst"`` or ``"subt"``.
        prior_mean: Prior mean, °C. Required — analogue mode cannot search for it.
        n_draws: Posterior draws per analogue location.
        seed: Seed for the draw.
        save_ensemble: Keep the (N, n_analogues, n_draws) ensemble.
        **kwargs: Passed to
            :func:`~baysparpy.bayspar_tex_analog.bayspar_tex_analog`.

    Returns:
        A :class:`~baysparpy.results.Prediction` carrying `analog_locs`.

    Raises:
        ValueError: if `prior_mean` is not given.
    """
    if prior_mean is None:
        raise ValueError(
            "analogue mode has no site location to search from, so prior_mean is "
            "required (the Wilson Lake demo uses 30 °C)"
        )
    return bayspar_tex_analog(tex, prior_mean, prior_std, search_tol, temptype,
                              n_draws=n_draws, seed=seed,
                              save_ensemble=save_ensemble, **kwargs)


def predict_tex(seatemp, lat, lon, temptype: str = "sst", *, seed=None,
                **kwargs) -> ForwardPrediction:
    """Forward-model TEX86 from temperature.

    Args:
        seatemp: (N,) temperatures in °C.
        lat: Latitude, or one per temperature.
        lon: Longitude, or one per temperature.
        temptype: ``"sst"`` or ``"subt"``.
        seed: Seed for the draw.
        **kwargs: Passed to :func:`~baysparpy.tex_forward.tex_forward`
            (``type``, ``search_tol``, ``n_out``, ...).

    Returns:
        A :class:`~baysparpy.results.ForwardPrediction`.
    """
    return tex_forward(lat, lon, seatemp, temptype, seed=seed, **kwargs)
