"""Small translations that several functions share.

Every one of these is a place where the obvious Python differs from MATLAB;
see the matching entry in ../../../docs/BAYSPARpy/PORTING.md.
"""
from __future__ import annotations

import numpy as np

from .constants import RUNNAME_ALIASES, RUNNAMES
from .errors import EnsembleSizeError


def as_column(x, name: str = "input") -> np.ndarray:
    """MATLAB ``x = x(:)`` — flatten to a 1-D float array (PORTING.md GEN-02).

    MATLAB flattens column-major and NumPy row-major, so anything with more than
    one non-trivial dimension is rejected rather than flattened differently.

    Args:
        x: Scalar or array-like.
        name: Name to use in the error message.

    Returns:
        A 1-D float64 array.

    Raises:
        ValueError: if `x` has more than one non-trivial dimension.
    """
    arr = np.asarray(x, dtype=float)
    if arr.ndim > 1 and sum(d > 1 for d in arr.shape) > 1:
        raise ValueError(
            f"{name} has shape {arr.shape}; MATLAB's (:) flattens column-major and "
            "NumPy row-major, so pass a 1-D sequence to make the order explicit"
        )
    return arr.reshape(-1)


def matlab_thinning(n_total: int, n_draws: int) -> np.ndarray:
    """MATLAB ``round(linspace(1, n_total, n_draws))``, 0-based (PORTING.md BT-06).

    MATLAB's ``round`` is half-away-from-zero; :func:`numpy.round` is half-to-even,
    which picks a different draw wherever the linear spacing lands on .5.

    Args:
        n_total: Draws in the store.
        n_draws: Draws wanted, spread across the whole chain.

    Returns:
        0-based indices, length `n_draws`.

    Raises:
        EnsembleSizeError: if more draws are requested than exist.
    """
    if n_draws > n_total:
        raise EnsembleSizeError(n_total, n_draws)
    if n_draws < 1:
        raise ValueError(f"n_draws must be at least 1, got {n_draws}")
    return np.floor(np.linspace(1, n_total, n_draws) + 0.5).astype(int) - 1


def prctile(a: np.ndarray, q, axis: int = -1) -> np.ndarray:
    """MATLAB ``prctile`` (PORTING.md BT-11).

    MATLAB uses Hazen plotting positions ``(i - 0.5) / n``; NumPy's default
    (linear, type 7) does not match and biases the tails toward the middle.

    Args:
        a: Sample array.
        q: Percentiles in 0-100.
        axis: Axis to reduce.

    Returns:
        Percentiles, with the reduced axis moved to the end.
    """
    out = np.percentile(a, q, axis=axis, method="hazen")
    return np.moveaxis(out, 0, -1)


def check_runname(runname: str) -> str:
    """Validate `runname` up front (PORTING.md BT-03).

    MATLAB falls through both branches on an unrecognised name and fails later
    with an unrelated message. Lower-case spellings are accepted as aliases of
    brews/baysparpy's ``temptype``.

    Args:
        runname: ``"SST"`` or ``"subT"`` (``"sst"``/``"subt"`` also accepted).

    Returns:
        The canonical spelling.

    Raises:
        ValueError: if the name is not one of the two fitted models.
    """
    try:
        return RUNNAME_ALIASES[str(runname)]
    except KeyError:
        raise ValueError(
            f"runname must be one of {RUNNAMES!r} (or their lower-case aliases), "
            f"got {runname!r}"
        ) from None
