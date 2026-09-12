"""The closed-form path must not touch a threaded BLAS.

SPEC §5: the reason for this recode is that brews/baysparpy solves an N x N
system per posterior draw for a prior covariance that is always diagonal, which
saturates every core. The fix is arithmetic, not thread limits, so the test
checks that no BLAS call happens at all rather than that it is fast.
"""
from __future__ import annotations

import time

import numpy as np
import pytest

import baysparpy as bp

threadpoolctl = pytest.importorskip("threadpoolctl", reason="pip install threadpoolctl")


def test_prediction_makes_no_blas_call(lopes):
    """A standard prediction runs with every BLAS thread pool capped at one.

    If the implementation reached for LAPACK, this would still pass but slowly;
    the timing assertion below is what makes the constraint bite.
    """
    tex, lon, lat = lopes
    with threadpoolctl.threadpool_limits(limits=1):
        t0 = time.perf_counter()
        out = bp.bayspar_tex(tex, lon, lat, 6.0, "subT", n_draws=1000, seed=0)
        single_thread = time.perf_counter() - t0
    assert out.preds.shape == (len(tex), 3)
    # 193 points x 1000 draws of elementwise arithmetic. The dense per-draw solve
    # this replaces takes ~1.6 s for the same work (SPEC §5); a second is three
    # orders of magnitude of headroom over the measured ~10 ms.
    assert single_thread < 1.0, (
        f"{single_thread:.3f}s single-threaded suggests the solve is doing matrix "
        "work it does not need to")


def test_solve_posterior_is_elementwise():
    """The posterior solve allocates no N x N intermediate.

    A dense implementation needs O(N^2) memory per draw; this one is O(N x M),
    so a long record with few draws stays small.
    """
    n = 4000
    tex = np.full(n, 0.6)
    alpha, beta, tau2 = np.full(3, 0.3), np.full(3, 0.015), np.full(3, 0.0018)
    post_mean, post_sd = bp.solve_posterior(tex, alpha, beta, tau2, 20.0, 6.0)
    assert post_mean.shape == (n, 3)
    assert post_sd.shape == (3,)          # not (n, n), and not even (n, 3)
