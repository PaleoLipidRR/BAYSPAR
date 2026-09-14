"""The package must reproduce audit/reference_trace.txt.

That file is what a reviewer compares against a MATLAB session, so if the code
drifts from it, either the code is wrong or the trace needs regenerating -- and
either way someone should look.
"""
from __future__ import annotations

import numpy as np
import pytest

import baysparpy as bp
from baysparpy import stores, utils

# Values as printed in docs/BAYSPARpy/audit/reference_trace.txt.
PRIOR_MEANS = {
    ("SST", "castaneda2010"): (21.7174001174, 22),
    ("SST", "lopes_santos2010"): (26.8674410752, 56),
    ("SST", "shevenell2011"): (0.0535030915, 97),
    ("subT", "castaneda2010"): (18.5711564752, 22),
    ("subT", "lopes_santos2010"): (19.1673884681, 56),
    ("subT", "shevenell2011"): (-0.4346589233, 93),
}
GRID_CELLS = {
    "castaneda2010": (30.0, 40.0),
    "lopes_santos2010": (-10.0, 0.0),
    "shevenell2011": (-70.0, -60.0),
}
SITES = {
    "castaneda2010": (34.0733, 31.6517),
    "lopes_santos2010": (-17.6635, 9.1660),
    "shevenell2011": (-64.2080, -64.8527),
}


@pytest.mark.parametrize("runname,series", sorted(PRIOR_MEANS))
def test_trace_C_prior_means(runname, series):
    """Section C: prior mean and the count of observations within 500 km."""
    expected_mean, expected_n = PRIOR_MEANS[(runname, series)]
    lon, lat = SITES[series]
    mean, n = stores.get_seatemp(runname).prior_mean(lon, lat, max_dist=500.0, min_num=1)
    assert n == expected_n
    assert mean == pytest.approx(expected_mean, abs=1e-9)


@pytest.mark.parametrize("series", sorted(GRID_CELLS))
def test_trace_C_grid_cells(series):
    """Section C: every demo site falls in exactly one cell."""
    lon, lat = SITES[series]
    draws = stores.get_draws("SST")
    cell = draws.cell_index(lon, lat, strict=True)      # strict: no boundary ambiguity
    assert tuple(draws.locs[cell]) == GRID_CELLS[series]


def test_trace_A_thinning():
    """Section A: the thinning index for five values of n_draws."""
    for n_draws, first, last in ((10, 1, 20000), (999, 1, 20000), (1000, 1, 20000),
                                 (1001, 1, 20000), (20000, 1, 20000)):
        ind = utils.matlab_thinning(20000, n_draws) + 1
        assert ind[0] == first and ind[-1] == last and ind.size == n_draws
    np.testing.assert_array_equal(utils.matlab_thinning(20000, 999)[:5] + 1,
                                  [1, 21, 41, 61, 81])


def test_trace_B_distances():
    """Section B: chordal distances, to the printed precision."""
    cases = [((0.0, 0.0), (0.0, 1.0), 111.318077888),
             ((0.0, 0.0), (1.0, 0.0), 111.318077888),
             ((-17.6635, 9.166), (-17.5, 9.5), 41.291074146),
             ((100.0, -45.0), (-100.0, 45.0), 12659.746604378)]
    for p, q, expected in cases:
        assert bp.earth_chord_distances(p, q)[0, 0] == pytest.approx(expected, abs=5e-9)


def test_trace_D_analog_selection(wilsonlake):
    """Section D: Wilson Lake selects 24 cells at std(dats)*2."""
    tol = wilsonlake.std(ddof=1) * 2
    assert wilsonlake.size == 52
    assert wilsonlake.mean() == pytest.approx(0.7742307692, abs=1e-9)
    assert tol == pytest.approx(0.1754030050, abs=1e-9)
    means = stores.get_coretops("SST").cell_means("tex86")
    sel = bp.find_analogs(means, float(wilsonlake.mean()), tol)
    assert sel.size == 24
    assert means.min() == pytest.approx(0.312907, abs=1e-6)
    assert means.max() == pytest.approx(0.839395, abs=1e-6)


def test_trace_E_posterior(lopes):
    """Section E: post_mean and post_sig for dats(1), first five thinned draws."""
    tex, lon, lat = lopes
    draws = stores.get_draws("subT")
    ind = utils.matlab_thinning(draws.n_draws, 1000)
    cell = draws.cell_index(lon, lat)
    post_mean, post_sd = bp.solve_posterior(
        tex, draws.alpha[cell, ind], draws.beta[cell, ind], draws.tau2[ind],
        19.1673884681, 6.0)
    expected_mean = [18.9957842624, 18.9528845464, 19.1093308376, 18.5569763711,
                     18.1852227316]
    expected_sd = [2.5632716642, 2.5861909363, 2.5574532133, 2.6952836733, 2.7369667635]
    np.testing.assert_allclose(post_mean[0, :5], expected_mean, atol=1e-8)
    np.testing.assert_allclose(post_sd[:5], expected_sd, atol=1e-8)


def test_trace_F_percentiles():
    """Section F: prctile(1:10, [5 50 95]) is 1.0, 5.5, 10.0."""
    np.testing.assert_allclose(utils.prctile(np.arange(1.0, 11.0), [5, 50, 95]),
                               [1.0, 5.5, 10.0])
