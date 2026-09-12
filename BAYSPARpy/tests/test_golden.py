"""Compare the port against values generated in MATLAB.

These are the only tests that can catch a *shared* misreading of the reference:
everything else checks the port against itself or against intermediates this
port produced. They need `docs/BAYSPARpy/audit/golden_matlab.mat`, which is made
by running `docs/BAYSPARpy/audit/generate_golden.m` in MATLAB from the
repository root. Until that file exists the whole module skips, and says so.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import baysparpy as bp
from baysparpy import stores, utils

scipy_io = pytest.importorskip("scipy.io")


def _golden_path() -> Path | None:
    for parent in Path(__file__).resolve().parents:
        candidate = parent / "docs" / "BAYSPARpy" / "audit" / "golden_matlab.mat"
        if candidate.is_file():
            return candidate
    return None


@pytest.fixture(scope="module")
def golden() -> dict:
    """The MATLAB reference values, or a skip explaining how to make them."""
    path = _golden_path()
    if path is None:
        pytest.skip(
            "golden file not yet generated from MATLAB — run "
            "docs/BAYSPARpy/audit/generate_golden.m from the repository root"
        )
    return scipy_io.loadmat(str(path), squeeze_me=True, struct_as_record=False)


def _field(struct, name):
    """Read a field from a squeezed MATLAB struct."""
    return getattr(struct, name)


def test_golden_thinning(golden):
    """MATLAB's round(linspace(...)) against matlab_thinning()."""
    thinning = golden["thinning"]
    for n_draws in (10, 999, 1000, 1001, 20000):
        expected = np.atleast_1d(_field(thinning, f"n{n_draws}")).astype(int)
        got = utils.matlab_thinning(20000, n_draws) + 1        # compare 1-based
        np.testing.assert_array_equal(got, expected)


def test_golden_distances(golden):
    """EarthChordDistances_2 against the Python port, to 1e-9 km."""
    pairs = np.atleast_2d(golden["distance_pairs"])
    expected = np.atleast_1d(golden["distances"])
    for (lon1, lat1, lon2, lat2), want in zip(pairs, expected):
        got = bp.earth_chord_distances((lon1, lat1), (lon2, lat2))[0, 0]
        assert got == pytest.approx(want, abs=1e-9, rel=1e-12)


@pytest.mark.parametrize("runname", ["SST", "subT"])
@pytest.mark.parametrize("series", ["castaneda2010", "lopes_santos2010", "shevenell2011"])
def test_golden_prior_mean_and_cell(golden, runname, series):
    """Prior mean, its observation count, and the grid cell — all exact."""
    ref = _field(golden["prior"], f"{runname}_{series}")
    lon, lat = float(_field(ref, "lon")), float(_field(ref, "lat"))
    mean, n = stores.get_seatemp(runname).prior_mean(lon, lat, max_dist=500.0, min_num=1)
    assert n == int(_field(ref, "n_below"))
    assert mean == pytest.approx(float(_field(ref, "prior_mean")), abs=1e-10, rel=1e-13)
    draws = stores.get_draws(runname)
    cell = draws.cell_index(lon, lat, strict=True)
    np.testing.assert_array_equal(draws.locs[cell], np.atleast_1d(_field(ref, "grid_loc")))


def test_golden_analog_selection(golden):
    """Wilson Lake: the cell means, the tolerance, and the selected set."""
    ref = golden["analog"]
    means = stores.get_coretops("SST").cell_means("tex86")
    np.testing.assert_allclose(means, np.atleast_1d(_field(ref, "cell_means")), rtol=1e-13)
    tol = float(_field(ref, "search_tol"))
    sel = bp.find_analogs(means, float(_field(ref, "mean_tex")), tol)
    expected = np.atleast_1d(_field(ref, "selected")).astype(int) - 1   # 1-based in MATLAB
    np.testing.assert_array_equal(sel, expected)


def test_golden_posterior(golden):
    """post_mean and post_sig for dats(1), the first five thinned draws."""
    ref = golden["posterior"]
    lopes = scipy_io.loadmat(stores.find_modeloutput() / "tex_testdata.mat")
    rec = lopes["lopes_santos2010"][0, 0]
    tex = np.ravel(rec["tex86"]).astype(float)
    lon, lat = float(np.ravel(rec["lon"])[0]), float(np.ravel(rec["lat"])[0])
    draws = stores.get_draws("subT")
    ind = utils.matlab_thinning(draws.n_draws, 1000)
    cell = draws.cell_index(lon, lat)
    post_mean, post_sd = bp.solve_posterior(
        tex, draws.alpha[cell, ind], draws.beta[cell, ind], draws.tau2[ind],
        float(_field(ref, "prior_mean")), float(_field(ref, "prior_std")))
    np.testing.assert_allclose(post_mean[0, :5],
                               np.atleast_1d(_field(ref, "post_mean")), rtol=1e-12)
    np.testing.assert_allclose(post_sd[:5],
                               np.atleast_1d(_field(ref, "post_sig")), rtol=1e-12)


def test_golden_prctile(golden):
    """prctile parity: Hazen plotting positions, not NumPy's default."""
    expected = np.atleast_1d(golden["prctile_1to10"])
    np.testing.assert_allclose(utils.prctile(np.arange(1.0, 11.0), [5, 50, 95]),
                               expected, rtol=1e-12)
    assert not np.allclose(np.percentile(np.arange(1.0, 11.0), [5, 50, 95]), expected)
