"""One test per entry in docs/BAYSPARpy/PORTING.md, named for its ID.

tests/test_porting_doc.py asserts the correspondence in both directions, so an
entry cannot be added to the record without a test, or a test renamed without
the record noticing.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest
import scipy.io as sio

import baysparpy as bp
from baysparpy import constants, stores, utils
from baysparpy.errors import (
    AmbiguousGridCellError, EnsembleSizeError, SearchToleranceError,
)

# ---------------------------------------------------------------- conventions


def test_GEN_01(modeloutput):
    """Index base: Inds_Stack is 1-based on disk, 0-based in memory."""
    raw = sio.loadmat(modeloutput / "Data_Input_SpatAg_SST.mat")["Data_Input"][0, 0]
    on_disk = np.ravel(raw["Inds_Stack"]).astype(int)
    ct = stores.get_coretops("SST")
    assert on_disk.min() == 1 and on_disk.max() == len(ct.locs)
    assert ct.cell_index.min() == 0 and ct.cell_index.max() == len(ct.locs) - 1
    np.testing.assert_array_equal(ct.cell_index, on_disk - 1)


def test_GEN_02():
    """Column vectors: (:) for 1-D, a refusal for genuinely 2-D input."""
    np.testing.assert_array_equal(utils.as_column([[1.0], [2.0], [3.0]]), [1, 2, 3])
    np.testing.assert_array_equal(utils.as_column(np.array([[1.0, 2.0, 3.0]])), [1, 2, 3])
    with pytest.raises(ValueError, match="column-major"):
        utils.as_column(np.ones((2, 3)))


def test_GEN_03(lopes):
    """RNG: a seed reproduces exactly; no seed does not."""
    tex, lon, lat = lopes
    a = bp.bayspar_tex(tex, lon, lat, 6.0, "subT", n_draws=200, seed=7)
    b = bp.bayspar_tex(tex, lon, lat, 6.0, "subT", n_draws=200, seed=7)
    c = bp.bayspar_tex(tex, lon, lat, 6.0, "subT", n_draws=200, seed=8)
    np.testing.assert_array_equal(a.preds, b.preds)
    assert not np.array_equal(a.preds, c.preds)


def test_GEN_04():
    """Types: float64 throughout, and the int16/uint8 fields cast before use."""
    d = stores.get_draws("SST")
    assert d.alpha.dtype == np.float64 and d.locs.dtype == np.float64
    ct = stores.get_coretops("SST")
    assert ct.cell_index.dtype.kind == "i"           # not uint8: it must go negative-safe
    assert np.issubdtype(ct.tex86.dtype, np.floating)


def test_GEN_05(lopes):
    """Broadcasting equals MATLAB's repmat, to float tolerance."""
    tex, _, _ = lopes
    rng = np.random.default_rng(0)
    alpha, beta = rng.normal(0.3, 0.01, 50), rng.normal(0.015, 0.001, 50)
    tau2 = rng.normal(0.0018, 1e-4, 50) ** 2 + 1e-6
    got_mean, got_sd = bp.solve_posterior(tex, alpha, beta, tau2, 20.0, 6.0)
    n, m = len(tex), len(alpha)
    pmu = np.tile(np.ones((n, 1)) * 20.0, (1, m))            # repmat, literally
    pinv = np.tile(6.0, (n, m)) ** -2
    sig = np.sqrt(tau2)
    num = pinv * pmu + np.tile(sig, (n, 1)) ** -2 * np.tile(beta, (n, 1)) * (
        tex[:, None] - np.tile(alpha, (n, 1)))
    den = pinv + np.tile(beta, (n, 1)) ** 2 * np.tile(sig, (n, 1)) ** -2
    np.testing.assert_allclose(got_mean, num / den, rtol=1e-12)
    np.testing.assert_allclose(np.broadcast_to(got_sd, (n, m)), np.sqrt(den ** -1), rtol=1e-12)


def test_GEN_06():
    """length(): cells and draws are named separately, never max(size(...))."""
    d = stores.get_draws("SST")
    assert d.alpha.shape == (162, 20000)
    assert d.n_draws == 20000                     # not 162, and not the documented 15000
    assert len(d.locs) == 162


def test_GEN_07(wilsonlake):
    """Error messages are MATLAB's, verbatim."""
    with pytest.raises(SearchToleranceError, match="^Your search tolerance is too narrow$"):
        bp.bayspar_tex_analog(wilsonlake, 30.0, 20.0, 1e-6, "SST", n_draws=10)
    with pytest.raises(ValueError, match="^To use analog mode, enter a search tolerance"):
        bp.tex_forward(0.0, 0.0, [20.0], "SST", type="analog")
    with pytest.raises(ValueError, match='^please enter "analog" to specify analog mode$'):
        bp.tex_forward(0.0, 0.0, [20.0], "SST", type="nonsense")


# ------------------------------------------------- EarthChordDistances_2.m


def test_ECD_01():
    """Earth radius is the MATLAB constant."""
    assert constants.EARTH_RADIUS_KM == 6378.137


def test_ECD_02():
    """Broadcasting reproduces the kron pairing for every (i, j)."""
    rng = np.random.default_rng(0)
    p1 = np.column_stack([rng.uniform(-180, 180, 4), rng.uniform(-90, 90, 4)])
    p2 = np.column_stack([rng.uniform(-180, 180, 6), rng.uniform(-90, 90, 6)])
    full = bp.earth_chord_distances(p1, p2)
    for i in range(4):
        for j in range(6):
            assert full[i, j] == pytest.approx(
                bp.earth_chord_distances(p1[i], p2[j])[0, 0], rel=1e-15)


def test_ECD_03():
    """The formula, against values computed from the MATLAB source."""
    assert bp.earth_chord_distances((0, 0), (0, 1))[0, 0] == pytest.approx(
        111.318077888, abs=1e-8)
    assert bp.earth_chord_distances((0, 0), (1, 0))[0, 0] == pytest.approx(
        111.318077888, abs=1e-8)
    assert bp.earth_chord_distances((100, -45), (-100, 45))[0, 0] == pytest.approx(
        12659.746604378, abs=1e-6)
    # the chord is shorter than the great-circle arc, by 0.026% at 500 km
    d = bp.earth_chord_distances((0, 0), (0, 4.4931))[0, 0]
    arc = constants.EARTH_RADIUS_KM * np.deg2rad(4.4931)
    assert d < arc and (arc - d) / arc == pytest.approx(2.6e-4, rel=0.2)


def test_ECD_04():
    """Output is (N, M): rows index points1, columns points2."""
    p1 = np.zeros((3, 2))
    p2 = np.zeros((7, 2))
    assert bp.earth_chord_distances(p1, p2).shape == (3, 7)


def test_ECD_05():
    """Identical points are zero; the dateline pair is float noise, as in MATLAB."""
    assert bp.earth_chord_distances((12.5, -33.25), (12.5, -33.25))[0, 0] == 0.0
    # 180E to 180W puts sin(pi) in the expression, which is 1.2e-16 rather than 0,
    # so the distance is ~1.6e-12 km. The MATLAB evaluates the same expression.
    assert bp.earth_chord_distances((180, 0), (-180, 0))[0, 0] == pytest.approx(0, abs=1e-9)
    assert bp.earth_chord_distances((1, 2), (3, 4)).shape == (1, 1)


# ------------------------------------------------------------ bayspar_tex.m


def test_BT_01(lopes):
    """Arguments: MATLAB's defaults, plus a real check on n_draws."""
    tex, lon, lat = lopes
    out = bp.bayspar_tex(tex, lon, lat, 6.0, "subT", seed=0)
    assert out.metadata["n_draws"] == constants.DEFAULT_N_DRAWS == 1000
    assert out.ensemble is None                         # ens_sel defaults to 0
    with pytest.raises(EnsembleSizeError):
        bp.bayspar_tex(tex, lon, lat, 6.0, "subT", n_draws=20001, seed=0)


def test_BT_02():
    """The store is cached and location-independent, unlike MATLAB's load()."""
    stores.get_draws("SST")
    before = stores._load_mat.cache_info()
    second = stores.get_draws("SST")
    after = stores._load_mat.cache_info()
    assert after.hits == before.hits + 1 and after.misses == before.misses
    assert np.shares_memory(second.alpha, stores.get_draws("SST").alpha)
    assert (stores.find_modeloutput() / "obsSST.mat").is_file()


def test_BT_03(lopes):
    """runname is validated up front; lower-case aliases are accepted."""
    tex, lon, lat = lopes
    a = bp.bayspar_tex(tex, lon, lat, 6.0, "subT", n_draws=50, seed=3)
    b = bp.bayspar_tex(tex, lon, lat, 6.0, "subt", n_draws=50, seed=3)
    np.testing.assert_array_equal(a.preds, b.preds)
    with pytest.raises(ValueError, match="runname must be one of"):
        bp.bayspar_tex(tex, lon, lat, 6.0, "SSTT")


def test_BT_04():
    """Grid half-space is 10 degrees, and the store tiles the globe with it."""
    assert constants.GRID_HALF_SPACE == 10.0
    locs = stores.get_draws("SST").locs
    assert sorted(np.unique(locs[:, 0])) == list(range(-170, 180, 20))
    assert sorted(np.unique(locs[:, 1])) == list(range(-80, 90, 20))


def test_BT_05(lopes):
    """The prior-mean search radius is MATLAB's by default and adjustable."""
    tex, lon, lat = lopes
    assert (constants.MAX_DIST_KM, constants.MIN_NUM) == (500.0, 1)
    wide = bp.bayspar_tex(tex, lon, lat, 6.0, "subT", n_draws=50, seed=0,
                          max_dist=1500.0)
    default = bp.bayspar_tex(tex, lon, lat, 6.0, "subT", n_draws=50, seed=0)
    assert wide.metadata["n_obs_within_max_dist"] > default.metadata["n_obs_within_max_dist"]
    assert wide.prior_mean != default.prior_mean


def test_BT_06():
    """Thinning spans the chain, with MATLAB's half-away-from-zero rounding."""
    ind = utils.matlab_thinning(20000, 1000) + 1          # report 1-based
    np.testing.assert_array_equal(ind[:5], [1, 21, 41, 61, 81])
    assert ind[-1] == 20000
    np.testing.assert_array_equal(utils.matlab_thinning(20000, 10) + 1,
                                  [1, 2223, 4445, 6667, 8889, 11112, 13334, 15556,
                                   17778, 20000])
    # linspace(1, 4, 3) lands on 2.5: MATLAB rounds away from zero (3, index 2),
    # np.round rounds to even (2, index 1). This is the whole point of the helper.
    assert utils.matlab_thinning(4, 3)[1] == 2
    assert int(np.round(np.linspace(1, 4, 3)[1])) - 1 == 1
    np.testing.assert_array_equal(utils.matlab_thinning(20000, 20000),
                                  np.arange(20000))


def test_BT_07(lopes):
    """Prior mean: value, count, and the fallback when nothing is in range."""
    _, lon, lat = lopes
    obs = stores.get_seatemp("subT")
    mean, n = obs.prior_mean(lon, lat, max_dist=500.0, min_num=1)
    assert mean == pytest.approx(19.1673884681, abs=1e-9)
    assert n == 56
    # nothing within 1 m of a mid-ocean point -> the single closest observation
    tiny, n0 = obs.prior_mean(lon, lat, max_dist=1e-3, min_num=1)
    assert n0 == 0
    d = bp.earth_chord_distances((lon, lat), obs.locs)[0]
    assert tiny == pytest.approx(obs.values[int(np.argmin(d))])


def test_BT_08():
    """Boundary sites: nearest centroid with a warning, or a refusal."""
    d = stores.get_draws("SST")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert d.locs[d.cell_index(-17.6635, 9.166)].tolist() == [-10.0, 0.0]
    with pytest.warns(UserWarning, match="cell boundary"):
        cell = d.cell_index(0.0, 0.0)
    assert abs(d.locs[cell][0]) == 10.0                    # one of the two neighbours
    with pytest.raises(AmbiguousGridCellError, match="cell boundary"):
        d.cell_index(0.0, 0.0, strict=True)
    with pytest.raises(ValueError, match="outside the calibration grid"):
        d.cell_index(200.0, 0.0)


def test_BT_09(lopes):
    """The conjugate update, against values traced from the MATLAB expression."""
    tex, lon, lat = lopes
    d = stores.get_draws("subT")
    ind = utils.matlab_thinning(d.n_draws, 1000)
    cell = d.cell_index(lon, lat)
    post_mean, post_sd = bp.solve_posterior(
        tex, d.alpha[cell, ind], d.beta[cell, ind], d.tau2[ind], 19.1673884681, 6.0)
    assert post_mean[0, 0] == pytest.approx(18.9957842624, abs=1e-8)
    assert post_sd[0] == pytest.approx(2.5632716642, abs=1e-8)
    assert post_mean[0, 4] == pytest.approx(18.1852227316, abs=1e-8)
    assert post_sd.shape == (1000,)          # denominator does not depend on the data


def test_BT_10(lopes):
    """The draw is post_mean + N(0,1)*post_sd, independent per point and draw."""
    tex, lon, lat = lopes
    out = bp.bayspar_tex(tex, lon, lat, 6.0, "subT", n_draws=4000,
                         save_ensemble=True, seed=11)
    d = stores.get_draws("subT")
    ind = utils.matlab_thinning(d.n_draws, 4000)
    cell = d.cell_index(lon, lat)
    post_mean, post_sd = bp.solve_posterior(
        tex, d.alpha[cell, ind], d.beta[cell, ind], d.tau2[ind], out.prior_mean, 6.0)
    z = (out.ensemble - post_mean) / post_sd
    assert abs(z.mean()) < 0.02 and z.std() == pytest.approx(1.0, abs=0.02)


def test_BT_11():
    """prctile uses Hazen plotting positions, not NumPy's default."""
    x = np.arange(1.0, 11.0)
    np.testing.assert_allclose(utils.prctile(x, [5, 50, 95]), [1.0, 5.5, 10.0])
    assert not np.allclose(np.percentile(x, [5, 50, 95]), [1.0, 5.5, 10.0])
    got = utils.prctile(np.arange(1.0, 5.0)[None, :], [25, 50, 75], axis=1)
    np.testing.assert_allclose(got[0], [1.5, 2.5, 3.5])


def test_BT_12(lopes):
    """The result carries MATLAB's fields, plus provenance."""
    tex, lon, lat = lopes
    out = bp.bayspar_tex(tex, lon, lat, 6.0, "subT", n_draws=100,
                         save_ensemble=True, seed=0)
    assert out.preds.shape == (len(tex), 3)
    np.testing.assert_array_equal(out.p50, out.preds[:, 1])
    assert (out.p5 <= out.p50).all() and (out.p50 <= out.p95).all()
    assert out.ensemble.shape == (len(tex), 100)
    d = out.to_dict()
    assert set(d) >= {"Preds", "SiteLoc", "GridLoc", "PriorMean", "PriorStd", "PredsEns"}
    assert out.metadata["thinning_index"].size == 100


# ----------------------------------------------------- bayspar_tex_analog.m


def test_BTA_01():
    """Analogue parameters come from the standard store, by cell index."""
    d = stores.get_draws("SST")
    ct = stores.get_coretops("SST")
    rows = stores.analog_cell_rows(d, ct)
    pa = sio.loadmat(stores.find_modeloutput() / "Output_SpatAg_SST" / "params_analog.mat")
    np.testing.assert_array_equal(d.alpha[rows], pa["alpha_samples"])
    np.testing.assert_array_equal(d.beta[rows], pa["beta_samples"])


def test_BTA_02(modeloutput):
    """Cell means match the MATLAB loop, and average observations not sites."""
    ct = stores.get_coretops("SST")
    di = sio.loadmat(modeloutput / "Data_Input_SpatAg_SST.mat")["Data_Input"][0, 0]
    inds = np.ravel(di["Inds_Stack"]).astype(int)
    obs = np.ravel(di["Obs_Stack"]).astype(float)
    expected = np.array([obs[inds == i + 1].mean() for i in range(80)])
    np.testing.assert_allclose(ct.cell_means("tex86"), expected, rtol=1e-15)


def test_BTA_03(wilsonlake):
    """Analogue selection is inclusive, and matches the demo's 24 cells."""
    ct = stores.get_coretops("SST")
    means = ct.cell_means("tex86")
    tol = wilsonlake.std(ddof=1) * 2
    sel = bp.find_analogs(means, wilsonlake.mean(), tol)
    assert sel.size == 24
    assert (means[sel] >= wilsonlake.mean() - tol).all()
    assert (means[sel] <= wilsonlake.mean() + tol).all()
    # inclusive at both ends: a tolerance landing exactly on a cell mean keeps it
    exact = bp.find_analogs(means, float(means[0]), 0.0)
    assert 0 in exact


def test_BTA_04():
    """An empty analogue set raises with MATLAB's message."""
    with pytest.raises(SearchToleranceError):
        bp.find_analogs(np.array([0.5, 0.6]), 0.9, 0.01)


def test_BTA_05(wilsonlake):
    """tau2 pairing: 'matlab' reproduces the mis-pairing, 'modern' does not.

    The difference on this record is within Monte-Carlo noise; the test records
    that rather than asserting the two are identical.
    """
    kw = dict(n_draws=800, seed=5, save_ensemble=True)
    modern = bp.bayspar_tex_analog(wilsonlake, 30.0, 20.0, 0.18, "SST", **kw)
    matlab = bp.bayspar_tex_analog(wilsonlake, 30.0, 20.0, 0.18, "SST", mode="matlab", **kw)
    n_an = len(modern.analog_locs)
    assert modern.metadata["mode"] == "modern" and matlab.metadata["mode"] == "matlab"
    # the pooled parameter sets are the same size; only the pairing differs
    assert modern.ensemble.reshape(len(wilsonlake), -1).shape == matlab.ensemble.shape
    assert np.abs(modern.preds - matlab.preds).max() < 0.5      # noise-scale, not a shift
    with pytest.raises(ValueError, match="mode must be"):
        bp.bayspar_tex_analog(wilsonlake, 30.0, 20.0, 0.18, "SST", mode="fast")
    assert n_an > 1                                   # the mis-pairing only bites here


def test_BTA_06(wilsonlake):
    """Pooling order: the location axis is explicit, not folded into the draws."""
    out = bp.bayspar_tex_analog(wilsonlake, 30.0, 20.0, 0.18, "SST", n_draws=64,
                                save_ensemble=True, seed=2)
    n_an = len(out.analog_locs)
    assert out.ensemble.shape == (len(wilsonlake), n_an, 64)
    flat = bp.bayspar_tex_analog(wilsonlake, 30.0, 20.0, 0.18, "SST", n_draws=64,
                                 save_ensemble=True, seed=2, mode="matlab")
    assert flat.ensemble.shape == (len(wilsonlake), n_an * 64)


def test_BTA_07(wilsonlake):
    """Analogue locations are returned in Data_Input.Locs order."""
    out = bp.bayspar_tex_analog(wilsonlake, 30.0, 20.0, 0.18, "SST", n_draws=32, seed=0)
    ct = stores.get_coretops("SST")
    np.testing.assert_array_equal(out.analog_locs, ct.locs[out.metadata["analog_cells"]])
    assert out.analog_locs.shape[1] == 2


def test_BTA_08(wilsonlake):
    """The ensemble is (N, analogues, draws) -- the shape MATLAB documents."""
    out = bp.bayspar_tex_analog(wilsonlake, 30.0, 20.0, 0.18, "SST", n_draws=40,
                                save_ensemble=True, seed=0)
    assert out.ensemble.shape == (len(wilsonlake), len(out.analog_locs), 40)
    # axis 1 lines up with analog_locs: a colder analogue gives colder predictions
    per_loc = np.median(out.ensemble, axis=(0, 2))
    assert per_loc.shape == (len(out.analog_locs),)


# ---------------------------------------------------------- TEX_forward.m


def test_TXF_01():
    """Arguments: (lat, lon, t) order, MATLAB's, and the analogue requirement."""
    out = bp.tex_forward(9.166, -17.6635, [25.0], "SST", n_draws=2000, seed=0)
    assert out.tex86.shape == (1, 1000)
    with pytest.raises(ValueError, match="search tolerance"):
        bp.tex_forward(9.166, -17.6635, [25.0], "SST", type="analog")


def test_TXF_02():
    """Mode is validated before any file is read."""
    with pytest.raises(ValueError, match='"analog"'):
        bp.tex_forward(0.0, 0.0, [20.0], "SST", type="Standard")


def test_TXF_03():
    """Several locations pair elementwise with temperatures, not as a grid."""
    ok = bp.tex_forward([9.166, 31.65], [-17.66, 34.07], [25.0, 22.0], "SST",
                        n_draws=2000, seed=0)
    assert ok.tex86.shape == (2, 1000)
    with pytest.raises(ValueError, match="one temperature per location"):
        bp.tex_forward([9.166, 31.65], [-17.66, 34.07], [25.0, 22.0, 20.0], "SST",
                       n_draws=2000, seed=0)
    single = bp.tex_forward(9.166, -17.6635, [25.0, 22.0, 20.0], "SST",
                            n_draws=2000, seed=0)
    assert single.tex86.shape == (3, 1000)


def test_TXF_04():
    """The forward analogue search is in degrees C, against Target_Stack."""
    ct = stores.get_coretops("SST")
    t_means = ct.cell_means("target_t")
    tex_means = ct.cell_means("tex86")
    assert t_means.max() > 20 and tex_means.max() < 1.0       # different fields entirely
    out = bp.tex_forward(0.0, 0.0, [25.0], "SST", type="analog", search_tol=2.0,
                         n_draws=2000, seed=0)
    sel = bp.find_analogs(t_means, 25.0, 2.0)
    assert out.tex86.shape == (1, 1000) and sel.size > 1


def test_TXF_05():
    """No thinning by default; thinning available."""
    full = bp.tex_forward(9.166, -17.6635, [25.0], "SST", seed=0)
    assert full.metadata["n_draws"] == 20000
    thin = bp.tex_forward(9.166, -17.6635, [25.0], "SST", n_draws=1500, seed=0)
    assert thin.metadata["n_draws"] == 1500


def test_TXF_06():
    """The draw is normal about alpha + beta*t with sd sqrt(tau2)."""
    d = stores.get_draws("SST")
    cell = d.cell_index(-17.6635, 9.166)
    expected = (d.alpha[cell] + d.beta[cell] * 25.0).mean()
    out = bp.tex_forward(9.166, -17.6635, [25.0], "SST", seed=4, n_out=20000)
    assert out.tex86.mean() == pytest.approx(expected, abs=0.01)


def test_TXF_07():
    """Output width follows n_out; MATLAB's 1000 is the default."""
    assert constants.FORWARD_N_OUT == 1000
    assert bp.tex_forward(9.166, -17.6635, [25.0], "SST", seed=0).tex86.shape[1] == 1000
    assert bp.tex_forward(9.166, -17.6635, [25.0], "SST", seed=0,
                          n_out=250).tex86.shape[1] == 250
    with pytest.raises(ValueError, match="exceeds"):
        bp.tex_forward(9.166, -17.6635, [25.0], "SST", n_draws=100, n_out=101, seed=0)


def test_TXF_08():
    """TEX86 is clipped to [0, 1], and the clipped share is reported."""
    cold = bp.tex_forward(-64.85, -64.21, [-30.0], "subT", seed=0, n_out=2000)
    assert cold.tex86.min() >= 0.0 and cold.tex86.max() <= 1.0
    assert cold.metadata["clipped_fraction"] > 0.0
    warm = bp.tex_forward(9.166, -17.6635, [25.0], "SST", seed=0, n_out=2000)
    assert warm.metadata["clipped_fraction"] == 0.0


# ------------------------------------------------------------------- stores


def test_STO_01():
    """params_analog is exactly the coretop rows of params_standard."""
    for runname in ("SST", "subT"):
        d = stores.get_draws(runname)
        ct = stores.get_coretops(runname)
        rows = stores.analog_cell_rows(d, ct)
        pa = sio.loadmat(stores.find_modeloutput() / f"Output_SpatAg_{runname}"
                         / "params_analog.mat")
        assert np.array_equal(d.alpha[rows], pa["alpha_samples"])
        assert np.array_equal(d.beta[rows], pa["beta_samples"])
        assert np.array_equal(d.tau2, np.ravel(pa["tau2_samples"]))


def test_STO_02():
    """NetCDF conversion with provenance attributes."""
    pytest.skip("not implemented: SPEC Phase 0, tools/convert_modeloutput.py")


def test_STO_03():
    """The two runs have different observation grids."""
    sst, subt = stores.get_seatemp("SST"), stores.get_seatemp("subT")
    assert sst.locs.shape == (37686, 2) and subt.locs.shape == (37105, 2)
    assert sst.values.size == 37686 and subt.values.size == 37105


def test_STO_04():
    """Coretops load as one row per observation, with the target and its error."""
    for runname, n_obs in (("SST", 903), ("subT", 906)):
        ct = stores.get_coretops(runname)
        assert ct.tex86.size == ct.target_t.size == ct.target_t_sd.size == n_obs
        assert ct.locs.shape == (80, 2)
        assert (ct.target_t_sd >= 0).all()
        assert 0.0 < ct.tex86.min() and ct.tex86.max() <= 1.0
    # 7 of the 903 SST targets carry an error SD of exactly zero, and none of the
    # subT ones do. An errors-in-variables likelihood divides by this, so the Stan
    # refit (SPEC 10.3) has to floor or special-case them. Pinned so the fix is
    # driven by data, not by a crash.
    assert int((stores.get_coretops("SST").target_t_sd == 0).sum()) == 7
    assert int((stores.get_coretops("subT").target_t_sd == 0).sum()) == 0


def test_STO_05():
    """Demo series exported as CSV, matching brews/baysparpy's example files."""
    pytest.skip("not implemented: SPEC Phase 0, demo data export")
