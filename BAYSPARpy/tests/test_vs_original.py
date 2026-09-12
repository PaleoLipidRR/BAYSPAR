"""Compare the port against the ORIGINAL MATLAB functions actually executing.

`test_golden.py` checks the arithmetic against intermediates recomputed inline;
this checks it against `bayspar_tex.m`, `bayspar_tex_analog.m` and
`TEX_forward.m` themselves, run unmodified. Deterministic fields must match
exactly; the draws are compared distributionally, since the two RNG streams
cannot be made to agree (PORTING.md GEN-03).

Generate the reference with `docs/BAYSPARpy/audit/verify_original.m` (MATLAB, or
Octave 8+ with the statistics package). Without it this module skips.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import baysparpy as bp

scipy_io = pytest.importorskip("scipy.io")

# The per-cell standard deviation is estimated from these same runs, so the test
# statistic is Student-t with len(SEEDS)-1 degrees of freedom, not normal. With
# nine seeds that mattered: t_8 gives mean |t| = 0.86 and P(|t| > 2) = 0.081,
# against 0.80 and 0.046 for a normal, which reads as over-dispersion when it is
# only the estimated denominator. Twenty-five brings t_24 close enough to normal
# (0.81 and 0.057) that the thresholds below mean what they say.
SEEDS = tuple(range(101, 126))

# One MATLAB run against the mean of len(SEEDS) Python runs, so the difference has
# sd = mcse * sqrt(1 + 1/n). Judged per cell at 3 sigma a 579-cell comparison
# would fail by chance more often than not; 5 sigma puts the family-wise false
# alarm near 3e-4 while still catching any real offset, which would show up as a
# shift across many cells rather than one.
Z_MAX = 5.0


def _reference_path() -> Path | None:
    for parent in Path(__file__).resolve().parents:
        candidate = parent / "docs" / "BAYSPARpy" / "audit" / "original_matlab.mat"
        if candidate.is_file():
            return candidate
    return None


@pytest.fixture(scope="module")
def original() -> dict:
    """Output of the original MATLAB functions, or a skip explaining how to make it."""
    path = _reference_path()
    if path is None:
        pytest.skip(
            "original_matlab.mat not found — run docs/BAYSPARpy/audit/verify_original.m "
            "from the repository root"
        )
    return scipy_io.loadmat(str(path), squeeze_me=True, struct_as_record=False)


def _mc_envelope(runs: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """Mean across seeds, and the seed-to-seed standard deviation."""
    stack = np.stack(runs)
    return stack.mean(axis=0), stack.std(axis=0, ddof=1)


def test_standard_deterministic_fields(original):
    """Prior mean and grid cell must match bayspar_tex.m exactly."""
    ref = original["standard"]
    inp = original["standard_input"]
    out = bp.bayspar_tex(np.atleast_1d(inp.tex86), float(inp.lon), float(inp.lat),
                         float(inp.prior_std), str(inp.runname),
                         n_draws=int(original["n_draws"]), seed=0)
    assert out.prior_mean == pytest.approx(float(ref.PriorMean), abs=1e-10, rel=1e-13)
    np.testing.assert_array_equal(np.asarray(out.grid_loc), np.atleast_1d(ref.GridLoc))
    np.testing.assert_allclose(np.asarray(out.site_loc), np.atleast_1d(ref.SiteLoc))


def _z_scores(runs: list[np.ndarray], matlab: np.ndarray, floor: float) -> np.ndarray:
    """Standardised difference between the Python mean and the single MATLAB run."""
    mean, mcse = _mc_envelope(runs)
    assert mean.shape == matlab.shape, f"{mean.shape} vs {matlab.shape}"
    sd = np.maximum(mcse, floor) * np.sqrt(1.0 + 1.0 / len(runs))
    return (mean - matlab) / sd


def _assert_agrees(z: np.ndarray, what: str) -> None:
    """No cell beyond Z_MAX, and the bulk distributed like standard normal.

    The second check is what would catch a small systematic offset -- a bias too
    small to push any single cell past 5 sigma still lifts the mean |z| well above
    the 0.8 a standard normal gives.
    """
    worst = np.unravel_index(np.abs(z).argmax(), z.shape)
    assert np.abs(z).max() < Z_MAX, (
        f"{what}: |z| = {np.abs(z).max():.2f} at {worst}, beyond {Z_MAX}")
    assert np.abs(z).mean() < 1.2, (
        f"{what}: mean |z| = {np.abs(z).mean():.2f}; sampling noise alone gives about "
        "0.81, so this looks like a systematic offset")


def test_standard_percentiles_within_monte_carlo_error(original):
    """The 5/50/95 must agree with bayspar_tex.m to within the sampling noise."""
    ref = original["standard"]
    inp = original["standard_input"]
    tex = np.atleast_1d(inp.tex86)
    runs = [bp.bayspar_tex(tex, float(inp.lon), float(inp.lat), float(inp.prior_std),
                           str(inp.runname), n_draws=int(original["n_draws"]),
                           seed=s).preds for s in SEEDS]
    z = _z_scores(runs, np.atleast_2d(ref.Preds), floor=0.01)
    _assert_agrees(z, "bayspar_tex")


def test_prctile_convention_on_the_reference_ensemble(original):
    """numpy 'hazen' over MATLAB's own ensemble must return MATLAB's own Preds.

    This is the sharpest test of PORTING.md BT-11 available without MATLAB in the
    loop: same numbers in, same percentiles out, no sampling noise anywhere. It
    reproduces all 579 values exactly, and NumPy's default method does not.
    """
    ref = original["standard"]
    ens = np.asarray(ref.PredsEns)
    saved = np.atleast_2d(ref.Preds)
    from baysparpy.utils import prctile
    np.testing.assert_array_equal(prctile(ens, [5, 50, 95], axis=1), saved)
    assert not np.allclose(np.percentile(ens, [5, 50, 95], axis=1).T, saved)


def test_analog_locations_match_exactly(original):
    """The analogue set must be identical -- it is a deterministic selection."""
    ref = original["analog"]
    inp = original["analog_input"]
    out = bp.bayspar_tex_analog(np.atleast_1d(inp.tex86), float(inp.prior_mean),
                                float(inp.prior_std), float(inp.search_tol),
                                str(inp.runname), n_draws=int(original["n_draws"]),
                                mode="matlab", seed=0)
    np.testing.assert_array_equal(out.analog_locs, np.atleast_2d(ref.AnLocs))


def test_analog_percentiles_within_monte_carlo_error(original):
    """Analogue-mode percentiles, in matlab pairing mode for like-for-like."""
    ref = original["analog"]
    inp = original["analog_input"]
    tex = np.atleast_1d(inp.tex86)
    runs = [bp.bayspar_tex_analog(tex, float(inp.prior_mean), float(inp.prior_std),
                                  float(inp.search_tol), str(inp.runname),
                                  n_draws=int(original["n_draws"]), mode="matlab",
                                  seed=s).preds for s in SEEDS]
    z = _z_scores(runs, np.atleast_2d(ref.Preds), floor=0.02)
    _assert_agrees(z, "bayspar_tex_analog")


def test_analog_ensemble_shape_confirms_BTA_08(original):
    """The reference really does return a flattened 2-D ensemble.

    PORTING.md BTA-08 says the MATLAB docstring promises
    (Nd, n_analogues, Nsamps) and the code returns (Nd, n_analogues * Nsamps).
    This asserts it against the file the reference implementation wrote, so the
    entry rests on execution rather than on reading.
    """
    ref = original["analog"]
    n_an = np.atleast_2d(ref.AnLocs).shape[0]
    n_draws = int(original["n_draws"])
    n_obs = np.atleast_2d(ref.Preds).shape[0]
    # the reference's own ensemble shape, recorded rather than stored: the array
    # itself is 10 MB and only its shape is under test
    assert tuple(np.atleast_1d(ref.PredsEns_size).astype(int)) == (n_obs, n_an * n_draws)
    out = bp.bayspar_tex_analog(np.atleast_1d(original["analog_input"].tex86), 30.0, 20.0,
                                float(original["analog_input"].search_tol), "SST",
                                n_draws=n_draws, save_ensemble=True, seed=0)
    assert out.ensemble.shape == (n_obs, n_an, n_draws)      # the documented shape


def test_forward_within_monte_carlo_error(original):
    """TEX_forward.m: per-temperature ensemble means."""
    ref = np.atleast_2d(original["forward"])
    inp = original["forward_input"]
    t = np.atleast_1d(inp.t)
    lat, lon = float(inp.lat), float(inp.lon)
    runs = [bp.tex_forward(np.full(t.size, lat), np.full(t.size, lon), t,
                           str(inp.runname), seed=s).tex86.mean(axis=1) for s in SEEDS]
    z = _z_scores(runs, ref.mean(axis=1), floor=0.0005)
    _assert_agrees(z, "TEX_forward")
