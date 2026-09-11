"""Prototype of the closed-form BAYSPAR standard prediction, benchmarked against
the dense per-draw solve used by brews/baysparpy. Backs ../SPEC.md section 5.

Run from the repository root:

    python docs/BAYSPARpy/audit/bench_closed_form.py

Requires numpy and scipy only. Timings are machine-dependent; the ratios are not.
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import scipy.io as sio

ROOT = Path(__file__).resolve().parents[3]
MO = ROOT / "ModelOutput"

RUNNAME = "subT"
SERIES = "lopes_santos2010"
PRIOR_STD = 6.0
N_DRAWS = 1000
EARTH_RADIUS_KM = 6378.137
GRID_HALF_SPACE = 10.0
MAX_DIST_KM = 500.0
MIN_NUM = 1


def earth_chord_distances(point, points):
    """Chordal distance in km, EarthChordDistances_2.m vectorised over `points`."""
    lon0, lat0 = point
    d2r = np.pi / 180.0
    half = np.arcsin(np.sqrt(
        np.sin((lat0 - points[:, 1]) * d2r / 2) ** 2
        + np.cos(lat0 * d2r) * np.cos(points[:, 1] * d2r)
        * np.sin(np.abs(lon0 - points[:, 0]) * d2r / 2) ** 2))
    return 2 * EARTH_RADIUS_KM * np.sin(half)


def matlab_thinning(n_total, n_draws):
    """round(linspace(1, n_total, n_draws)) with MATLAB's half-away-from-zero round."""
    return np.floor(np.linspace(1, n_total, n_draws) + 0.5).astype(int) - 1


def prior_mean_for(lon, lat, locs_obs, st_obs):
    """bayspar_tex.m prior mean: all obs within 500 km, else the single closest."""
    d = earth_chord_distances((lon, lat), locs_obs)
    order = np.argsort(d, kind="stable")
    n_below = int(np.searchsorted(d[order], MAX_DIST_KM, side="left"))
    take = n_below if n_below > MIN_NUM else MIN_NUM
    return st_obs[order[:take]].mean(), n_below


def closed_form(tex, alpha, beta, tau2, prior_mean, prior_std, rng):
    """Vectorised conjugate solve — no BLAS, one core. Returns (N, M)."""
    pinv_cov = prior_std ** -2
    den = pinv_cov + beta ** 2 / tau2                     # (M,)
    num = pinv_cov * prior_mean + (beta / tau2) * (tex[:, None] - alpha)
    post_mean = num / den
    post_sd = np.sqrt(1.0 / den)
    return post_mean + rng.standard_normal(post_mean.shape) * post_sd


def dense_solve(tex, alpha, beta, tau2, prior_mean, prior_std, rng, n_ens):
    """brews/baysparpy's target_timeseries_pred: an (N, N) solve + Cholesky per draw."""
    nd = len(tex)
    mu = np.ones(nd) * prior_mean
    inv_cov = np.eye(nd) * prior_std ** -2
    out = np.empty((nd, n_ens))
    for j in range(n_ens):
        inv_post_cov = inv_cov + beta[j] ** 2 / tau2[j] * np.eye(nd)
        post_cov = np.linalg.solve(inv_post_cov, np.eye(nd))
        sqrt_post_cov = np.linalg.cholesky(post_cov).T
        mean_first = inv_cov @ mu + (1 / tau2[j]) * beta[j] * (tex - alpha[j])
        out[:, j] = post_cov @ mean_first + sqrt_post_cov @ rng.standard_normal(nd)
    return out


def main() -> None:
    params = sio.loadmat(MO / f"Output_SpatAg_{RUNNAME}" / "params_standard.mat")
    obs = sio.loadmat(MO / f"obs{RUNNAME}.mat")
    ts = sio.loadmat(MO / "tex_testdata.mat")[SERIES][0, 0]

    tex = np.asarray(ts["tex86"]).ravel()
    lat = float(np.asarray(ts["lat"]).ravel()[0])
    lon = float(np.asarray(ts["lon"]).ravel()[0])
    print(f"{SERIES} ({RUNNAME}): N={tex.size}, site ({lon}, {lat}), mean TEX86 {tex.mean():.3f}")

    prior_mean, n_below = prior_mean_for(lon, lat, obs["locs_st_obs"],
                                         obs["st_obs_ave_vec"].ravel())
    print(f"prior mean {prior_mean:.3f} degC from {n_below} WOA cells within {MAX_DIST_KM:.0f} km")

    ind = matlab_thinning(params["tau2_samples"].size, N_DRAWS)
    tau2 = params["tau2_samples"].ravel()[ind]
    locs = params["Locs_Comp"].astype(float)
    cells = np.flatnonzero((np.abs(locs[:, 0] - lon) <= GRID_HALF_SPACE)
                           & (np.abs(locs[:, 1] - lat) <= GRID_HALF_SPACE))
    print(f"grid cell(s) matched: {locs[cells].tolist()}")
    alpha = params["alpha_samples_comp"][cells[0], ind]
    beta = params["beta_samples_comp"][cells[0], ind]

    rng = np.random.default_rng(0)
    t0 = time.perf_counter()
    preds = closed_form(tex, alpha, beta, tau2, prior_mean, PRIOR_STD, rng)
    t_cf = time.perf_counter() - t0
    # MATLAB prctile == Hazen plotting positions, NOT numpy's default.
    pct = np.percentile(preds, [5, 50, 95], axis=1, method="hazen").T
    print(f"\nclosed form: {t_cf:.4f} s; median spans "
          f"{pct[:, 1].min():.2f} to {pct[:, 1].max():.2f} degC")
    print("first three rows of 5/50/95:")
    print(np.round(pct[:3], 3))

    print(f"\nscaling, n_draws={N_DRAWS} (this machine):")
    print(f"  {'N':>5}  {'closed form':>12}  {'dense solve':>12}  {'speed-up':>9}")
    for nd in (50, 200, 500):
        fake = np.full(nd, tex.mean())
        t0 = time.perf_counter(); closed_form(fake, alpha, beta, tau2, prior_mean, PRIOR_STD, rng)
        t1 = time.perf_counter(); dense_solve(fake, alpha, beta, tau2, prior_mean, PRIOR_STD, rng, N_DRAWS)
        t2 = time.perf_counter()
        print(f"  {nd:5d}  {t1 - t0:11.4f}s  {t2 - t1:11.3f}s  {(t2 - t1) / (t1 - t0):8.0f}x")

    # How many cores does the dense path actually occupy?
    nd, n_ens = 300, 400
    fake = np.full(nd, tex.mean())
    w0, c0 = time.perf_counter(), time.process_time()
    dense_solve(fake, alpha, beta, tau2, prior_mean, PRIOR_STD, rng, n_ens)
    wall, cpu = time.perf_counter() - w0, time.process_time() - c0
    print(f"\ndense solve (N={nd}, {n_ens} draws): wall {wall:.2f}s, CPU {cpu:.2f}s "
          f"-> ~{cpu / wall:.1f} cores busy")
    print("(re-run with OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 to see it drop to ~1.0)")


if __name__ == "__main__":
    main()
