"""Deterministic intermediates of the BAYSPAR demos, for review against MATLAB.

Everything printed here is computed *before* any random draw, so MATLAB and
Python must agree to the digit. Run the MATLAB snippet quoted under each block
and compare — a mismatch localises the porting error to one step.

    python docs/BAYSPARpy/audit/reference_trace.py > docs/BAYSPARpy/audit/reference_trace.txt
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import scipy.io as sio

ROOT = Path(__file__).resolve().parents[3]
MO = ROOT / "ModelOutput"

EARTH_RADIUS_KM = 6378.137
GRID_HALF_SPACE = 10.0
MAX_DIST_KM = 500.0
MIN_NUM = 1

np.set_printoptions(legacy="1.25")


def earth_chord_distances(point, points):
    lon0, lat0 = point
    d2r = np.pi / 180.0
    half = np.arcsin(np.sqrt(
        np.sin((lat0 - points[:, 1]) * d2r / 2) ** 2
        + np.cos(lat0 * d2r) * np.cos(points[:, 1] * d2r)
        * np.sin(np.abs(lon0 - points[:, 0]) * d2r / 2) ** 2))
    return 2 * EARTH_RADIUS_KM * np.sin(half)


def matlab_thinning(n_total: int, n_draws: int) -> np.ndarray:
    """MATLAB ind_s = round(linspace(1, n_total, n_draws)); 0-based on return."""
    return np.floor(np.linspace(1, n_total, n_draws) + 0.5).astype(int) - 1


def hr(title: str) -> None:
    print(f"\n{'=' * 78}\n{title}\n{'=' * 78}")


hr("A. Thinning index  --  MATLAB: ind_s = round(linspace(1, 20000, Nsamps))")
print("  Reported 1-based, as MATLAB prints them.")
for n_draws in (10, 999, 1000, 1001, 20000):
    ind = matlab_thinning(20000, n_draws) + 1
    print(f"  Nsamps={n_draws:>5}: first 5 {ind[:5].tolist()}  last 3 {ind[-3:].tolist()}  "
          f"unique={len(np.unique(ind))}")

hr("B. Chordal distance  --  MATLAB: EarthChordDistances_2([lon lat], [lon2 lat2])")
pairs = [((0.0, 0.0), (0.0, 1.0)), ((0.0, 0.0), (1.0, 0.0)),
         ((-17.6635, 9.166), (-17.5, 9.5)), ((100.0, -45.0), (-100.0, 45.0)),
         ((180.0, 0.0), (-180.0, 0.0))]
for p, q in pairs:
    d = earth_chord_distances(p, np.array([q], dtype=float))[0]
    print(f"  ({p[0]:>9.4f}, {p[1]:>8.4f}) -> ({q[0]:>9.4f}, {q[1]:>8.4f})  =  {d:.9f} km")

hr("C. Prior mean and grid cell, per demo series")
print("  MATLAB (inside bayspar_tex.m, after line 126 and line 140):")
print("      Output_Struct.PriorMean, Output_Struct.GridLoc\n")
series = sio.loadmat(MO / "tex_testdata.mat")
for run in ("SST", "subT"):
    obs = sio.loadmat(MO / f"obs{run}.mat")
    locs_obs, st_obs = obs["locs_st_obs"], obs["st_obs_ave_vec"].ravel()
    locs_comp = sio.loadmat(MO / f"Output_SpatAg_{run}" / "params_standard.mat")["Locs_Comp"].astype(float)
    print(f"  --- runname = '{run}' ---")
    for name in ("castaneda2010", "lopes_santos2010", "shevenell2011"):
        rec = series[name][0, 0]
        lat = float(np.asarray(rec["lat"]).ravel()[0])
        lon = float(np.asarray(rec["lon"]).ravel()[0])
        tex = np.asarray(rec["tex86"]).ravel()
        d = earth_chord_distances((lon, lat), locs_obs)
        order = np.argsort(d, kind="stable")
        n_below = int(np.searchsorted(d[order], MAX_DIST_KM, side="left"))
        take = n_below if n_below > MIN_NUM else MIN_NUM
        prior_mean = st_obs[order[:take]].mean()
        cells = np.flatnonzero((np.abs(locs_comp[:, 0] - lon) <= GRID_HALF_SPACE)
                               & (np.abs(locs_comp[:, 1] - lat) <= GRID_HALF_SPACE))
        print(f"  {name:<18} N={tex.size:>4}  site=({lon:>9.4f}, {lat:>8.4f})  "
              f"mean TEX86={tex.mean():.6f}")
        print(f"  {'':<18} PriorMean={prior_mean:.10f}  from {n_below} obs within "
              f"{MAX_DIST_KM:.0f} km")
        print(f"  {'':<18} GridLoc={locs_comp[cells].astype(int).tolist()}  "
              f"(cells matched: {cells.size})")

hr("D. Analogue selection  --  Demo_AnalogPrediction.m, wilsonlake, runname='SST'")
print("  MATLAB: search_tol = std(dats)*2;  then bayspar_tex_analog(...);")
print("          Output_Struct.AnLocs\n")
wl = sio.loadmat(MO / "wilsonlake.mat")["wilsonlake"][0, 0]
tex = np.asarray(wl["tex86"]).ravel()
di = sio.loadmat(MO / "Data_Input_SpatAg_SST.mat")["Data_Input"][0, 0]
locs = di["Locs"].astype(int)
inds = np.ravel(di["Inds_Stack"])
obs_stack = np.ravel(di["Obs_Stack"])
spatial_mean = np.array([obs_stack[inds == i + 1].mean() for i in range(len(locs))])
tol = tex.std(ddof=1) * 2          # MATLAB std() normalises by N-1
print(f"  N={tex.size}  mean TEX86={tex.mean():.10f}  std(dats)={tex.std(ddof=1):.10f}  "
      f"search_tol={tol:.10f}")
sel = np.flatnonzero((spatial_mean >= tex.mean() - tol) & (spatial_mean <= tex.mean() + tol))
print(f"  cell-mean TEX86 range across the 80 calibration cells: "
      f"{spatial_mean.min():.6f} to {spatial_mean.max():.6f}")
print(f"  analogues selected: {sel.size}")
for i in sel:
    print(f"    cell {i + 1:>3}  centroid [{locs[i,0]:>5}, {locs[i,1]:>4}]  "
          f"mean TEX86 {spatial_mean[i]:.10f}")

hr("E. Per-draw posterior, before noise  --  bayspar_tex.m lines 149-152")
print("  MATLAB: post_mean and post_sig for dats(1), draws 1..5 of the thinned set.")
print("  Demo_StandardPrediction.m as shipped: lopes_santos2010, prior_std=6, runname='subT'\n")
rec = series["lopes_santos2010"][0, 0]
tex = np.asarray(rec["tex86"]).ravel()
lat = float(np.asarray(rec["lat"]).ravel()[0])
lon = float(np.asarray(rec["lon"]).ravel()[0])
params = sio.loadmat(MO / "Output_SpatAg_subT" / "params_standard.mat")
obs = sio.loadmat(MO / "obssubT.mat")
d = earth_chord_distances((lon, lat), obs["locs_st_obs"])
order = np.argsort(d, kind="stable")
n_below = int(np.searchsorted(d[order], MAX_DIST_KM, side="left"))
prior_mean = obs["st_obs_ave_vec"].ravel()[order[:n_below]].mean()
prior_std = 6.0
ind = matlab_thinning(params["tau2_samples"].size, 1000)
locs_comp = params["Locs_Comp"].astype(float)
cell = np.flatnonzero((np.abs(locs_comp[:, 0] - lon) <= GRID_HALF_SPACE)
                      & (np.abs(locs_comp[:, 1] - lat) <= GRID_HALF_SPACE))[0]
alpha = params["alpha_samples_comp"][cell, ind]
beta = params["beta_samples_comp"][cell, ind]
tau2 = params["tau2_samples"].ravel()[ind]
pinv = prior_std ** -2
den = pinv + beta ** 2 / tau2
num = pinv * prior_mean + (beta / tau2) * (tex[0] - alpha)
print(f"  dats(1)={tex[0]:.10f}  PriorMean={prior_mean:.10f}  prior_std={prior_std}")
print(f"  {'draw':>6} {'alpha':>14} {'beta':>14} {'tau2':>14} {'post_mean':>14} {'post_sig':>12}")
for j in range(5):
    print(f"  {ind[j] + 1:>6} {alpha[j]:>14.10f} {beta[j]:>14.10f} {tau2[j]:>14.10f} "
          f"{num[j] / den[j]:>14.10f} {np.sqrt(1 / den[j]):>12.10f}")

hr("F. Percentile convention  --  MATLAB: prctile(x, [5 50 95])")
print("  MATLAB uses Hazen plotting positions (i-0.5)/n. NumPy's default does not.\n")
x = np.arange(1.0, 11.0)
for method in ("hazen", "linear", "midpoint"):
    q = np.percentile(x, [5, 50, 95], method=method)
    mark = "  <- matches MATLAB" if method == "hazen" else ""
    print(f"  np.percentile(1..10, [5,50,95], method={method!r:<10}) = "
          f"{np.round(q, 6).tolist()}{mark}")
print("\n  MATLAB check:  prctile(1:10, [5 50 95])   % expect 1.0  5.5  10.0")
