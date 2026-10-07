"""Compute every number the dossier charts, from the three implementations.

Run from BAYSPARpy/ with `uv run python ../docs/BAYSPARpy/audit/dossier/gen_data.py`.
Writes data.json next to this file. Needs brews/baysparpy installed (`uv sync` does that).
"""
import base64, json, time
from pathlib import Path

import numpy as np
import scipy.io as sio

import baysparpy as ours
import bayspar as brews
from bayspar.observations.core import chord_distance as brews_chord

ROOT = Path(__file__).resolve().parents[4]          # the BAYSPAR repository
OUT = Path(__file__).with_name("data.json")
REF = sio.loadmat(str(ROOT / "docs/BAYSPARpy/audit/original_matlab.mat"),
                  squeeze_me=True, struct_as_record=False)
N = 1000
SEED = 0
r = lambda a, k=4: np.round(np.asarray(a, dtype=float), k).tolist()
D = {"meta": {"engine": str(REF["engine"]), "generated": str(REF["generated"]),
              "numpy": np.__version__, "ours": ours.__version__}}

# ---------------------------------------------------------------- standard mode
si = REF["standard_input"]
tex = np.atleast_1d(si.tex86).astype(float)
lon, lat, pstd, run = float(si.lon), float(si.lat), float(si.prior_std), str(si.runname)
m_pred = np.atleast_2d(REF["standard"].Preds)
o_out = ours.bayspar_tex(tex, lon, lat, pstd, run, n_draws=N, seed=SEED)
np.random.seed(SEED)
b_out = brews.predict_seatemp(tex, lat=lat, lon=lon, prior_std=pstd, temptype="sst", nens=N)
b_pred = np.percentile(b_out.ensemble, [5, 50, 95], axis=1, method="nearest").T
mo = sio.loadmat(str(ROOT / "ModelOutput/tex_testdata.mat"))
age = np.ravel(mo["lopes_santos2010"][0, 0]["age"]).astype(float)
o = np.argsort(age)
D["standard"] = {
    "site": [lon, lat], "prior_std": pstd, "runname": run, "n": int(tex.size),
    "grid": r(REF["standard"].GridLoc, 1), "grid_ours": list(o_out.grid_loc),
    "prior_mean": {"MATLAB": float(REF["standard"].PriorMean), "BAYSPARpy": o_out.prior_mean,
                   "brews": float(b_out.prior_mean)},
    "age": r(age[o], 2), "tex": r(tex[o]),
    "MATLAB": r(m_pred[o], 3), "BAYSPARpy": r(o_out.preds[o], 3), "brews": r(b_pred[o], 3),
    "summary": {k: r(np.median(p, axis=0), 4) for k, p in
                [("MATLAB", m_pred), ("BAYSPARpy", o_out.preds), ("brews", b_pred)]},
    "maxdiff": {"BAYSPARpy": float(np.abs(o_out.preds - m_pred).max()),
                "brews": float(np.abs(b_pred - m_pred).max())},
    "meandiff50": {"BAYSPARpy": float(np.abs(o_out.preds - m_pred)[:, 1].mean()),
                   "brews": float(np.abs(b_pred - m_pred)[:, 1].mean())},
}
# Monte-Carlo floor: ours vs ours at another seed
o2 = ours.bayspar_tex(tex, lon, lat, pstd, run, n_draws=N, seed=1)
D["standard"]["mc_floor50"] = float(np.abs(o2.preds - o_out.preds)[:, 1].mean())

# ---------------------------------------------------------------- prior-mean search
obs = ours.get_seatemp(run)
d = ours.earth_chord_distances((lon, lat), obs.locs)[0]
box = (obs.locs[:, 0] >= -32) & (obs.locs[:, 0] <= 2) & (obs.locs[:, 1] >= -12) & (obs.locs[:, 1] <= 22)
D["prior"] = {"pts": [[float(a), float(b), round(float(v), 2), round(float(dd), 1)]
                      for (a, b), v, dd in zip(obs.locs[box], obs.values[box], d[box])],
              "n_below": int((d < 500).sum()), "max_dist": 500.0}

# ---------------------------------------------------------------- ocean mask (1° cells with SST obs)
mask = np.zeros((180, 360), dtype=np.uint8)
L = obs.locs
mask[(89.5 - L[:, 1]).astype(int), (L[:, 0] + 179.5).astype(int)] = 1
D["ocean"] = base64.b64encode(np.packbits(mask.ravel()).tobytes()).decode()

# ---------------------------------------------------------------- draws, thinning
dr = ours.get_draws(run)
cell = dr.cell_index(lon, lat)
spread = ours.matlab_thinning(dr.n_draws, N)
first = np.arange(N)
brewsidx = np.arange(10000, 10000 + N)   # brews ships draws 10001-20000 and takes the first nens
step = 10
D["draws"] = {
    "n_total": int(dr.n_draws), "cell": cell, "cell_loc": dr.locs[cell].tolist(),
    "alpha_trace": r(dr.alpha[cell, ::step], 4), "step": step,
    "spread_head": (spread[:6] + 1).tolist(), "spread_tail": (spread[-3:] + 1).tolist(),
    "alpha_mean": {"all": float(dr.alpha[cell].mean()), "spread": float(dr.alpha[cell, spread].mean()),
                   "first": float(dr.alpha[cell, first].mean()),
                   "brews": float(dr.alpha[cell, brewsidx].mean())},
    "median": {"alpha": float(np.median(dr.alpha[cell])), "beta": float(np.median(dr.beta[cell])),
               "tau2": float(np.median(dr.tau2))},
    "locs": dr.locs.astype(int).tolist(),
}

def preds_for(index, method, seed=SEED):
    rng = np.random.default_rng(seed)
    pm, ps = ours.solve_posterior(tex, dr.alpha[cell, index], dr.beta[cell, index],
                                  dr.tau2[index], o_out.prior_mean, pstd)
    ens = pm + rng.standard_normal(pm.shape) * ps
    return np.percentile(ens, [5, 50, 95], axis=1, method=method).T

base = preds_for(spread, "hazen")
eff = lambda p: [float(np.abs(p[:, 1] - base[:, 1]).mean()), float(np.abs(p[:, 1] - base[:, 1]).max())]
D["effects"] = {"first_draws": eff(preds_for(first, "hazen")),
                "brews_draws": eff(preds_for(brewsidx, "hazen")),
                "brews_both": eff(preds_for(brewsidx, "nearest")),
                "nearest": eff(preds_for(spread, "nearest")),
                "both": eff(preds_for(first, "nearest")),
                "mc": eff(preds_for(spread, "hazen", seed=1))}

# ---------------------------------------------------------------- analogue mode
ai, aref = REF["analog_input"], REF["analog"]
wl = np.atleast_1d(ai.tex86).astype(float)
tol, pm_a, ps_a = float(ai.search_tol), float(ai.prior_mean), float(ai.prior_std)
m_an = np.atleast_2d(aref.Preds)
o_an = ours.bayspar_tex_analog(wl, pm_a, ps_a, tol, "SST", n_draws=N, mode="matlab", seed=SEED)
o_an_mod = ours.bayspar_tex_analog(wl, pm_a, ps_a, tol, "SST", n_draws=N, mode="modern", seed=SEED)
o_an_mod2 = ours.bayspar_tex_analog(wl, pm_a, ps_a, tol, "SST", n_draws=N, mode="modern", seed=1)
np.random.seed(SEED)
b_an = brews.predict_seatemp_analog(wl, prior_std=ps_a, temptype="sst", search_tol=tol,
                                    prior_mean=pm_a, nens=N, progressbar=False)
b_an_p = np.percentile(b_an.ensemble, [5, 50, 95], axis=(1, 2), method="nearest").T
ct = ours.get_coretops("SST")
cm = ct.cell_means("tex86")
sel = set(ours.find_analogs(cm, wl.mean(), tol).tolist())
depth = np.ravel(sio.loadmat(str(ROOT / "ModelOutput/wilsonlake.mat"))["wilsonlake"][0, 0]["depth"]).astype(float)
oa = np.argsort(depth)
D["analog"] = {
    "tex_mean": float(wl.mean()), "tol": tol, "prior_mean": pm_a, "prior_std": ps_a, "n": int(wl.size),
    "cells": [[float(a), float(b), round(float(v), 4), int(i in sel)]
              for i, ((a, b), v) in enumerate(zip(ct.locs, cm))],
    "n_sel": len(sel), "depth": r(depth[oa], 2),
    "MATLAB": r(m_an[oa], 3), "BAYSPARpy": r(o_an.preds[oa], 3), "brews": r(b_an_p[oa], 3),
    "summary": {k: r(np.median(p, axis=0), 4) for k, p in
                [("MATLAB", m_an), ("BAYSPARpy", o_an.preds), ("brews", b_an_p)]},
    "pairing_diff": r(np.abs(o_an.preds - o_an_mod.preds).mean(axis=0), 4),
    "pairing_mc": r(np.abs(o_an_mod2.preds - o_an_mod.preds).mean(axis=0), 4),
    "ens_shape": {"MATLAB": [int(x) for x in aref.PredsEns_size],
                  "BAYSPARpy": list(o_an_mod.ensemble.shape) if o_an_mod.ensemble is not None else None,
                  "brews": list(b_an.ensemble.shape)},
}
o_an_full = ours.bayspar_tex_analog(wl, pm_a, ps_a, tol, "SST", n_draws=N, save_ensemble=True, seed=SEED)
o_an_mat = ours.bayspar_tex_analog(wl, pm_a, ps_a, tol, "SST", n_draws=N, save_ensemble=True,
                                   mode="matlab", seed=SEED)
D["analog"]["ens_shape"]["BAYSPARpy"] = list(o_an_full.ensemble.shape)
D["analog"]["ens_shape"]["BAYSPARpy_matlab"] = list(o_an_mat.ensemble.shape)

# ---------------------------------------------------------------- forward model
fi = REF["forward_input"]
t_in = np.atleast_1d(fi.t).astype(float)
m_fw = np.asarray(REF["forward"], dtype=float)
o_fw = ours.tex_forward(np.repeat(float(fi.lat), 3), np.repeat(float(fi.lon), 3), t_in, "SST", seed=SEED)
np.random.seed(SEED)
b_fw = brews.predict_tex(t_in, lat=float(fi.lat), lon=float(fi.lon), temptype="sst", nens=N).ensemble
bins = np.linspace(0.3, 1.0, 36)
D["forward"] = {
    "t": t_in.tolist(), "bins": r(bins, 3),
    "hist": {k: [np.histogram(a[i], bins=bins)[0].tolist() for i in range(3)]
             for k, a in [("MATLAB", m_fw), ("BAYSPARpy", o_fw.tex86), ("brews", b_fw)]},
    "stats": {k: [[float(np.mean(a[i])), float(np.std(a[i])), float(np.min(a[i])), float(np.max(a[i]))]
                  for i in range(3)]
              for k, a in [("MATLAB", m_fw), ("BAYSPARpy", o_fw.tex86), ("brews", b_fw)]},
    "shape": {"MATLAB": list(m_fw.shape), "BAYSPARpy": list(o_fw.tex86.shape), "brews": list(b_fw.shape)},
}
# forward curve over temperature, our port, plus clipping fraction for subT at the same site
tg = np.arange(-2, 36.1, 1.0)
fc = ours.tex_forward(lat, lon, tg, "SST", seed=SEED, mode="modern")
pre = None
D["forward"]["curve"] = {"t": tg.tolist(),
                         "p": r(np.percentile(fc.tex86, [5, 50, 95], axis=1).T, 4),
                         "clipped": fc.metadata["clipped_fraction"]}
hot = ours.tex_forward(lat, lon, [34.0, 36.0, 38.0], "SST", seed=SEED)
D["forward"]["clip_hot"] = hot.metadata["clipped_fraction"]

# ---------------------------------------------------------------- distances
pts = obs.locs[::997]
d_ours = ours.earth_chord_distances((lon, lat), pts)[0]
d_brews = np.ravel(brews_chord(np.array([[lat, lon]]), pts[:, ::-1]))
D["distance"] = {"max_abs_diff_vs_brews": float(np.abs(d_ours - d_brews).max()),
                 "n_pairs": int(pts.shape[0])}

# ---------------------------------------------------------------- timing
def timed(fn, *a, **k):
    best = 1e9
    for _ in range(2):
        t0 = time.perf_counter(); fn(*a, **k); best = min(best, time.perf_counter() - t0)
    return best

Ns = [25, 50, 100, 193, 400]
rng = np.random.default_rng(3)
timing = {"N": Ns, "BAYSPARpy": [], "brews": []}
for n in Ns:
    x = rng.choice(tex, size=n)
    timing["BAYSPARpy"].append(timed(ours.bayspar_tex, x, lon, lat, pstd, run, n_draws=N, seed=0))
    timing["brews"].append(timed(brews.predict_seatemp, x, lat=lat, lon=lon, prior_std=pstd,
                                 temptype="sst", nens=N))
D["timing"] = timing

OUT.write_text(json.dumps(D, separators=(",", ":")))
print("wrote", OUT, OUT.stat().st_size, "bytes")
print(json.dumps({k: D["standard"][k] for k in ["prior_mean", "summary", "maxdiff", "meandiff50", "mc_floor50", "grid", "grid_ours"]}, indent=1))
print(json.dumps({k: D["analog"][k] for k in ["summary", "n_sel", "pairing_diff", "pairing_mc", "ens_shape"]}))
print(json.dumps(D["forward"]["stats"])); print(D["forward"]["shape"], D["forward"]["curve"]["clipped"], D["forward"]["clip_hot"])
print(D["effects"], D["draws"]["alpha_mean"], D["draws"]["median"], D["prior"]["n_below"], len(D["prior"]["pts"]))
print(D["distance"], D["timing"])
