"""Audit of ModelOutput/ backing the factual claims in ../SPEC.md §3.

Run from the repository root:

    python docs/BAYSPARpy/audit/audit_modeloutput.py

Requires numpy and scipy only. Prints a report; exits non-zero if any claim the
spec relies on fails to hold.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import scipy.io as sio

ROOT = Path(__file__).resolve().parents[3]
MO = ROOT / "ModelOutput"
RUNS = ("SST", "subT")

failures: list[str] = []


def check(claim: str, ok: bool) -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {claim}")
    if not ok:
        failures.append(claim)


def section(title: str) -> None:
    print(f"\n{title}\n{'-' * len(title)}")


section("1. Parameter store inventory")
for run in RUNS:
    ps = sio.loadmat(MO / f"Output_SpatAg_{run}" / "params_standard.mat")
    pa = sio.loadmat(MO / f"Output_SpatAg_{run}" / "params_analog.mat")
    print(f"  {run}: standard alpha {ps['alpha_samples_comp'].shape}, "
          f"analog alpha {pa['alpha_samples'].shape}, "
          f"tau2 {ps['tau2_samples'].shape}, Locs_Comp {ps['Locs_Comp'].shape}")
    check(f"{run}: 20000 draws are shipped (docs say the cap is 15000)",
          ps["tau2_samples"].size == 20000)
    check(f"{run}: Locs_Comp is the full 18x9 tiling",
          ps["Locs_Comp"].shape == (162, 2)
          and sorted(np.unique(ps["Locs_Comp"][:, 0])) == list(range(-170, 180, 20))
          and sorted(np.unique(ps["Locs_Comp"][:, 1])) == list(range(-80, 90, 20)))

section("2. params_analog is a subset of params_standard (SPEC 3.2, finding 2)")
for run in RUNS:
    ps = sio.loadmat(MO / f"Output_SpatAg_{run}" / "params_standard.mat")
    pa = sio.loadmat(MO / f"Output_SpatAg_{run}" / "params_analog.mat")
    di = sio.loadmat(MO / f"Data_Input_SpatAg_{run}.mat")["Data_Input"][0, 0]
    locs_a = di["Locs"].astype(int)
    locs_c = ps["Locs_Comp"].astype(int)
    idx = np.array([np.flatnonzero((locs_c[:, 0] == lo) & (locs_c[:, 1] == la))[0]
                    for lo, la in locs_a])
    check(f"{run}: all {len(locs_a)} analog cells are in Locs_Comp", len(idx) == len(locs_a))
    check(f"{run}: alpha_samples == alpha_samples_comp[Data_Input.Locs]",
          np.array_equal(ps["alpha_samples_comp"][idx], pa["alpha_samples"]))
    check(f"{run}: beta_samples == beta_samples_comp[Data_Input.Locs]",
          np.array_equal(ps["beta_samples_comp"][idx], pa["beta_samples"]))
    check(f"{run}: tau2_samples identical in both files",
          np.array_equal(ps["tau2_samples"], pa["tau2_samples"]))

section("3. Data_Input structure (SPEC 3.3)")
for run in RUNS:
    di = sio.loadmat(MO / f"Data_Input_SpatAg_{run}.mat")["Data_Input"][0, 0]
    obs_all, obs_locs, target = di["ObsAll"], di["ObsLocs"], di["Target"]
    n_slot, n_cell = obs_all.shape
    n_obs = sum(np.size(obs_all[j, c]) for j in range(n_slot) for c in range(n_cell))
    per_slot = [sum(1 for c in range(n_cell) if np.size(obs_all[j, c])) for j in range(n_slot)]
    n_locs_each = {np.shape(obs_locs[j, c])[0]
                   for j in range(n_slot) for c in range(n_cell) if np.size(obs_locs[j, c])}
    print(f"  {run}: ObsAll {obs_all.shape} (site slot x cell), {n_obs} observations, "
          f"cells populated per slot {per_slot[:4]}...{per_slot[-1]}")
    check(f"{run}: ObsAll observations == len(Obs_Stack)", n_obs == di["Obs_Stack"].size)
    check(f"{run}: exactly one 1-degree location per (slot, cell)", n_locs_each == {1})
    check(f"{run}: slot 0 populated for every cell", per_slot[0] == n_cell)
    check(f"{run}: each Target slot carries Field and Err_Stds per cell",
          all(np.shape(target[0, j]["Field"]) == (n_cell, 1)
              and np.shape(target[0, j]["Err_Stds"]) == (n_cell, 1) for j in range(n_slot)))
    check(f"{run}: Inds_Stack indexes the {n_cell} cells",
          di["Inds_Stack"].min() == 1 and di["Inds_Stack"].max() == n_cell)

section("4. Store size, as float64 and float32 (SPEC 8)")
import io as _io  # noqa: E402

totals = {"float64": 0.0, "float32": 0.0}
for run in RUNS:
    ps = sio.loadmat(MO / f"Output_SpatAg_{run}" / "params_standard.mat")
    for name, dt in (("float64", np.float64), ("float32", np.float32)):
        buf = _io.BytesIO()
        np.savez_compressed(buf, **{k: (v.astype(dt) if v.dtype.kind == "f" else v)
                                    for k, v in ps.items() if not k.startswith("__")})
        mb = buf.getbuffer().nbytes / 1e6
        totals[name] += mb
        print(f"  {run} params_standard as {name}: {mb:6.1f} MB compressed")
print(f"  total without the redundant analog files: "
      f"{totals['float64']:.0f} MB float64 / {totals['float32']:.0f} MB float32 "
      f"(vs 145 MB of .mat today)")
a = sio.loadmat(MO / "Output_SpatAg_SST" / "params_standard.mat")["alpha_samples_comp"]
rel = np.abs(a.astype(np.float32) - a).max() / np.abs(a).max()
print(f"  worst float32 round-trip error on alpha: {rel:.2e} relative")

section("Result")
if failures:
    print(f"{len(failures)} claim(s) FAILED:")
    for f in failures:
        print(f"  - {f}")
    sys.exit(1)
print("all claims hold")
