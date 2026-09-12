"""Reading the fitted parameters and the calibration inputs.

The MATLAB functions ``load`` their ``.mat`` files on every call, resolved
against the current working directory, which is why the demos only run from the
repository root. These stores resolve the directory once and cache each file for
the life of the process (PORTING.md BT-02).

Set ``BAYSPAR_MODELOUTPUT`` to point at a ``ModelOutput`` directory elsewhere.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
import scipy.io as sio

from .constants import GRID_HALF_SPACE
from .distance import earth_chord_distances
from .errors import AmbiguousGridCellError
from .utils import check_runname, matlab_thinning


def find_modeloutput() -> Path:
    """Locate the ``ModelOutput`` directory.

    Search order: ``BAYSPAR_MODELOUTPUT``, then a ``ModelOutput`` directory in
    any parent of this file (which finds it in a checkout of the MATLAB repo).

    Returns:
        Path to the directory.

    Raises:
        FileNotFoundError: if it is not found.
    """
    env = os.environ.get("BAYSPAR_MODELOUTPUT")
    if env:
        p = Path(env).expanduser()
        if (p / "obsSST.mat").is_file():
            return p
        raise FileNotFoundError(
            f"BAYSPAR_MODELOUTPUT is set to {p}, which has no obsSST.mat in it"
        )
    for parent in Path(__file__).resolve().parents:
        candidate = parent / "ModelOutput"
        if (candidate / "obsSST.mat").is_file():
            return candidate
    raise FileNotFoundError(
        "no ModelOutput directory found above "
        f"{Path(__file__).resolve()}; set BAYSPAR_MODELOUTPUT to point at one"
    )


@lru_cache(maxsize=8)
def _load_mat(path_str: str) -> dict:
    return sio.loadmat(path_str)


@dataclass(frozen=True)
class Draws:
    """Posterior draws of the regression parameters, on the 20° × 20° grid.

    Attributes:
        runname: ``"SST"`` or ``"subT"``.
        alpha: (162, n_draws) intercept draws, one row per grid cell.
        beta: (162, n_draws) slope draws.
        tau2: (n_draws,) residual variance draws, shared across cells.
        locs: (162, 2) cell centroids as (lon, lat), float.
    """

    runname: str
    alpha: np.ndarray
    beta: np.ndarray
    tau2: np.ndarray
    locs: np.ndarray

    @property
    def n_draws(self) -> int:
        """Draws in the store (20,000 as shipped, not the documented 15,000)."""
        return self.tau2.size

    def thin(self, n_draws: int) -> "Draws":
        """Take `n_draws` draws spread across the chain (PORTING.md BT-06).

        Args:
            n_draws: How many to keep.

        Returns:
            A new :class:`Draws` holding the thinned columns.
        """
        ind = matlab_thinning(self.n_draws, n_draws)
        return Draws(self.runname, self.alpha[:, ind], self.beta[:, ind],
                     self.tau2[ind], self.locs)

    def cell_index(self, lon: float, lat: float, *, strict: bool = False) -> int:
        """Index of the grid cell containing a site (PORTING.md BT-08).

        The MATLAB test is inclusive at ±10°, so a site exactly on a boundary
        matches two cells (a corner, four). MATLAB then errors for N > 1 and
        silently returns one row per cell for N = 1; here the nearest centroid
        wins and a warning names the alternatives.

        Args:
            lon: Longitude, -180 to 180.
            lat: Latitude, -90 to 90.
            strict: Raise on a boundary site instead of choosing.

        Returns:
            Row index into :attr:`locs`.

        Raises:
            AmbiguousGridCellError: on a boundary site when `strict`.
            ValueError: if no cell contains the site.
        """
        import warnings

        hits = np.flatnonzero((np.abs(self.locs[:, 0] - lon) <= GRID_HALF_SPACE)
                              & (np.abs(self.locs[:, 1] - lat) <= GRID_HALF_SPACE))
        if hits.size == 0:
            raise ValueError(
                f"({lon}, {lat}) is outside the calibration grid; longitude must be "
                "in [-180, 180] and latitude in [-90, 90]"
            )
        if hits.size == 1:
            return int(hits[0])
        d = earth_chord_distances((lon, lat), self.locs[hits])[0]
        nearest = int(hits[int(np.argmin(d))])
        matched = self.locs[hits].astype(int).tolist()
        if strict:
            raise AmbiguousGridCellError(
                f"({lon}, {lat}) lies on a cell boundary and matches {matched}"
            )
        warnings.warn(
            f"({lon}, {lat}) lies on a cell boundary and matches {matched}; using the "
            f"nearest centroid {self.locs[nearest].astype(int).tolist()}. Pass "
            "strict_grid=True to raise instead.",
            UserWarning, stacklevel=3,
        )
        return nearest


@dataclass(frozen=True)
class SeaTempObs:
    """The 1° instrumental climatology used for the prior mean.

    Attributes:
        runname: ``"SST"`` or ``"subT"`` — the two have different grids, because
            subT drops cells where the 0-200 m average is undefined.
        locs: (n, 2) locations as (lon, lat).
        values: (n,) temperatures in °C.
    """

    runname: str
    locs: np.ndarray
    values: np.ndarray

    def prior_mean(self, lon: float, lat: float, *, max_dist: float,
                   min_num: int) -> tuple[float, int]:
        """Mean temperature near a site (PORTING.md BT-07).

        All observations strictly within `max_dist`, or the closest `min_num`,
        whichever is the larger set — MATLAB's ``find(vals < max_dist, 1, 'last')``
        on sorted distances is a count, and its empty result falls through to the
        closest point.

        Args:
            lon: Longitude of the site.
            lat: Latitude of the site.
            max_dist: Search radius in km.
            min_num: Fallback count when nothing is within `max_dist`.

        Returns:
            The prior mean in °C, and how many observations were within
            `max_dist`.
        """
        d = earth_chord_distances((lon, lat), self.locs)[0]
        # stable, so tied distances keep their original order as MATLAB's sort does
        order = np.argsort(d, kind="stable")
        n_below = int(np.searchsorted(d[order], max_dist, side="left"))
        take = n_below if n_below > min_num else min_num
        return float(self.values[order[:take]].mean()), n_below


@dataclass(frozen=True)
class CoreTops:
    """The coretop calibration inputs, from ``Data_Input_SpatAg_*.mat``.

    Attributes:
        runname: ``"SST"`` or ``"subT"``.
        locs: (80, 2) centroids of the cells that hold coretop data.
        cell_index: (n_obs,) 0-based index into :attr:`locs` per observation.
        tex86: (n_obs,) observed TEX86.
        target_t: (n_obs,) target temperature at each observation's site, °C.
        target_t_sd: (n_obs,) the target's error standard deviation.
    """

    runname: str
    locs: np.ndarray
    cell_index: np.ndarray
    tex86: np.ndarray
    target_t: np.ndarray
    target_t_sd: np.ndarray

    def cell_means(self, field: str) -> np.ndarray:
        """Mean of a field within each of the 80 cells (PORTING.md BTA-02).

        Note this averages *observations*, not sites: a cell with eight cores at
        one 1° location and one at another weights the first eight times. That is
        MATLAB's behaviour.

        Args:
            field: ``"tex86"`` (as ``bayspar_tex_analog``) or ``"target_t"``
                (as ``TEX_forward``) — see PORTING.md TXF-04.

        Returns:
            (80,) array of cell means.
        """
        values = getattr(self, field)
        n = len(self.locs)
        counts = np.bincount(self.cell_index, minlength=n)
        sums = np.bincount(self.cell_index, weights=values, minlength=n)
        return sums / counts   # every cell holds >= 1 observation (audited)


def get_draws(runname: str) -> Draws:
    """Load the fitted parameters for a model.

    Args:
        runname: ``"SST"`` or ``"subT"``.

    Returns:
        The full 20,000-draw :class:`Draws`.
    """
    runname = check_runname(runname)
    mat = _load_mat(str(find_modeloutput() / f"Output_SpatAg_{runname}" / "params_standard.mat"))
    return Draws(runname,
                 np.asarray(mat["alpha_samples_comp"], dtype=float),
                 np.asarray(mat["beta_samples_comp"], dtype=float),
                 np.asarray(mat["tau2_samples"], dtype=float).ravel(),
                 np.asarray(mat["Locs_Comp"], dtype=float))


def get_seatemp(runname: str) -> SeaTempObs:
    """Load the instrumental climatology for a model.

    Args:
        runname: ``"SST"`` or ``"subT"``.

    Returns:
        The :class:`SeaTempObs` for that model.
    """
    runname = check_runname(runname)
    stem = "obsSST" if runname == "SST" else "obssubT"
    mat = _load_mat(str(find_modeloutput() / f"{stem}.mat"))
    return SeaTempObs(runname,
                      np.asarray(mat["locs_st_obs"], dtype=float),
                      np.asarray(mat["st_obs_ave_vec"], dtype=float).ravel())


def get_coretops(runname: str) -> CoreTops:
    """Load the coretop calibration inputs for a model.

    Args:
        runname: ``"SST"`` or ``"subT"``.

    Returns:
        The :class:`CoreTops` for that model.
    """
    runname = check_runname(runname)
    mat = _load_mat(str(find_modeloutput() / f"Data_Input_SpatAg_{runname}.mat"))
    di = mat["Data_Input"][0, 0]
    return CoreTops(
        runname,
        np.asarray(di["Locs"], dtype=float),
        np.asarray(di["Inds_Stack"], dtype=int).ravel() - 1,   # 1-based on disk
        np.asarray(di["Obs_Stack"], dtype=float).ravel(),
        np.asarray(di["Target_Stack"], dtype=float).ravel(),
        np.asarray(di["Target_Err_Stds_Stack"], dtype=float).ravel(),
    )


def analog_cell_rows(draws: Draws, coretops: CoreTops) -> np.ndarray:
    """Rows of `draws` corresponding to the coretop cells (PORTING.md STO-01).

    ``params_analog.mat`` is bit-identical to these rows of
    ``params_standard.mat``, in ``Data_Input.Locs`` order, so analogue mode reads
    the standard store by index rather than loading a second 24 MB file.

    Args:
        draws: The full parameter store.
        coretops: The calibration inputs, whose `locs` give the order.

    Returns:
        (80,) row indices into ``draws.alpha`` / ``draws.beta``.
    """
    rows = []
    for lon, lat in coretops.locs:
        hit = np.flatnonzero((draws.locs[:, 0] == lon) & (draws.locs[:, 1] == lat))
        if hit.size != 1:
            raise ValueError(f"coretop cell ({lon}, {lat}) is not a grid centroid")
        rows.append(int(hit[0]))
    return np.asarray(rows, dtype=int)
