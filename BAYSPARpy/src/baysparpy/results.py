"""Result containers — the Python form of MATLAB's ``Output_Struct``."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass
class Prediction:
    """Temperatures predicted from TEX86 (PORTING.md BT-12).

    Attributes:
        preds: (N, 3) the 5th, 50th and 95th percentiles, °C.
        prior_mean: Prior mean used, °C.
        prior_std: Prior standard deviation used, °C.
        site_loc: (lon, lat) of the site, standard mode only.
        grid_loc: (lon, lat) of the grid cell used, standard mode only.
        analog_locs: (n_analogues, 2) centroids of the analogue cells,
            analogue mode only. Aligned with axis 1 of :attr:`ensemble`.
        ensemble: The draws, if `save_ensemble` was set — (N, n_draws) in
            standard mode, (N, n_analogues, n_draws) in analogue mode
            (PORTING.md BTA-08). ``None`` otherwise.
        metadata: Provenance — runname, n_draws, mode, seed, thinning indices.
    """

    preds: np.ndarray
    prior_mean: float
    prior_std: float
    site_loc: tuple[float, float] | None = None
    grid_loc: tuple[float, float] | None = None
    analog_locs: np.ndarray | None = None
    ensemble: np.ndarray | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def p5(self) -> np.ndarray:
        """The 5th percentile, °C."""
        return self.preds[:, 0]

    @property
    def p50(self) -> np.ndarray:
        """The median, °C."""
        return self.preds[:, 1]

    @property
    def p95(self) -> np.ndarray:
        """The 95th percentile, °C."""
        return self.preds[:, 2]

    def to_dict(self) -> dict[str, Any]:
        """The MATLAB field names, for diffing against a saved ``Output_Struct``.

        Returns:
            A dict keyed as MATLAB's structure is.
        """
        out: dict[str, Any] = {"Preds": self.preds, "PriorMean": self.prior_mean,
                               "PriorStd": self.prior_std}
        if self.site_loc is not None:
            out["SiteLoc"] = np.asarray(self.site_loc)
        if self.grid_loc is not None:
            out["GridLoc"] = np.asarray(self.grid_loc)
        if self.analog_locs is not None:
            out["AnLocs"] = self.analog_locs
        if self.ensemble is not None:
            out["PredsEns"] = self.ensemble
        return out

    def to_dataframe(self):
        """The percentiles as a DataFrame.

        Returns:
            A pandas DataFrame with columns ``p5``, ``p50``, ``p95``.

        Raises:
            ImportError: if pandas is not installed.
        """
        try:
            import pandas as pd
        except ImportError as exc:
            raise ImportError("to_dataframe() needs pandas: pip install pandas") from exc
        return pd.DataFrame(self.preds, columns=["p5", "p50", "p95"])


@dataclass
class ForwardPrediction:
    """TEX86 forward-modelled from temperature (PORTING.md TXF-07).

    Attributes:
        tex86: (N, n_out) simulated TEX86, clipped to [0, 1].
        metadata: Provenance, including ``clipped_fraction`` — the share of the
            ensemble pinned at a bound, which for cold subT draws is not small.
    """

    tex86: np.ndarray
    metadata: dict[str, Any] = field(default_factory=dict)

    def percentiles(self, q=(5, 50, 95)) -> np.ndarray:
        """Percentiles across the ensemble, MATLAB's convention.

        Args:
            q: Percentiles in 0-100.

        Returns:
            (N, len(q)) array.
        """
        from .utils import prctile
        return prctile(self.tex86, list(q), axis=1)
