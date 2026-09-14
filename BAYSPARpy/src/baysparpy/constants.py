"""Constants hard-coded in the MATLAB source, gathered in one place.

Each names the MATLAB line it comes from, so a reviewer can check it against the
reference without searching.
"""

#: Radius of the earth in km. EarthChordDistances_2.m:10 (PORTING.md ECD-01).
EARTH_RADIUS_KM = 6378.137

#: Half-width of a calibration grid cell, in degrees. bayspar_tex.m:73,
#: TEX_forward.m:52 — "grid spacing is hard-coded here, will never change".
#: Not a parameter: the 162-cell tiling in the parameter file depends on it.
GRID_HALF_SPACE = 10.0

#: Search radius for the modern prior mean, in km. bayspar_tex.m:79 (BT-05).
MAX_DIST_KM = 500.0

#: Minimum number of observations averaged for the prior mean. bayspar_tex.m:76.
MIN_NUM = 1

#: Number of columns TEX_forward.m returns, whatever the ensemble size
#: (TEX_forward.m:92, PORTING.md TXF-07).
FORWARD_N_OUT = 1000

#: Default posterior draws. MATLAB's default; brews/baysparpy uses 5000.
DEFAULT_N_DRAWS = 1000

#: The two fitted models.
RUNNAMES = ("SST", "subT")

#: Lower-case aliases, for compatibility with brews/baysparpy's `temptype`.
RUNNAME_ALIASES = {"sst": "SST", "subt": "subT", "SST": "SST", "subT": "subT"}
