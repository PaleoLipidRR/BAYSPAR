"""Shared fixtures. The stores read ModelOutput/ from the MATLAB repo directly."""
from __future__ import annotations

import numpy as np
import pytest
import scipy.io as sio

from baysparpy.stores import find_modeloutput


@pytest.fixture(scope="session")
def modeloutput():
    """Path to the MATLAB ModelOutput directory."""
    return find_modeloutput()


@pytest.fixture(scope="session")
def lopes(modeloutput):
    """The lopes_santos2010 demo series: (tex, lon, lat)."""
    rec = sio.loadmat(modeloutput / "tex_testdata.mat")["lopes_santos2010"][0, 0]
    return (np.ravel(rec["tex86"]).astype(float),
            float(np.ravel(rec["lon"])[0]), float(np.ravel(rec["lat"])[0]))


@pytest.fixture(scope="session")
def wilsonlake(modeloutput):
    """The Wilson Lake PETM demo series (TEX86 only)."""
    rec = sio.loadmat(modeloutput / "wilsonlake.mat")["wilsonlake"][0, 0]
    return np.ravel(rec["tex86"]).astype(float)
