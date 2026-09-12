"""Chordal distances on a sphere — a port of EarthChordDistances_2.m."""
from __future__ import annotations

import numpy as np

from .constants import EARTH_RADIUS_KM


def earth_chord_distances(points1, points2) -> np.ndarray:
    """Chordal distance in km between two sets of (lon, lat) points.

    A port of ``EarthChordDistances_2.m``, term for term (PORTING.md ECD-01..05).
    This is the straight-line chord through the sphere, ``2R sin(a/2)``, not the
    great-circle arc ``R a`` — 0.026% shorter at 500 km.

    MATLAB builds an N·M × 4 matrix of every pair; this broadcasts instead, which
    matters because the prior-mean search pairs one site against ~37,000 cells.

    Args:
        points1: (N, 2) array of (lon, lat) in degrees, or a single (lon, lat).
        points2: (M, 2) array of (lon, lat) in degrees, or a single (lon, lat).

    Returns:
        (N, M) array; entry (i, j) is the distance from `points1[i]` to
        `points2[j]`, in km.
    """
    p1 = np.atleast_2d(np.asarray(points1, dtype=float))
    p2 = np.atleast_2d(np.asarray(points2, dtype=float))
    d2r = np.pi / 180.0
    lon1, lat1 = p1[:, 0, None], p1[:, 1, None]
    lon2, lat2 = p2[None, :, 0], p2[None, :, 1]
    # abs() on the longitude difference is a no-op (sin is odd, and it is squared);
    # kept because the MATLAB has it and removing it changes nothing.
    half_angles = np.arcsin(np.sqrt(
        np.sin((lat1 - lat2) * d2r / 2) ** 2
        + np.cos(lat1 * d2r) * np.cos(lat2 * d2r)
        * np.sin(np.abs(lon1 - lon2) * d2r / 2) ** 2))
    return 2 * EARTH_RADIUS_KM * np.sin(half_angles)
