"""
geod_to_geoc.py

Retrieves the geocentric latitude phi and the geodetic-direction altitude
h_phi (above a) for a given geodetic latitude phi_geod.

Note: h_geod returned here is != R_e (the Earth radius); see Vermeille's
algorithm reference below.

Examples (from the original MATLAB docstring)
-----------------------------------------------
    geod_to_geoc(60.4350632585099, Earth.r1_km, 0.8 * Earth.r1_km)
        -> (phi_geoc, h_geod) == (50, 831.969398327804)

    geod_to_geoc([0, 90], a, b)
        -> (phi_geoc, h_geod) == ([0, 90], [0, a - b])

Translated from the original MATLAB implementation. Only "Derivation #2"
is implemented, since "Derivation #1" was inside an `if 0` (dead-code /
reference-only) branch in the original source; it is kept below purely as
documentation.

Derivation #1 (reference only, not executed, valid only for h0 = 0)
---------------------------------------------------------------------
Using the identity  1 - e^2 sin(phi)^2 == a^2 / R_N^2  to simplify:

    x = (R_N + h_geod) * cos(phi_geod)                  # == x_sat
    z = (R_N * (1 - e2) + h_geod) * sin(phi_geod)        # == z_sat
    r = hypot(x, z)                                      # == a
    # [r, norm(R_sat), a] should all match

    h = h_geod  (sought solution), solving:
    a^2 == (R_N + h)^2 * cos(phi_geod)^2
         + (R_N * (1 - e2) + h)^2 * sin(phi_geod)^2
    for h via the quadratic form (B, C below), then:
    [R_sat, phi_geoc] = geod_to_pos(phi_geod, 0, h_geod, a, b)
"""

import numpy as np
from SOL_Tools.AstroConstants import *  # all variales
from SOL_Tools.Math_tools import * 

from geod_to_pos import geod_to_pos

def geod_to_geoc(phi_geod, a=6378.137, b=6356.75231424518, h0=0.0):
    """
    Parameters
    ----------
    phi_geod : array_like
        Geodetic latitude(s) [deg].
    a : float, optional
        Semi-major axis. Defaults to WGS84 Earth value (km).
    b : float, optional
        Semi-minor axis. Defaults to WGS84 Earth value (km).
    h0 : float, optional
        Altitude above `a` (i.e. rho = a + h0 is the radius of the
        reference sphere onto which the geodetic-normal line is
        projected). Default = 0.

    Returns
    -------
    phi_geoc : ndarray
        Geocentric latitude(s) [deg].
    h_geod : ndarray
        Geodetic-direction altitude(s) h_phi (Vermeille's algorithm).
        Note h_geod != R_e.
    """

    phi_geod = np.atleast_1d(np.asarray(phi_geod, dtype=float))
    rho = a + h0

    # --- Derivation #2 (valid for any h0) ----------------------------------

    # a) Sub-satellite point on the ellipsoid surface (h = 0)
    R_sub, _, _, R_N, *_ = geod_to_pos(phi_geod, 0.0, 0.0, a, b)
    x_sub = R_sub[0, :]
    z_sub = R_sub[2, :]

    # b) Normal line perpendicular to the ellipsoid at the sub-satellite point
    #    (direct solution, equivalent to the tangential-slope derivation
    #    commented out in the original MATLAB code)
    m = cotd(90 - phi_geod)
    b_ = z_sub - m * x_sub  # y-intercept (unrelated to the ellipsoid parameter b)

    # c) Intersect the normal line at R_sub with the sphere of radius rho
    #    m*xx + b_ == sqrt(rho^2 - xx^2)  @ x = x_sat
    x_sat = (-b_ * m + np.sqrt(rho**2 * (m**2 + 1) - b_**2)) / (m**2 + 1)
    z_sat = m * x_sat + b_

    h_geod = x_sat / cosd(phi_geod) - R_N          # h_phi (Vermeille's algorithm)
    phi_geoc = atan2d(z_sat, x_sat)                 # phi

    # Handle special/singular cases (poles and equator)
    k = np.isin(phi_geod, [-90.0, 0.0, 90.0])
    phi_geoc = np.where(k, phi_geod, phi_geoc)
    h_geod = np.where(phi_geod == 0.0, h0, h_geod)
    h_geod = np.where(np.abs(phi_geod) == 90.0, a + h0 - b, h_geod)

    # NOTE (from original source): could also be solved in (R, phi):
    #   y = R*sin(phi) + z0 ; y = sqrt((a+h0)^2 - R^2*cos(phi)^2)
    #   R^2*c1 + c2 == sqrt(c3 - R^2*c4)
    #   (R^2*c1 + c2)^2 == c3 - R^2*c4   -> quadratic solve for R
    # ==> verify +/- solution depending on quadrant, as in the original.

    return phi_geoc, h_geod


if __name__ == "__main__":
    # Test case from the original MATLAB docstring
    phi_geod = 50.0
    h0 = 1000.0

    a = 6378.137
    b = 6356.7523142

    phi_geoc, h_geod = geod_to_geoc(phi_geod, a, b, h0)
    print(f"phi_geoc = {phi_geoc[0]:.10f}  (expected ~49.8364099440923)")
    print(f"h_geod   = {h_geod[0]:.8f}     (expected ~1012.51026593283)")

    # Second documented example
    phi_geoc2, h_geod2 = geod_to_geoc(60.4350632585099, a, 0.8 * a)
    print(f"\nphi_geoc = {phi_geoc2[0]:.6f}  (expected 50)")
    print(f"h_geod   = {h_geod2[0]:.6f}  (expected 831.969398327804)")

    # Poles / equator edge cases
    phi_geoc3, h_geod3 = geod_to_geoc([0.0, 90.0], a, b)
    print(f"\nphi_geoc = {phi_geoc3}  (expected [0, 90])")
    print(f"h_geod   = {h_geod3}  (expected [0, a-b] = [0, {a - b}])")