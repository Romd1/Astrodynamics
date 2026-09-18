"""
geod_to_pos.py

Converts geodetic latitude phi [deg] to a position vector [x; y; z] and to
geocentric latitude, at a given altitude h_phi above the ellipsoid.

Reference: Vallado, 2013, Eq. (3-7), r_delta & r_K;
from Hedgley, 1976, "An Exact Transformation from Geocentric to Geodetic
Coordinates for Nonzero Altitudes".

Note: geod_to_pos(lambda, phi_geoc, h_geoc, a, b) == Sph2Cart(lambda,
phi_geoc, a + h_geoc) when b == a.

Translated from the original MATLAB implementation.
"""

import numpy as np

from SOL_Tools.AstroConstants import *  # all variales
from SOL_Tools.Math_tools import * 


def geod_to_pos(phi_geod, lambda_=0.0, h_geod=0.0, a=Earth.r1_km, b=Earth.r2_km):
    """
    Parameters
    ----------
    phi_geod : array_like
        Geodetic latitude(s) [deg].
    lambda_ : array_like, optional
        Longitude(s) [deg]. Default = 0. Must be <= 360.
    h_geod : array_like, optional
        Geodetic altitude(s) h_phi above the ellipsoid. Default = 0.
    a : float, optional
        Semi-major axis. Defaults to WGS84 Earth value (km).
    b : float, optional
        Semi-minor axis. Defaults to WGS84 Earth value (km).

    Returns
    -------
    R_sat : ndarray, shape (3, N)
        ECEF position vector(s) [x; y; z].
    phi_geoc : ndarray
        Geocentric latitude(s) [deg].
    R_e : ndarray
        Distance from the ellipsoid center to its surface, measured along
        the geocentric radial direction, R_e(phi).
    R_N : ndarray
        Radius of curvature in the prime vertical (terminated by the
        minor axis), R_N(phi) -- independent of h_geod.
        R_N(phi=0) == a ; R_N(phi=90) == a^2/b ~= 2a - b.
    R_M : ndarray
        Radius of curvature in the meridian, R_M(phi).
        R_M(phi=0) == a*(1-e2) ; R_M(phi=90) == a^2/b.
    x0 : ndarray
        x-coordinate auxiliary quantity (see derivation in the original
        source): r_delta - z_sub * r_delta / r_K_.
    z0 : ndarray
        z-coordinate auxiliary quantity: z_sub - r_K_.

    Notes
    -----
    Two disabled ("if 0") blocks in the original MATLAB source -- one
    plotting the geodetic/geocentric latitude difference vs. altitude,
    one computing various mean-Earth-radius definitions -- are not
    executable code paths and are not translated here; see the comments
    near the bottom of this docstring for what they contained.

    Disabled plotting block (reference only):
        For phi_geod in [0, 90] deg and h_geod in {0, 500, ..., 2000} km,
        plots phi_geod - phi_geoc vs. phi_geod, to show how much
        difference results from using an oblate spheroid.

    Disabled mean-radius block (reference only):
        R1 = (2*a + b) / 3            # WGS84 (1984) mean radius
        R2 = (a**2 * b) ** (1/3)      # Wertz 2001 Sec. 9.1.5 (same-volume sphere)
        R_Gaussian = sqrt(R_N * R_M)
        R_Eulerian = R_N * R_M / (R_N * cosd(az)**2 + R_M * sind(az)**2)
    """

    phi_geod = np.atleast_1d(np.asarray(phi_geod, dtype=float))
    lambda_ = np.atleast_1d(np.asarray(lambda_, dtype=float))
    h_geod = np.atleast_1d(np.asarray(h_geod, dtype=float))

    e2 = 1 - (b / a) ** 2

    if np.any(lambda_ > 360):
        raise ValueError(">>> geod_to_pos: second argument (lambda) out of bound?!")

    # R_N(phi): radius of curvature in the prime vertical (terminated by
    # the minor axis) -- independent of h_geod
    R_N = a / np.sqrt(1 + ((b / a) ** 2 - 1) * sind(phi_geod) ** 2)
    # == a / sqrt(1 - e2 * sind(phi_geod)**2)
    # == a**2 / sqrt(a**2 * cosd(phi_geod)**2 + b**2 * sind(phi_geod)**2)
    # R_N(phi=0)  == a
    # R_N(phi=90) == a**2/b ~= 2a - b

    # R_N(phi0):
    # phi_geoc0 = atand(tand(phi_geod) * (b/a)**2)
    # num = a * sqrt(a**4 * sind(phi_geoc0)**2 + b**4 * cosd(phi_geoc0)**2)
    # den = b * sqrt(a**2 * sind(phi_geoc0)**2 + b**2 * cosd(phi_geoc0)**2)
    # R_N = num / den
    # R_N = (a * sqrt(1 - (2*e2 - e2**2) * cosd(phi_geoc0)**2)) / \
    #       (sqrt(1 - e2) * sqrt(1 - e2 * cosd(phi_geoc0)**2))

    # --- (x, y, z) = f(phi) -------------------------------------------------
    x = (R_N + h_geod) * cosd(phi_geod) * cosd(lambda_)   # r_delta @ h_phi
    y = (R_N + h_geod) * cosd(phi_geod) * sind(lambda_)
    z = ((b / a) ** 2 * R_N + h_geod) * sind(phi_geod)    # r_K @ h_phi
    # -------------------------------------------------------------------------

    R_sat = np.vstack([x, y, z])
    # Rho = norm(R_sat)

    # [~, phi_geoc] = Cart2Sph(R_sat)
    # phi_geoc = atan2d(z, hypot(x, y)) == asind(z / norm(R_sat))

    # --- phi, h_phi = f(R(phi)) ----------------------------------------------
    with np.errstate(divide="ignore", invalid="ignore"):
        phi_geoc = atand(tand(phi_geod) * (R_N * (b / a) ** 2 + h_geod) / (R_N + h_geod))
    # -------------------------------------------------------------------------
    # == atand(tand(phi_geod) * (b/a)**2) only when h_phi == 0 -> phi0
    # h_geod == hypot(x, y) / cosd(phi_geod) - R_N

    z_sub = (b / a) ** 2 * R_N * sind(phi_geod)

    r_delta = R_N * cosd(phi_geod)
    r_K_ = R_N * sind(phi_geod)  # extending down to z0
    z0 = z_sub - r_K_
    with np.errstate(divide="ignore", invalid="ignore"):
        x0 = r_delta - z_sub * r_delta / r_K_

    r_K = R_N * (1 - e2) * sind(phi_geod)  # noqa: F841 (kept for parity with source)
    # r_sub = sqrt(r_delta**2 + r_K**2) == R_e

    # R_e(phi) = distance from the ellipsoid center to its surface, measured
    # along the geocentric radial direction
    R_e = a * np.sqrt(
        ((1 - e2) ** 2 * sind(phi_geod) ** 2 + cosd(phi_geod) ** 2)
        / (1 - e2 * sind(phi_geod) ** 2)
    )
    # == a * sqrt((1 - (2*e2 - e2**2) * sind(phi_geod)**2) / (1 - e2 * sind(phi_geod)**2))
    # == R_N * sqrt(cosd(phi_geod)**2 + (b/a)**4 * sind(phi_geod)**2)
    # == sqrt((a**4*cosd(phi_geod)**2 + b**4*sind(phi_geod)**2) / (a**2*cosd(phi_geod)**2 + b**2*sind(phi_geod)**2))
    # == norm(R_sub)

    # R_e(phi_geoc):
    # phi_geoc0 = atand(tand(phi_geod) * (b/a)**2)
    # R_e = a / sqrt(1 + ((a/b)**2 - 1) * sind(phi_geoc0)**2)
    #    == (a * sqrt(1 - e2)) / sqrt(1 - e2 * cosd(phi_geoc0)**2)
    # R_e = (a * b) / sqrt(a**2 * sind(phi_geoc0)**2 + b**2 * cosd(phi_geoc0)**2)

    # --- Radius of curvature in the meridian (Meridian Radius of Curvature) ---
    R_M = a * (1 - e2) / (1 - e2 * sind(phi_geod) ** 2) ** 1.5  # R_M(phi)
    # == R_N * (1 - e2) / (1 - e2 * sind(phi_geod)**2)
    # R_M(phi=0)  == a * (1 - e2)
    # R_M(phi=90) == a**2/b

    return R_sat, phi_geoc, R_e, R_N, R_M, x0, z0


if __name__ == "__main__":
    R_sat, phi_geoc, R_e, R_N, R_M, x0, z0 = geod_to_pos(50.0, 0.0, 1000.0)
    print("R_sat    =", R_sat.ravel())
    print("phi_geoc =", phi_geoc[0])
    print("R_e      =", R_e[0])
    print("R_N      =", R_N[0])
    print("R_M      =", R_M[0])
    print("x0       =", x0[0])
    print("z0       =", z0[0])