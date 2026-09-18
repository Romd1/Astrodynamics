"""
geoc_to_geod.py

(ECEF) Geodetic coordinate conversion: (x, y, z) -> (phi, h)
h_geod = deviation with respect to a perfect sphere
         (0 km at the equator, ~21.38 km at the poles for Earth/WGS84)

Note: the largest difference between geocentric (phi) and geodetic (phi_geod)
latitude is at 45 deg.

References
----------
Shu, 2010, "An iterative algorithm to compute geodetic coordinates"
    (Table 1: comparison of arithmetic operations - latitude and height)
Escobal, Pedro Ramon: "Methods of Orbit Determination",
    John Wiley & Sons, Inc., 1965, p. 23.
Heikkinen, 1982 (closed-form, non-iterative)
Zhu, 1993, "Direct, exact transformation from ECEF to geodetic coordinates"
Borkowski, 1987, "Transformation of geocentric to geodetic coordinates
    without approximations"
Bowring, iterative method (see CelesTrak, https://celestrak.org/columns/v02n03/)
Clynch, 2006, "Geodetic Coordinate Conversions"

Translated from the original MATLAB implementation.
"""

import numpy as np
from SOL_Tools.AstroConstants import *  # all variales
from SOL_Tools.Math_tools import * 


def geoc_to_geod(R_ECEF, a=Earth.r1_km, b=Earth.r2_km, method=2, tol=1e-14):
    """
    Parameters
    ----------
    R_ECEF : (3, N) array_like
        ECEF Cartesian coordinates [x; y; z] (any consistent length unit).
    a : float, optional
        Semi-major axis. Defaults to WGS84 Earth value (km).
    b : float, optional
        Semi-minor axis. Defaults to WGS84 Earth value (km).
    method : int, optional
        Conversion method to use (see below). Default = 2 (Heikkinen,
        exact non-iterative, recommended).
    tol : float, optional
        Convergence tolerance for iterative methods (must be < 1e-15... in
        practice values around 1e-14 work well).

    Returns
    -------
    phi_geod : ndarray
        Geodetic latitude [deg].
    h_geod : ndarray
        Geodetic height / altitude above the reference ellipsoid.
    Lon : ndarray
        Geocentric longitude [deg] (== geodetic longitude).
    Lat : ndarray
        Geocentric latitude [deg].
    Rho : ndarray
        Geocentric radius.

    Methods
    -------
    2  : Heikkinen's Algorithm (1982)      -- exact, non-iterative (default, BEST)
    3  : Zhu (1993)                        -- exact, non-iterative
    4  : Borkowski (1987)                  -- exact, non-iterative
    10 : Bowring, classic iterative scheme
    11 : Bowring, fixed-point iteration (fast convergence, ~2-3 iters)
    12 : Clynch (2006), iterative
    13 : Fukushima-style Halley iteration

    Notes
    -----
    Methods 0 (MATLAB Aerospace Toolbox) and 1 (Long 1974 / NASA TN D-7522
    fast approximation) from the original MATLAB code depend on external
    functions (`geoc2geod`, `Geodetic_approx`) not translated here, and are
    not implemented in this Python port. Likewise, methods 14 ("Nievergelt
    & Keeler"), 20 (series expansion), 990 and 991 were incomplete /
    "study" branches in the original MATLAB file (marked TODO) and are
    not implemented.
    """

    R_ECEF = np.asarray(R_ECEF, dtype=float)
    if R_ECEF.ndim == 1:
        R_ECEF = R_ECEF.reshape(3, 1)


    e2 = 1 - (b / a) ** 2     # first eccentricity squared
    e_2 = (a / b) ** 2 - 1    # second eccentricity squared == e2 / (1 - e2)

    # --- 0) Geocentric longitude, latitude, and radius --------------------
    Lon, Lat, Rho = Cart2Sph(R_ECEF)

    # --- 1) Geodetic latitude and altitude ---------------------------------
    h_geod = None

    if method in (0, 1) and a != 6378.137:
        raise ValueError(
            "geoc_to_geod: Methods 0 & 1 are for Earth only (not applicable "
            "for normalized radii) => use method 2 instead."
        )

    r = np.hypot(R_ECEF[0, :], R_ECEF[1, :])   # distance from Z axis
    z = R_ECEF[2, :]
    m = r.size

    if method == 0:
        raise NotImplementedError(
            "Method 0 relies on MATLAB's Aerospace Toolbox (geoc2geod); "
            "not available in this Python port. Use method 2 instead."
        )

    elif method == 1:
        raise NotImplementedError(
            "Method 1 relies on 'Geodetic_approx' (Long 1974, NASA "
            "TN D-7522), not translated here. Use method 2 instead."
        )

    elif method == 2:
        # Heikkinen's Algorithm (1982) - closed form, non-iterative
        # Sub-millimeter accuracy, stable at poles/center.
        F = 54 * b**2 * z**2
        G = r**2 + (1 - e2) * z**2 - e2 * (a**2 - b**2)
        c = (e2**2 * F * r**2) / (G**3)

        s = (1 + c + np.sqrt(c**2 + 2 * c)) ** (1 / 3)
        k = s + 1 + 1 / s
        P = F / (3 * k**2 * G**2)
        Q = np.sqrt(1 + 2 * e2**2 * P)

        term1 = -(P * e2 * r) / (1 + Q)
        term2 = 0.5 * a**2 * (1 + 1 / Q)
        term3 = (P * (1 - e2) * z**2) / (Q * (1 + Q))
        term4 = 0.5 * P * r**2

        r0 = term1 + np.sqrt(np.maximum(0.0, term2 - term3 - term4))

        U = np.sqrt((r - e2 * r0) ** 2 + z**2)
        V = np.sqrt((r - e2 * r0) ** 2 + (1 - e2) * z**2)

        z0 = (b**2 * z) / (a * V)

        h_geod = U * (1 - b**2 / (a * V))
        phi_geod = atan2d(z + e_2 * z0, r)

    elif method == 3:
        # Zhu (1993) - exact, closed form
        l = e2 / 2
        m_ = (r / a) ** 2
        n = ((1 - e2) * z / b) ** 2
        i = -(2 * l**2 + m_ + n) / 2
        k = l**2 * (l**2 - m_ - n)
        q = (m_ + n - 4 * l**2) ** 3 / 216 + m_ * n * l**2
        D = np.sqrt((2 * q - m_ * n * l**2) * m_ * n * l**2)
        beta = i / 3 - np.cbrt(q + D) - np.cbrt(q - D)
        t = (np.sqrt(np.sqrt(beta**2 - k) - (beta + i) / 2)
             - np.sign(m_ - n) * np.sqrt((beta - i) / 2))
        w1 = r / (t + l)
        z1 = (1 - e2) * z / (t - l)

        phi_geod = atan2d(z1, (1 - e2) * w1)
        h_geod = np.sign(t - 1 + l) * np.sqrt((r - w1) ** 2 + (z - z1) ** 2)

    elif method == 4:
        # Borkowski (1987) - exact, non-iterative (element-wise, uses a
        # cubic-equation solve per point, as in the original MATLAB loop)
        h_geod = np.zeros(m)
        phi_geod = np.zeros(m)

        for idx in range(m):
            zi = z[idx]
            ri = r[idx]

            if ri < 1e-12:
                # Special polar case
                phi_geod[idx] = np.sign(zi) * 90.0
                h_geod[idx] = abs(zi) - b
                continue

            E = ((zi + b) * b / a - a) / ri
            F = ((zi - b) * b / a + a) / ri

            P = 4 / 3 * (E * F + 1)
            Q = 2 * (E**2 - F**2)
            Dd = P**3 + Q**2

            if Dd >= 0:
                v = np.cbrt(np.sqrt(Dd) - Q) - np.cbrt(np.sqrt(Dd) + Q)
            else:
                v = 2 * np.sqrt(-P) * np.cos(1 / 3 * np.arccos(Q / (P * np.sqrt(-P))))

            if abs(v**2) < abs(P):
                v = -(v**3 + 2 * Q) / (3 * P)

            G = 0.5 * (np.sqrt(E**2 + v) + E)
            t = np.sqrt(G**2 + (F - v * G) / (2 * G - E)) - G

            phi = atan2d(a * (1 - t**2), 2 * b * t)
            phi_geod[idx] = phi

            sinLat = sind(phi)
            R_N = a / np.sqrt(1 - e2 * sinLat**2)
            h_geod[idx] = ri * cosd(phi) + (zi + e2 * R_N * sinLat) * sinLat - R_N

    elif method == 10:
        # Bowring, classic iterative scheme
        h_geod = np.zeros(m)
        phi_geod = np.zeros(m)

        for idx in range(m):
            phi = Lat[idx]
            ri = r[idx]
            zi = z[idx]

            phi_geod_ = np.inf
            R_N = a  # will be reassigned in the loop
            while abs(phi_geod_ - phi) > tol:
                phi_geod_ = phi
                R_N = a / np.sqrt(1 - e2 * sind(phi_geod_) ** 2)
                phi = atan2d(zi + e2 * R_N * sind(phi_geod_), ri)

            phi_geod[idx] = phi
            h_geod[idx] = ri / cosd(phi) - R_N

    elif method == 11:
        # Fixed-point iteration of Bowring's formula (fast, ~2-3 iterations)
        phi_geod = np.zeros(m)

        for idx in range(m):
            phi = Lat[idx]
            ri = r[idx]
            zi = z[idx]

            phi_geod_ = phi
            while True:
                beta = atand(tand(phi_geod_) * b / a)
                phi_geod_new = atand((zi + b * e_2 * sind(beta) ** 3) /
                                      (ri - a * e2 * cosd(beta) ** 3))
                if abs(phi_geod_new - phi) < tol:
                    phi_geod_ = phi_geod_new
                    break
                phi = phi_geod_new
                phi_geod_ = phi_geod_new

            phi_geod[idx] = phi_geod_
        h_geod = None  # computed generically below

    elif method == 12:
        # Clynch (2006), iterative
        phi_geod = np.zeros(m)
        h_geod = np.zeros(m)

        for idx in range(m):
            ri = r[idx]
            zi = z[idx]
            phi = Lat[idx]
            phi0 = phi
            h = 0.0

            for _ in range(10):  # typically converges in <= 5 iterations
                R_N = a / np.sqrt(1 - e2 * sind(phi) ** 2)
                h = ri / cosd(phi) - R_N
                phi = atand(zi / (ri * (1 - e2 * R_N / (R_N + h))))
                if abs(phi - phi0) < tol:
                    break
                phi0 = phi

            phi_geod[idx] = phi
            h_geod[idx] = h

    elif method == 13:
        # Fukushima-style Halley iteration
        phi_geod = np.zeros(m)
        h_geod = np.zeros(m)

        for idx in range(m):
            ri = r[idx]
            zi = z[idx]

            if ri < 1e-15:
                phi_geod[idx] = np.sign(zi) * 90.0
                h_geod[idx] = abs(zi) - b
                continue

            lat = np.arctan2(zi, ri * (1 - e2))  # initial guess (radians)

            for _ in range(6):  # 6 iters -> ~1e-14 accuracy
                R_N = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
                h = ri / np.cos(lat) - R_N
                lat_new = np.arctan2(zi, ri * (1 - e2 * R_N / (R_N + h)))
                if abs(lat_new - lat) < tol:
                    lat = lat_new
                    break
                lat = lat_new

            R_N = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
            h_geod[idx] = ri / np.cos(lat) - R_N
            phi_geod[idx] = DEG(lat)

    else:
        raise ValueError(f"geoc_to_geod: Unknown or unimplemented method #{method}")

    phi_geod = np.atleast_1d(np.asarray(phi_geod, dtype=float))

    R_N = a / np.sqrt(1 - e2 * sind(phi_geod) ** 2)  # R_N(phi_geod)

    if h_geod is None:
        h_geod = np.hypot(R_ECEF[0, :], R_ECEF[1, :]) / cosd(phi_geod) - R_N
    else:
        h_geod = np.atleast_1d(np.asarray(h_geod, dtype=float))

    return phi_geod, h_geod, Lon, Lat, Rho


# ----------------------------------------------------------------------
# Simple self-test / example (mirrors the MATLAB docstring test case)
# ----------------------------------------------------------------------

if __name__ == "__main__":

    phi_geoc = 50.0
    rho = Earth.r1_km + 1000.0
    R_ECEF = Sph2Cart(0.0, phi_geoc, rho)

    phi_geod, h_geod, Lon, Lat, Rho = geoc_to_geod(R_ECEF)

    print(f"phi_geod = {phi_geod[0]:.10f} deg   (expected ~50.1634243856166)")
    print(f"h_geod   = {h_geod[0]:.8f} km       (expected ~1012.57038273569)")
    print(f"Rho      = {Rho[0]:.8f}             (expected == rho = {rho})")

    # Quick cross-check across all implemented exact/iterative methods
    print("\nCross-check across methods for the same point:")
    for method in (2, 3, 4, 10, 11, 12, 13):
        phi, h, *_ = geoc_to_geod(R_ECEF, method=method)
        print(f"  method {method:2d}: phi_geod = {phi[0]:.10f}  h_geod = {h[0]:.8f}")
