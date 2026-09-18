"""
geoc_to_pos.py

Converts geocentric latitude phi & longitude lambda [deg] to position
vector [x; y; z], optionally at a given altitude h above the ellipsoid
in the geocentric direction (h_phi).

Reference: Vallado, 2013, Eq. (3-10), r_delta & r_K.

Although not used as often as geod_to_pos, it has certain advantages in
some situations.
-> Result is approximate when h_geod (not measured directly along the
   geocentric radius) is provided instead of h_geoc (but still very good,
   since h_phi ~= h_phi_geodetic).

Note: geoc_to_pos(lambda, phi_geoc, h_geoc, a, b) == Sph2Cart(lambda,
phi_geoc, a + h_geoc) when b == a.

Translated from the original MATLAB implementation.
"""

import numpy as np

from SOL_Tools.AstroConstants import *  # all variales
from SOL_Tools.Math_tools import * 


def geoc_to_pos(phi_geoc, lambda_=0.0, h_geoc=0.0, a=Earth.r1_km, b=Earth.r2_km):
    """
    Parameters
    ----------
    phi_geoc : array_like
        Geocentric latitude(s) [deg].
    lambda_ : array_like, optional
        Geocentric longitude(s) [deg]. Default = 0.
    h_geoc : array_like, optional
        Altitude(s) above the ellipsoid, measured along the geocentric
        radius direction. Default = 0.
    a : float, optional
        Semi-major axis. Defaults to WGS84 Earth value (km).
    b : float, optional
        Semi-minor axis. Defaults to WGS84 Earth value (km).

    Returns
    -------
    R : ndarray, shape (3, N)
        ECEF position vector(s) [x; y; z].
    r_geoc : ndarray
        Geocentric radius of the ellipsoid surface (h = 0) at phi_geoc,
        i.e. distance from the ellipsoid center to its surface
        (r_geoc + h_geoc == Rho, the geocentric radius to the point).
    """

    phi_geoc = np.atleast_1d(np.asarray(phi_geoc, dtype=float))
    lambda_ = np.atleast_1d(np.asarray(lambda_, dtype=float))
    h_geoc = np.atleast_1d(np.asarray(h_geoc, dtype=float))

    # Radius of the ellipsoid R(phi) from its center to its surface (h = 0)
    # in function of phi; Vallado 2013, Eq. (3-10)
    r_geoc = a / np.sqrt(1 + ((a / b) ** 2 - 1) * sind(phi_geoc) ** 2)  # != R_N
    
    # = a / sqrt(1 + e_2 * sind(phi_geoc)**2)          # e_2 = (a/b)^2 - 1  (second eccentricity squared)
    # = a * sqrt((1 - e2) / (1 - e2 * cosd(phi_geoc)**2))   # <- Vallado (3-10)
    # r_geoc = GeocentricRadius(phi_geoc, a, b)
    #        = a / sqrt(1 + ((a/b)^2 - 1) * sind(phi_geoc)**2)
    #        = a / sqrt(1 + e2 / (1 - e2) * sind(phi_geoc)**2)
    #        = a / sqrt(1 + e_2 * sind(phi_geoc)**2)
    #        = a * sqrt(1 / (1 + (1/(1 - f)**2 - 1) * sind(phi_geoc)**2))
    # r_geoc + h_geoc == Rho

    # x = (r_geoc + h_geoc) * cosd(phi_geoc) * cosd(lambda_)
    # y = (r_geoc + h_geoc) * cosd(phi_geoc) * sind(lambda_)
    # z = (r_geoc + h_geoc) * sind(phi_geoc)

    R = Sph2Cart(lambda_, phi_geoc, r_geoc + h_geoc)
    return R, r_geoc


if __name__ == "__main__":
    # Quick sanity check
    R, r_geoc = geoc_to_pos(50.0, 0.0, 1000.0)
    print("R =", R.ravel())
    print("r_geoc =", r_geoc[0])
    print("Rho (r_geoc + h_geoc) =", r_geoc[0] + 1000.0)