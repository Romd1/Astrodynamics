import numpy as np
from SOL_Tools.Math_tools import * 
from SOL_Tools.Orbit_tools import * 

from geod_to_pos import geod_to_pos
from geoc_to_geod import geoc_to_geod

class Observer:
    """
    Represents the observer of the target.

    Inputs
    ------
    lat_geod : float
        Geodetic latitude [deg]
    lon : float
        Longitude [deg].
    h : float
        Geodetic altitude [km].

    Attributes
    --------------------
    Lat_geod : ndarray
        Geodetic latitude [deg].
    Lon : ndarray
        Longitude [deg].
    h : Geodetic altitude [km].
    pos : ndarray
        Observer ECEF position vector [x, y, z].
    Lat_geoc : ndarray
        Geocentric latitude [deg].
    R_e : ndarray
        Geocentric distance to the ellipsoid surface [km].
    z0 : ndarray
        Auxiliary z-coordinate [km].
    normal_vector : ndarray
        Unit normal vector to the ellipsoid surface.
    """

    def __init__(self, lat_geod, lon, h):
        self.Lat_geod = lat_geod  # [deg]
        self.Lon = lon            # [deg]
        self.h = h                # [km]

        R_geoc = geod_to_pos(self.Lat_geod, self.Lon, self.h)

        self.pos = np.asarray(R_geoc[0], dtype=float).reshape(3) # ECEF position [km]
        self.Lat_geoc = R_geoc[1]            # Geocentric latitude [deg]
        self.R_e = R_geoc[2]                 # Geocentric radius [km]
        self.z0 = R_geoc[6]                  # Auxiliary z-coordinate [km]

        self.normal_vector = RotZ(self.Lon) @ RotY(-self.Lat_geod) @ np.array([1, 0, 0]) # Unit normal vector, RHR for rot


    def range_vector(self, target):
        """Return ECEF range vector from observer to target [km]."""
        return target.pos - self.pos

    def aer(self, target):
        """Return Azimuth [deg], Elevation [deg], Range [km], SEZ [km]."""
        RV_ECEF = self.range_vector(target)

        return ECEF_to_AER(RV_ECEF, self.Lat_geod, self.Lon)



class Target:
    """
    Represents a target from its ECEF position.

    Inputs
    ------
    R_ECEF : array_like
        ECEF target position [x, y, z] [km].

    Attributes
    ----------
    pos : ndarray
        ECEF position vector [km].
    Lat_geoc : ndarray
        Geocentric latitude [deg].
    Lat_geod : ndarray
        Geodetic latitude [deg].
    h_geod : ndarray
        Geodetic altitude [km].
    Lon : ndarray
        Longitude [deg].
    Rho : ndarray
        Geocentric radius. [km]
    """

    def __init__(self, R_ECEF):
        self.pos = np.asarray(R_ECEF, dtype=float).reshape(3)

        # ECEF -> Geocentric and geodetic
        self.Lat_geod, self.h_geod, self.Lon, self.Lat_geoc, self.Rho = geoc_to_geod(R_ECEF)

    def subsatellite_point(self):
        """
        Computes the subsatellite point.

        Returns
        -------
        Lat_geod : ndarray
            Geodetic latitude [deg].
        h_geod : ndarray
            Geodetic altitude [km].
        R_subsat : ndarray
            ECEF coordinates [km]
        """
        Pos=geod_to_pos(self.Lat_geod, self.Lon)
        R_subsat=Pos[0]

        return self.Lat_geod, self.Lon, R_subsat