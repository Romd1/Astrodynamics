# ┌───────────────────────────┐
# │  SOL_Tools\Math_tools.py  │ 
# └───────────────────────────┘
#  © SpaceOrbitLAB

import numpy as np
from geopy.distance import great_circle
#import math

def DEG(a):
# Transforms an angle in radians into degrees
    return np.rad2deg(a)
def RAD(a):
# Transforms an angle in degrees into radians
    return np.deg2rad(a)

def sind(angle):
    return np.sin(np.deg2rad(angle))
def cosd(angle):
    return np.cos(np.deg2rad(angle))
def tand(angle):
    return np.tan(np.deg2rad(angle))
def cotd(x):
    return 1.0 / np.tan(np.radians(x))

def asind(x):
    return np.rad2deg(np.arcsin(x))
def acosd(x):
    return np.rad2deg(np.arccos(x))
def atan2d(y, x):
    return np.rad2deg(np.arctan2(y, x))
def atand(x):
    return np.degrees(np.arctan(x))

def norm(v):
# Computes the norm of a vector
    return np.linalg.norm(v)

def RotX(angle):
# Rotation matrix around X axis (3×3)  X′ = R * X
# As standard mathematical convention, https://en.wikipedia.org/wiki/Rotation_matrix
#
# like Battin 1999, Chobotov 2002 (4.108), but unlike Bate 1971 (Eq.2.6-8) and Curtis 2020 (Eq. 4.32)
# who use a different convention (Vallado 2013, Section 3.4.1 p.162-168)

# The transformation matrix R corresponding to a single rotation of the coordinate frame about the positive X axis
# through a positive angle t_deg, counterclockwise, or Right-Hand-Rule (RHR) is:

    c = cosd(angle)
    s = sind(angle)

    return np.array([
        [1.0, 0.0, 0.0],
        [0.0,   c,  -s],
        [0.0,   s,   c]
    ])

def RotY(angle):
# Rotation matrix around Y axis (3×3), counterclockwise (RHR)
# https://en.wikipedia.org/wiki/Rotation_matrix
    
    c = cosd(angle)
    s = sind(angle)

    return np.array([
        [  c, 0.0,   s],
        [0.0, 1.0, 0.0],
        [ -s, 0.0,   c]
    ])

def RotZ(angle):
# Rotation matrix around Z axis (3×3), counterclockwise (RHR)
# https://en.wikipedia.org/wiki/Rotation_matrix

    c = cosd(angle)
    s = sind(angle)

    return np.array([
        [ c, -s, 0.0],
        [ s,  c, 0.0],
        [0.0, 0.0, 1.0]
    ])


def Sph2Cart(az, el, r = 1.0, stack = True):
    """Transform spherical to Cartesian coordinates, angles in degrees.

    Parameters
    ----------
    az : array_like
        Counterclockwise angle (alpha) in the xy plane, measured from the positive x axis. Degrees.
    el : array_like
        Elevation angle (delta) from the xy plane. Degrees.
    r : array_like, default 1.0
        Radius. Scalar, or broadcastable against `az` and `el`.
    stack : bool, default False
        If False, return the tuple ``(x, y, z)``.
        If True, return a single ``(3, n)`` array whose rows are x, y, z.
        This replaces the MATLAB ``nargout`` switch, which has no Python
        equivalent: the caller must say which form it wants.

    Returns
    -------
    (x, y, z) : tuple of ndarray, shape (n,)
        When ``stack`` is False.
    xyz : ndarray, shape (3, n)
        When ``stack`` is True.

    Examples
    --------
    >>> x, y, z = Sph2Cart(0.0, 90.0)
    >>> round(float(z[0]), 12)
    1.0
    >>> Sph2Cart([0, 90], [0, 0], 1, stack=True).shape
    (3, 2)
    """
    az = np.ravel(np.asarray(az, dtype=float))
    el = np.ravel(np.asarray(el, dtype=float))
    r = np.asarray(r, dtype=float)

    rcoselev = r * cosd(el)
    x = rcoselev * cosd(az)
    y = rcoselev * sind(az)
    z = r * sind(el)

    if stack:
        return np.vstack((x, y, z))
    else:
        return x, y, z


def Cart2Sph(xyz):
    """
    Convert Cartesian unit-vector coordinates to spherical longitude/latitude.

    Parameters
    ----------
    xyz : ndarray
        3 x N Cartesian coordinates.

    Returns
    -------
    lon : ndarray
        Longitude / right ascension / azimuth [deg], in [0, 360).

    lat : ndarray
        Latitude / declination / elevation [deg].

    r : ndarray
        Range
    """

    x = xyz[0, :]
    y = xyz[1, :]
    z = xyz[2, :]

    r = np.linalg.norm(xyz, axis=0)

    lon = np.degrees(np.arctan2(y, x)) % 360.0
    lat = np.degrees(np.arcsin(z / r))

    return lon, lat, r


def GreatCircle2(lat1, lon1, lat2, lon2, t=None):

    """
    if True:
        # Calculate great-circle distance
        distance_miles = great_circle(nyc, london).miles
        distance_km = great_circle(nyc, london).kilometers
    else:
        do_the_math

    Generate 3D great-circle coordinates on a unit sphere passing through two points.

    Parameters
    ----------
    lat1, lon1 : float
        Latitude and longitude of first point [deg].

    lat2, lon2 : float
        Latitude and longitude of second point [deg].

    t : ndarray, optional
        Angles along the great-circle arc [deg].
        If omitted, t goes from 0 to the angular separation between the two points.

    Returns
    -------
    GC : dict (OLD) → now class
        
    """

    # Angular separation Δσ
    cos_arclen = (sind(lat1) * sind(lat2) + cosd(lat1) * cosd(lat2) * cosd(lon1 - lon2))

    # Numerical protection against tiny roundoff outside [-1, 1]
    cos_arclen = np.clip(cos_arclen, -1.0, 1.0)

    GC_arclen = acosd(cos_arclen)

    # Initial azimuth angle of the great circle
    GC_az = atan2d(cosd(lat2) * sind(lon2 - lon1), cosd(lat1) * sind(lat2) - sind(lat1) * cosd(lat2) * cosd(lon2 - lon1))

    # Default sampling along the arc from point 1 to point 2 in degrees
    if t is None:
        N = 100
        t = np.linspace(0.0, GC_arclen, N)
    else:
        t = np.asarray(t, dtype=float)

    # Endpoint unit vectors
    print('lon1, lat1 = ', lon1, lat1)
    x1, y1, z1 = Sph2Cart(lon1, lat1, 1.0, False)
    print('xyz1 = ', x1, y1, z1)
    x2, y2, z2 = Sph2Cart(lon2, lat2, 1.0, False)

    denom = sind(GC_arclen)

    if np.isclose(denom, 0.0):
        raise ValueError(
            "The two points are identical or antipodal. "
            "The great circle is not uniquely defined."
        )

    # Parametric formulas for intermediate points along a great circle
    c1 = sind(GC_arclen - t) / denom
    c2 = sind(t) / denom

    GC_xyz = np.vstack((
        c1 * x1 + c2 * x2,
        c1 * y1 + c2 * y2,
        c1 * z1 + c2 * z2,
    ))

    GC_lon, GC_lat, _ = Cart2Sph(GC_xyz)

    # Great-circle navigation quantities
    lon12 = lon2 - lon1

    GC_alpha1 = atan2d(cosd(lat2) * sind(lon12), cosd(lat1) * sind(lat2) - sind(lat1) * cosd(lat2) * cosd(lon12))
    GC_alpha2 = atan2d(cosd(lat1) * sind(lon12), -cosd(lat2) * sind(lat1) + sind(lat2) * cosd(lat1) * cosd(lon12))

    # Longitude at the node alpha_0
    GC_alpha0 = atan2d(sind(GC_alpha1) * cosd(lat1), np.sqrt(cosd(GC_alpha1)**2 + (sind(GC_alpha1) * sind(lat1))**2))

    GC_lat01 = atan2d(tand(lat1), cosd(GC_alpha1))
    GC_lon01 = atan2d(sind(GC_alpha0) * sind(GC_lat01), cosd(GC_lat01))

    GC_lon0 = lon1 - GC_lon01
    GC_lat0 = 0.0

    # Stores data in output class
    class GC:
        arclen = GC_arclen
        az = GC_az
        xyz = GC_xyz
        lon = GC_lon
        lat = GC_lat
        alpha1 = GC_alpha1
        alpha2 = GC_alpha2
        alpha0 = GC_alpha0
        lat01 = GC_lat01
        lon01 = GC_lon01
        lon0 = GC_lon0
        lat0 = GC_lat0
      
    return GC
        
# end GreatCircle()


def SmallCircle(az, el, alpha, t=None):
    """
    Generates a unit 3D Small Circle on a sphere.

    Parameters
    ----------
    az : float
        RA / longitude / azimuth of the circle center [°].

    el : float
        Declination / latitude / elevation of the circle center [°].

    alpha : float
        Central angle radius of the small circle [°].

    t : array_like, optional
        Angles along the circle [°].
        Default is 100 points from 0 to 360°.
        Origin is at East, measured counterclockwise.

    Returns
    -------
    lat : ndarray
        Latitude / declination / elevation along the small circle [deg].

    lon : ndarray
        Longitude / RA / azimuth along the small circle [deg].

    xyz : ndarray
        3 x N Cartesian unit vectors along the small circle.
    """

    if t is None:
        t = np.linspace(0.0, 360.0, 100)

    t = np.asarray(t, dtype=float)

    # Expanded algebraic formula
    x = (-sind(alpha) * sind(el) * cosd(az) * sind(t) - sind(alpha) * sind(az) * cosd(t) + cosd(alpha) * cosd(el) * cosd(az))
    y = (-sind(alpha) * sind(el) * sind(az) * sind(t) + sind(alpha) * cosd(az) * cosd(t) + cosd(alpha) * cosd(el) * sind(az))
    z = (sind(alpha) * cosd(el) * sind(t) + cosd(alpha) * sind(el))

    xyz = np.vstack((x, y, z))

    lon, lat, _ = Cart2Sph(xyz)

    return lat, lon, xyz

# end of SmallCircle


def Angle(R1, R2):
# Computes the angle between 2 vectors in degrees
#   R1,2 = 3×n or 1×3  → works in 2D and 3D
# Returns the acute or smallest angle: 0 ≤ α ≤ 180°

    # Angle: input vectors shall be of size 3×n.');
    #end

    alpha = np.arccos(np.dot(R1, R2) / (norm(R1) * norm(R2)))  # normalization when using non unit vectors
    #alpha = real(alpha);  # for numerical instabilities for quasi-colinear vectors 

    return DEG(alpha)

# end of Angle()