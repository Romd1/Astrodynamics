# ┌────────────────────────────┐
# │  SOL_Tools\Orbit_tools.py  │ 
# └────────────────────────────┘
#  © SpaceOrbitLAB

import numpy as np
import math
from scipy.integrate import solve_ivp
from datetime import datetime, timezone, timedelta

from SOL_Tools.AstroConstants import Earth
from SOL_Tools.Math_tools import *


def OrbitPropagation_2BN(SV, time_s):
# Computes orbit solving Newton's equation on time vector

    # → see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html

      return ECI_pos_2BN, ECI_vel_2BN

def Rates_Newton(t, z):
    return
  

def OrbitPropagation_2BK(COE, t_array, mu):
# Two-Body Keplerian Orbit Propagation from COEs
    """
    Propagate an elliptic orbit using the two-body Keplerian method.

    Inputs:
        a       : semi-major axis [km]
        e       : eccentricity [-]
        i       : inclination [°]
        RAAN    : right ascension of ascending node [°]
        w       : argument of perigee [°]
        M0      : Mean Anomaly at epoch [°]
        t_array : time since epoch [s]
        mu      : gravitational parameter µ [km³/s²]

    Outputs:
        ECI_pos : position history, shape (N, 3) [km]
        ECI_vel : velocity history, shape (N, 3) [km/s]
        nu_all  : true anomaly history [rad]
    """

    if e >= 1.0:
        raise ValueError(">>> OrbitPropagation_2BK: This implementation is for elliptic orbits only: e < 1.")

    M0_rad = np.deg2rad(M0)

    # Mean motion
    n = np.sqrt(mu / a**3)  # [rad/s]
    # print('n = ', n)

    # Pre-allocate output arrays  (use float instead of double to save space)
    N = len(t_array)
    ECI_pos = np.zeros((3, N), dtype=float)
    ECI_vel = np.zeros((3, N), dtype=float)
    nu_all  = np.zeros(N, dtype=float)

    # Point-by-point calculation along the orbit
    # KeplerSolver is assumed to operate on scalar mean anomalies
    for k, t in enumerate(t_array):  # (do not use i as a loop counter as here it represents inclination)

        # Mean anomaly at time t
        M = M0_rad + n * t  # [rad]

        # Solve Kepler equation
        E, nu = KeplerSolver(M, e)  # [rad]
        nu = np.rad2deg(nu)  # [°]

        if False: # k == 0:
            print('M = ', np.rad2deg(M))
            print('E = ', np.rad2deg(E))
            print('ν = ', nu)

        # Convert COEs to inertial state
        R, V = COE_to_SV(a, e, i, RAAN, w, nu, mu)   

        # Store directly as columns in 3 x N arrays
        ECI_pos[:, k] = np.asarray(R, dtype=float).reshape(3)
        ECI_vel[:, k] = np.asarray(V, dtype=float).reshape(3)

        # Store true anomaly
        nu_all[k] = nu
        # (could also store Eccentric Anomaly and return as a vector)

    return ECI_pos, ECI_vel, nu_all


def KeplerSolver(M, e, tol=1e-12, max_iter=50):
    """
    Solve Kepler's equation:  M = E - e sin(E) 

    Inputs:
        M : mean anomaly [rad]
        e : eccentricity [n.u.]

    Output:
        E : eccentric anomaly [rad]
        ν : True Anomaly [rad]

    *** All angles are in radians here ***
    """

    # Wraps M to improve convergence
    M = np.mod(M, 2*np.pi)

    # Initial guess
    if e < 0.8:
        E = M
    else:
        E = np.pi

    #--- Newton’s method to solve Kepler Equation
    for _ in range(max_iter):
        f = E - e * np.sin(E) - M
        fp = 1.0 - e * np.cos(E)

        dE = -f / fp
        E = E + dE

        if abs(dE) < tol:
            break

    #--- Converts Eccentric Anomaly E to True Anomaly ν
    nu = 2 * np.atan(np.tan(E/2) * np.sqrt((1 + e) / (1 - e)));  # [rad]

    return E, nu


def COE_to_SV(COE, mu):
    """
    Convert classical orbital elements to inertial position and velocity (ECI)

    Inputs:
        a    : semi-major axis [km]
        e    : eccentricity [-]
        i    : inclination [°]
        RAAN : right ascension of ascending node [°]
        w    : argument of perigee [°]
        nu   : true anomaly [°]
        mu   : gravitational parameter [km³/s²]

    Outputs:
        R : inertial position vector [km]
        V : inertial velocity vector [km/s]
    """

    p = COE.a * (1.0 - COE.e**2)  # ellipse semi-latus rectum
    rho = p / (1.0 + COE.e * cosd(COE.nu))  # radius in ellipse at ν

    # State vectors in perifocal frame
    r_PQW = np.array([rho * cosd(COE.nu), rho * sind(COE.nu), 0.0])
    v_PQW = np.sqrt(mu / p) * np.array([-sind(COE.nu), COE.e + cosd(COE.nu), 0.0])

    #--- Rotation from perifocal PQW frame to inertial ECI frame
    Rot_PQW_to_ECI = RotZ(COE.RAAN) @ RotX(COE.i) @ RotZ(COE.w)

    R = Rot_PQW_to_ECI @ r_PQW
    V = Rot_PQW_to_ECI @ v_PQW  

    return R, V

# end of COE_to_SV()


def JulianDate(epoch):
    """
    Compute Julian Date from a Python datetime.

    Parameters
    ----------
    epoch : datetime or float
        If datetime, it is interpreted as UTC unless timezone-aware.
        If float, it is assumed to already be a Julian Date.

    Returns
    -------
    JD : float
        Julian Date [days]
    """

    # If epoch is already a Julian Date
    if isinstance(epoch, (int, float)):
        return float(epoch)

    if not isinstance(epoch, datetime):
        raise TypeError("epoch must be a datetime object or a Julian Date float.")

    # Treat naive datetime as UTC
    if epoch.tzinfo is None:
        epoch = epoch.replace(tzinfo=timezone.utc)
    else:
        epoch = epoch.astimezone(timezone.utc)

    year = epoch.year
    month = epoch.month
    day = epoch.day

    hour = epoch.hour
    minute = epoch.minute
    second = epoch.second + epoch.microsecond * 1e-6

    # Fraction of day
    frac_day = (hour + minute / 60.0 + second / 3600.0) / 24.0

    # Gregorian calendar correction
    if month <= 2:
        year -= 1
        month += 12

    A = math.floor(year / 100)
    B = 2 - A + math.floor(A / 4)

    JD = (
        math.floor(365.25 * (year + 4716))
        + math.floor(30.6001 * (month + 1))
        + day
        + frac_day
        + B
        - 1524.5
    )

    return JD

# end of JulianDate


def ComputeGMST(epoch, dut1=0.0):
    """
    Compute Greenwich Mean Sidereal Time (GMST) in degrees using the IAU-82 formula.

    GMST is the Earth rotation angle measured from the mean equinox
    / First Point of Aries to the Greenwich meridian.

    Parameters
    ----------
    epoch : datetime or float
        Epoch of interest. If datetime, it is interpreted as UTC unless timezone-aware.
        If float, it is assumed to be Julian Date UTC/UT1.
    dut1 : float, optional
        UT1 - UTC correction [s]. Default is 0.0.
        For high precision, supply the EOP value of DUT1.

    For precise astrodynamics, the input time should be UT1, not simply UTC. The optional dut1 argument accounts for that through
    UT1=UTC+(UT1−UTC).

    Returns
    -------
    GMST : float
        Greenwich Mean Sidereal Time [deg], in the range [0, 360).
    """

    # Convert UTC epoch to approximate UT1 epoch
    if isinstance(epoch, datetime):
        epoch_ut1 = epoch + timedelta(seconds=dut1)
        JD_UT1 = JulianDate(epoch_ut1)
    else:
        # If a Julian Date is supplied, assume it is already JD_UT1
        JD_UT1 = float(epoch)

    # Julian centuries from J2000.0, using UT1
    T_UT1 = (JD_UT1 - 2451545.0) / 36525.0

    # Vallado 2013, Eq. (3-47), IAU-82
    GMST_s = (67310.54841 + (876600.0 * 3600.0 + 8640184.812866) * T_UT1 + 0.093104 * T_UT1**2 - 6.2e-6 * T_UT1**3)

    # Convert sidereal seconds to degrees and wrap to [0, 360)
    GMST = (360.0 * GMST_s / (24 * 60**2)) % 360.0

    return GMST

# end of ComputeGMST()
    

"""Greenwich Mean Sidereal Time, IAU 1982 formulation."""

# Julian date of the J2000.0 epoch (2000 Jan 1, 12h TT)
JD_J2000 = 2451545.0
DAYS_PER_JULIAN_CENTURY = 36525.0
SECONDS_PER_DAY = 86400.0

def compute_gmst(epoch):
    """Greenwich Mean Sidereal Time in degrees.

    GMST is the rotation angle from the First Point of Aries to the Greenwich
    meridian (longitude 0) at the given epoch. It is the angle needed to
    rotate an ECEF position into ECI, ignoring precession, nutation and
    polar motion.

    Parameters
    ----------
    epoch : datetime, or array_like of datetime
        UT1 epoch. Scalar or sequence.

    Returns
    -------
    float or ndarray
        GMST in degrees, wrapped to [0, 360).

    Notes
    -----
    IAU 1982 expression, Vallado 2013 Section 3.5.2 Eq. (3-47), pp. 184-188,
    originally from the US Naval Observatory's Astronomical Almanac (1984).
    Accurate to roughly 0.1 arcsecond over a century.

    The argument is nominally UT1, not UTC. The difference (DUT1) is bounded
    by 0.9 s, which is up to 3.75 mdeg of Earth rotation, or about 400 m at
    the equator. Adequate for plotting; not for precise orbit determination.

    Examples
    --------
    >>> from datetime import datetime, timezone
    >>> round(float(compute_gmst(datetime(2000, 1, 1, 12, tzinfo=timezone.utc))), 4)
    280.4606
    """
    jd_ut1 = julian_date(epoch)
    t_ut1 = (jd_ut1 - JD_J2000) / DAYS_PER_JULIAN_CENTURY

    gmst_s = (67310.54841
              + (876600.0 * 3600.0 + 8640184.812866) * t_ut1
              + 0.093104 * t_ut1 ** 2
              - 6.2e-6 * t_ut1 ** 3)

    return np.mod(360.0 * gmst_s / SECONDS_PER_DAY, 360.0)


def julian_date(epoch):
    """Julian date from a datetime, or an array of them.

    Naive datetimes are taken to be UTC.
    """
    epoch = np.asarray(epoch)
    scalar = epoch.ndim == 0

    unix = np.array(
        [_to_unix(e) for e in np.atleast_1d(epoch).ravel()],
        dtype=float,
    )
    jd = unix / SECONDS_PER_DAY + 2440587.5

    return jd.item() if scalar else jd.reshape(np.shape(epoch))


def _to_unix(dt):
    """Seconds since 1970-01-01 UTC, assuming UTC for naive datetimes."""
    from datetime import timezone

    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return dt.timestamp()



def Sun_Ephemerides(epoch):
    """Apparent RA and declination of the Sun, geocentric, mean equinox of date.
    Apparent geocentric right ascension and declination of the Sun.

    Low-precision analytic model (Vallado, Algorithm 29, Section 5.1), derived from the Astronomical Almanac's abridged series.
    Accurate to roughly 0.01° in RA/Dec over 1950-2050 -- fine for lighting a globe, sizing an eclipse cone, or a first-cut beta-angle calculation.
    Not adequate where sub-arcsecond astrometry is required; use a JPL DE ephemeris for that.

    Parameters
    ----------
    epoch : datetime or array_like of datetime
        UT1 epoch. Naive datetimes are taken as UTC.
    degrees : bool, default True
        Return angles in degrees. If False, radians.

    Returns
    -------
    ra : float or ndarray
        Right ascension, wrapped to [0, 360) degrees or [0, 2*pi) radians.
    dec : float or ndarray
        Declination, in [-90, +90] degrees or [-pi/2, +pi/2] radians.
    r_au : float or ndarray
        Sun-Earth distance in astronomical units.

    Examples
    --------
    >>> from datetime import datetime, timezone
    >>> ra, dec, r = sun_radec(datetime(2026, 6, 21, 12, tzinfo=timezone.utc))
    >>> bool(23.0 < dec < 23.5)          # near the June solstice
    True
    """
    JD = julian_date(epoch)
    t = (JD - JD_J2000) / DAYS_PER_JULIAN_CENTURY

    # Mean longitude and mean anomaly of the Sun, degrees.
    lam_mean = np.mod(280.460 + 36000.771 * t, 360.0)
    m_sun = np.radians(np.mod(357.5291092 + 35999.05034 * t, 360.0))

    # Ecliptic longitude, corrected for the equation of centre.
    lam_ecl = np.radians(lam_mean + 1.914666471 * np.sin(m_sun) + 0.019994643 * np.sin(2.0 * m_sun))

    # Obliquity of the ecliptic, mean of date.
    eps = np.radians(23.439291 - 0.0130042 * t)

    # Ecliptic -> equatorial. The Sun's ecliptic latitude is taken as zero,
    # which costs at most about 1 arcsecond.
    sin_lam = np.sin(lam_ecl)
    RA = np.mod(np.arctan2(np.cos(eps) * sin_lam, np.cos(lam_ecl)), 2.0 * np.pi)
    Dec = np.arcsin(np.sin(eps) * sin_lam)

    # Sun-Earth distance, AU.
    R_AU = 1.000140612 - 0.016708617 * np.cos(m_sun) - 0.000139589 * np.cos(2.0 * m_sun)
    Dist_km = R_AU / Earth.AU

    return DEG(RA), DEG(Dec), R_AU


def sun_vector_eci(epoch, km=True):
    """Geocentric position vector of the Sun in ECI (mean equator of date).

    Returns a (3,) array for a scalar epoch, or (3, n) for a sequence,
    matching the row-stacked convention used elsewhere in this project.
    """
    ra, dec, r_au = sun_radec(epoch, degrees=False)
    r = r_au * AU_KM if km else r_au

    return np.array([r * np.cos(dec) * np.cos(ra),
                     r * np.cos(dec) * np.sin(ra),
                     r * np.sin(dec)])

def ECEF_to_AER(RV_ECEF, lat_geod, lon):
    """
    Converts an ECEF vector into Azimuth, Elevation, Range and SEZ.

    Inputs
    ------
    RV_ECEF : ndarray
        Position vector in ECEF frame [km].
    lat_geod : float
        Geodetic latitude [deg].
    lon : float
        Longitude [deg].

    Returns
    -------
    Azim : float
        Azimuth [deg].
    Elev : float
        Elevation [deg].
    Range : float
        Distance [km].
    SEZ : ndarray
        Vector expressed in the local South-East-Zenith frame.
    """

    # ECEF -> SEZ transformation
    SEZ = RotY(lat_geod-90) @ RotZ(-lon) @ RV_ECEF

    # Range
    Range = np.linalg.norm(SEZ)

    # Azimuth
    Azim = atan2d(SEZ[1], -SEZ[0]) % 360

    # Elevation
    Elev = asind(SEZ[2] / Range)

    return Azim, Elev, Range, SEZ


def geod_to_pos(phi_geod, lam=0.0, h_geod=0.0, a=Earth.r1_km, b=Earth.r2_km):
    """
    Convert geodetic latitude to Cartesian position vector and
    geocentric latitude at a given geodetic height above an ellipsoid.

    Parameters
    ----------
    phi_geod : float or array_like
        Geodetic latitude [deg].
    lam : float or array_like, optional
        Longitude [deg], east-positive.
        Default = 0 deg.
    h_geod : float or array_like, optional
        Geodetic height above the reference ellipsoid [km].
        Default = 0 km.
    a : float, optional
        Semi-major (equatorial) radius of the ellipsoid [km].
        Default = 6378.137 km (WGS84).
    b : float, optional
        Semi-minor (polar) radius of the ellipsoid [km].
        Default = 6356.75231424518 km (WGS84).

    Returns
    -------
    R_sat : ndarray
        Cartesian position vector [x, y, z] [km].
        Shape is (3,) for scalar input and (3, N) for vector input.
    phi_geoc : float or ndarray
        Geocentric latitude [deg].

    Notes
    -----
    Based on the same equations as the MATLAB geod_to_pos function.
    Angles are supplied in degrees.
    """
    phi_geod = np.asarray(phi_geod, dtype=float)
    lam = np.asarray(lam, dtype=float)
    h_geod = np.asarray(h_geod, dtype=float)

    if np.any(lam > 360.0):
        raise ValueError("geod_to_pos: second argument lambda is out of bounds.")

    # Prime vertical radius of curvature (only quantity needed here)
    R_N, _, _, _, _ = ellipsoid_radii(phi_geod, a, b)

    sin_phi = sind(phi_geod)
    cos_phi = cosd(phi_geod)
    cos_lam = cosd(lam)
    sin_lam = sind(lam)

    # Cartesian coordinates:
    #
    # x = (R_N + h) cos(phi) cos(lambda)
    # y = (R_N + h) cos(phi) sin(lambda)
    # z = ((b/a)^2 R_N + h) sin(phi)
    x = (R_N + h_geod) * cos_phi * cos_lam
    y = (R_N + h_geod) * cos_phi * sin_lam
    z = ((b / a) ** 2 * R_N + h_geod) * sin_phi

    R_sat = np.array([x, y, z])

    # Geocentric latitude
    phi_geoc = atand(
        tand(phi_geod) * (R_N * (b / a) ** 2 + h_geod) / (R_N + h_geod)
    )

    return R_sat, phi_geoc


def geoc_to_pos(phi_geoc, lam=0.0, h_geoc=0.0, a=Earth.r1_km, b=Earth.r2_km):
    """
    Convert geocentric latitude and longitude to Cartesian position vector.
 
    Parameters
    ----------
    phi_geoc : float or array_like
        Geocentric latitude [deg].
 
    lam : float or array_like, optional
        Longitude [deg], east-positive.
        Default = 0 deg.
 
    h_geoc : float or array_like, optional
        Altitude above the ellipsoid measured along the
        geocentric radial direction [km].
        Default = 0 km.
 
    a : float, optional
        Semi-major (equatorial) radius of the ellipsoid [km].
        Default = 6378.137 km.
 
    b : float, optional
        Semi-minor (polar) radius of the ellipsoid [km].
        Default = 6356.75231424518 km.
 
    Returns
    -------
    R : ndarray
        Cartesian position vector [x, y, z] [km].
        Shape is (3,) for scalar input and (3, N) for vector input.
 
    r_geoc : float or ndarray
        Distance from the center of the ellipsoid to its surface
        along the specified geocentric latitude [km].
 
    Notes
    -----
    The altitude h_geoc is measured along the geocentric radial
    direction, not along the ellipsoid normal.
 
    For a spherical body (a == b), this reduces to the usual
    spherical-to-Cartesian transformation with radius a + h_geoc.
    """
 
    # Convert inputs to NumPy arrays for vectorized calculations
    phi_geoc = np.asarray(phi_geoc, dtype=float)
    lam      = np.asarray(lam, dtype=float)
    h_geoc   = np.asarray(h_geoc, dtype=float)
 
    # --------------------------------------------------------------
    # Radius of ellipsoid along the geocentric direction
    #
    # r_geoc = a / sqrt(1 + ((a/b)^2 - 1) sin^2(phi_geoc))
    r_geoc = a / np.sqrt(1.0 + ((a / b) ** 2 - 1.0) * sind(phi_geoc) ** 2)
 
    # --- Cartesian coordinates via the generic spherical -> Cartesian
    # transform (radius = r_geoc + h_geoc, along phi_geoc/lam)
    rho = r_geoc + h_geoc
 
    R = Sph2Cart(lam, phi_geoc, rho)
 
    return R, r_geoc


def pos_to_geoc(R_ECEF):
    """
    Obtain the geocentric quantities used throughout `pos_to_geod` from a
    set of ECEF Cartesian coordinates: geocentric longitude, latitude,
    radius, and the auxiliary cylindrical radius / height.
 
    Parameters
    ----------
    R_ECEF : (3,) or (3, N) array_like
        ECEF Cartesian coordinates [x; y; z] (any consistent length unit).
 
    Returns
    -------
    r : ndarray
        Distance from the Z axis, sqrt(x^2 + y^2).
    z : ndarray
        Z coordinate.
    Lon : ndarray
        Geocentric longitude [deg] (== geodetic longitude).
    Lat : ndarray
        Geocentric latitude [deg].
    Rho : ndarray
        Geocentric radius.
    """
    R_ECEF = np.asarray(R_ECEF, dtype=float)
    if R_ECEF.ndim == 1:
        R_ECEF = R_ECEF.reshape(3, 1)
 
    Lon, Lat, Rho = Cart2Sph(R_ECEF)
    Lon = ((Lon + 180) % 360) - 180
    r = np.hypot(R_ECEF[0, :], R_ECEF[1, :])   # distance from Z axis
    z = R_ECEF[2, :]
 
    return r, z, Lon, Lat, Rho
 
 
def pos_to_geod(R_ECEF, a=Earth.r1_km, b=Earth.r2_km, method=2, tol=1e-14):
    """
    (ECEF) Cartesian position -> geodetic coordinates: (x, y, z) -> (phi, h).
 
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
 
    See Also
    --------
    pos_to_geoc : obtains Lon, Lat, Rho, r, z from R_ECEF (used internally).
 
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
 
    # --- 0) Geocentric longitude, latitude, radius, and auxiliary r, z ----
    r, z, Lon, Lat, Rho = pos_to_geoc(R_ECEF)
 
    e2 = 1 - (b / a) ** 2     # first eccentricity squared
    e_2 = (a / b) ** 2 - 1    # second eccentricity squared == e2 / (1 - e2)
 
    # --- 1) Geodetic latitude and altitude ---------------------------------
    h_geod = None
 
    if method in (0, 1) and a != 6378.137:
        raise ValueError(
            "pos_to_geod: Methods 0 & 1 are for Earth only (not applicable "
            "for normalized radii) => use method 2 instead."
        )
 
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
        h_geod = None  # completed generically below (see step 2)
 
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
        raise ValueError(f"pos_to_geod: Unknown or unimplemented method #{method}")
 
    # --- 2) Finalize the geodetic latitude/height pair ---------------------
    #
    # Ensures phi_geod is a 1-D array, and - for methods that only produce
    # the latitude directly (h_geod is None, e.g. method 11) - completes
    # the geodetic height generically from the prime vertical radius of
    # curvature R_N(phi_geod):
    #
    #   R_N    = a / sqrt(1 - e^2 sin^2(phi_geod))
    #   h_geod = r / cos(phi_geod) - R_N
    phi_geod = np.atleast_1d(np.asarray(phi_geod, dtype=float))
 
    R_N = a / np.sqrt(1 - e2 * sind(phi_geod) ** 2)  # R_N(phi_geod)
 
    if h_geod is None:
        h_geod = r / cosd(phi_geod) - R_N
    else:
        h_geod = np.atleast_1d(np.asarray(h_geod, dtype=float))
 
    return phi_geod, h_geod, Lon, Lat, Rho


def geod_to_geoc(phi_geod, a=Earth.r1_km, b=Earth.r2_km, h0=0.0):
    """
    Compute the geocentric latitude and Vermeille's geodetic-direction altitude for a given geodetic latitude.
 
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
    x_sat, z_sat = line_circle_intersection(m, b_, rho)
 
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
