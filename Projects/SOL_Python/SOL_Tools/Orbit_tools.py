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

    t_span = np.array([time_s[0], time_s[-1]])

    # Solve using RK45 (Python's ode45 equivalent)
    sol = solve_ivp(Rates_Newton, t_span, SV, t_eval=time_s, method='RK45', rtol=1e-10, atol=1e-12)

    # (here ".y" does not mean the 2nd 3D dimension)
    ECI_pos_2BN = sol.y[0:3, :]  # first 3 elements
    ECI_vel_2BN = sol.y[3:6, :]  # last 3 elements

    return ECI_pos_2BN, ECI_vel_2BN

def Rates_Newton(t, z):

    r_vec = z[0:3]
    v_vec = z[3:6]

    r = norm(r_vec)

    dzdt = np.zeros(6)
    dzdt[0:3] = v_vec
    dzdt[3:6] = -Earth.mu * r_vec / r**3

    return dzdt


def OrbitPropagation_2BK(a, e, i, RAAN, w, M0, t_array, mu):
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


def COE_to_SV(a, e, i, RAAN, w, nu, mu):
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

    p = a * (1.0 - e**2)  # ellipse semi-latus rectum
    rho = p / (1.0 + e * cosd(nu))  # radius in ellipse at ν

    # State vectors in perifocal frame
    r_PQW = np.array([rho * cosd(nu), rho * sind(nu), 0.0])
    v_PQW = np.sqrt(mu / p) * np.array([-sind(nu), e + cosd(nu), 0.0])

    #--- Rotation from perifocal PQW frame to inertial ECI frame
    Rot_PQW_to_ECI = RotZ(RAAN) @ RotX(i) @ RotZ(w)

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



"""Apparent geocentric right ascension and declination of the Sun.

Low-precision analytic model (Vallado, Algorithm 29, Section 5.1), derived
from the Astronomical Almanac's abridged series. Accurate to roughly
0.01 degrees in RA/Dec over 1950-2050 -- fine for lighting a globe, sizing
an eclipse cone, or a first-cut beta-angle calculation. Not adequate where
sub-arcsecond astrometry is required; use a JPL DE ephemeris for that.
"""


JD_J2000 = 2451545.0
DAYS_PER_JULIAN_CENTURY = 36525.0
SECONDS_PER_DAY = 86400.0
AU_KM = 149597870.700


def sun_radec(epoch, degrees=True):
    """Apparent RA and declination of the Sun, geocentric, mean equinox of date.

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
    jd = julian_date(epoch)
    t = (jd - JD_J2000) / DAYS_PER_JULIAN_CENTURY

    # Mean longitude and mean anomaly of the Sun, degrees.
    lam_mean = np.mod(280.460 + 36000.771 * t, 360.0)
    m_sun = np.radians(np.mod(357.5291092 + 35999.05034 * t, 360.0))

    # Ecliptic longitude, corrected for the equation of centre.
    lam_ecl = np.radians(
        lam_mean
        + 1.914666471 * np.sin(m_sun)
        + 0.019994643 * np.sin(2.0 * m_sun)
    )

    # Sun-Earth distance, AU.
    r_au = (1.000140612
            - 0.016708617 * np.cos(m_sun)
            - 0.000139589 * np.cos(2.0 * m_sun))

    # Obliquity of the ecliptic, mean of date.
    eps = np.radians(23.439291 - 0.0130042 * t)

    # Ecliptic -> equatorial. The Sun's ecliptic latitude is taken as zero,
    # which costs at most about 1 arcsecond.
    sin_lam = np.sin(lam_ecl)
    ra = np.mod(np.arctan2(np.cos(eps) * sin_lam, np.cos(lam_ecl)), 2.0 * np.pi)
    dec = np.arcsin(np.sin(eps) * sin_lam)

    if degrees:
        return np.degrees(ra), np.degrees(dec), 0, r_au
    return ra, dec, 0, r_au


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
