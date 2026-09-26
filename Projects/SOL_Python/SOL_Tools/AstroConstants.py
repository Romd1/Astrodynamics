# ┌───────────────────────────────┐
# │  SOL_Tools\AstroConstants.py  │ 
# └───────────────────────────────┘
#  © SpaceOrbitLAB 

#TODO: → Progressively import parameters as needed from AstroConstants.m & DefineConstants.m

import numpy as np
#from typing import Final

# Basic global constants
# MY_PI: Final = 3.14159
#DATABASE_URL: Final = "sqlite:///my_database.db"
#MAX_LOGIN_ATTEMPTS: Final = 5

EXPAND_FACTOR = 1.005  # factor to make radius vectors slightly larger so that he line can be seen at the surface of the Earth

JD_J2000 = 2451545.0
DAYS_PER_JULIAN_CENTURY = 36525.0
SECONDS_PER_DAY = 86400.0
AU_KM = 149597870.700

SUN_DISTANCE = 1.0e9  # [km] Only the direction matters for a directional light, but a large value keeps the geometry unambiguous.


class Earth:
    
    mu = 398600.4418  # [km³/s²]  Earth gravitational parameter

    r1_km = 6378.137  # [km]  Earth equatorial radius ; Radius (at sea level) at equator
    r2_km = 6356.75231424518  # [km]  Earth polar radius (at sea level)

    if False: # artificially shrinks Earth's polar radius for geometry validation
        print('*** Artificially scaling down Earth’s polar radius ***')
        r2_km *= 0.8
    sf = r2_km / r1_km  # scale factor

    # Eccentricity of the Earth's oblate ellipsoidal shape  e ≡ √(1 – (b/a)²) = sin(θ)
    ecc = np.sqrt(1 - (r2_km / r1_km)**2)  # [n.u.] 

    # TODO: implement EarthObliquity.py
    Obliquity = 23.43582299  # [°]

    # The astronomical unit (AU) is a unit of length derived from the Earth’s orbit. It is the average distance the
    # Earth gets from the Sun on the long axis of the ellipse. Its definition is: the length of the semi-major axis
    # of the Earth’s elliptical orbit around the Sun.
    AU = 149597870.700;  # [km] Astronomical Unit; "average" distance from Earth to the Sun (definition)
                        # .700 used by Horizons
                        # → distance between the Sun and the Earth's orbital center (not to surface)
                        # == (Sun.mu * (24 * 60^2 * Earth.Sidereal_year / (2*pi))^2)^(1/3)
    a = 149598023;      # [km] Semi-major axis (www.glyphweb.com/esky/concepts/semimajoraxis.html)
    # TODO: implement EarthOrbitEccentricity
    e = 0.016698  # [n.u.] Earth's orbit eccentricity e ≈ 0.016698 in 2025
