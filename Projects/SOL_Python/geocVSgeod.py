"""
geod_geoc_compare.py

Script comparing, for a range of geodetic latitudes phi (0 to 90 deg, by
15 deg steps), the corresponding geocentric latitude phi_geoc, the
altitude discrepancy between the geodetic-normal and geocentric-radius
directions, and North-South / East-West arc-length distances -- with a
cross-check of meridian-arc length via the analytical elliptic-integral
formula.

Translated from the original MATLAB script. MATLAB Mapping Toolbox calls
(`wgs84Ellipsoid`, `distance`, `meridianarc`) were only used inside
disabled (`if 0`) validation blocks in the original script and are not
reproduced here; the analytical elliptic-integral cross-check
(`segment_dist1`, using `ellipticE`) IS translated, using
`scipy.special.ellipeinc`.
"""

import numpy as np
from SOL_Tools.AstroConstants import *  # all variales
from SOL_Tools.Math_tools import * 

np.seterr(divide="ignore", invalid="ignore")
from scipy.special import ellipeinc

from geod_to_geoc import geod_to_geoc
from geod_to_pos import geod_to_pos
from geoc_to_pos import geoc_to_pos


def main():
    a = Earth.r1_km
    b = Earth.r2_km
    e2 = 1 - (b / a) ** 2  # == Earth.ecc**2

    # NOTE: MATLAB's wgs84Ellipsoid('kilometer') / distance() / meridianarc()
    # (Mapping Toolbox) are only used below in disabled (`if 0`) validation
    # blocks; they are not reproduced here (kept as comments for reference).

    phi_geod = np.sort(np.arange(0, 91, 15, dtype=float))
    # phi_geod = np.sort(np.arange(0, 91, 5, dtype=float))

    lambda_ECEF = 0.0
    h0 = 0.0
    # h0 = 1000.0
    Rho = a + h0

    n = len(phi_geod)
    phi_geoc0 = np.zeros(n)
    segment_dist = np.full(n, np.nan)
    segment_dist1 = np.full(n, np.nan)
    segment_dist2 = np.full(n, np.nan)  # not computed (Mapping Toolbox 'distance' unavailable)

    print(f"h0 = {h0:g} km")
    print("  phi [deg] │   phi_geoc [deg]  │  Delta = phi-phi_geoc  │ Delta [arcmin] │ "
          "Delta [arcsec] │    Delta h_perp  │ N-S dist │  arc [km]")
    print("-" * 118)

    for i in range(n):
        phi_geoc_i, h_geod_i = geod_to_geoc(phi_geod[i], a, b, h0)
        phi_geoc_i = phi_geoc_i[0]
        h_geod_i = h_geod_i[0]

        phi_geoc0[i] = atand(tand(phi_geod[i]) * (b / a) ** 2)
        _, R_geoc = geoc_to_pos(phi_geoc0[i], lambda_ECEF, 0.0, a, b)
        R_geoc = R_geoc[0]
        h_geoc = Rho - R_geoc

        delta_lat = phi_geod[i] - phi_geoc_i  # Delta_lat = phi - phi_geoc [deg]
        delta_h = h_geoc - h_geod_i           # h_phi - h_phi_geoc [km] (cm-order)

        # Here the computed radii R_e, R_N and R_M are independent of h
        # (R_sub would not be)
        _, _, _, R_N, R_M, *_ = geod_to_pos(phi_geod[i], lambda_ECEF, 0.0, a, b)
        R_N = R_N[0]
        R_M = R_M[0]

        NS_dist_km = R_M * RAD(delta_lat)  # arc length

        # --- Validation block (disabled in the original MATLAB script) ---
        # if 0:
        #     md = meridianarc(RAD(phi_geoc_i), RAD(phi_geod[i]), ellipsoid)  # [m]
        #     arclen, az = distance(phi_geoc_i, 0, phi_geod[i], 0, ellipsoid)  # @ lambda = 0 deg
        #     print([md / 1000, arclen / 1000, NS_dist_km])

        if i > 0:
            # R_M is computed independent of h
            phi_mid = (phi_geod[i] + phi_geod[i - 1]) / 2
            _, _, _, _, R_M_mid, *_ = geod_to_pos(phi_mid, lambda_ECEF, 0.0, a, b)
            R_M_mid = R_M_mid[0]
            segment_dist[i] = R_M_mid * RAD(phi_geod[i] - phi_geod[i - 1])  # arc length

            # Analytical solution (incomplete elliptic integral of the 2nd kind,
            # parameter convention m = e2, matching MATLAB's ellipticE(phi, m))
            segment_dist1[i] = a * (
                ellipeinc(RAD(phi_geoc0[i]), e2) - ellipeinc(RAD(phi_geoc0[i - 1]), e2)
            )

            # MATLAB Mapping Toolbox function (not available/translated here):
            # segment_dist2[i] = distance(phi_geod[i], 0, phi_geod[i-1], 0, ellipsoid)

        delta_lon = delta_lat
        EW_dist_km = R_N * cosd(phi_geod[i]) * RAD(delta_lon)  # ds = (R_N cos(phi)) dlambda

        # --- Validation block (disabled in the original MATLAB script) ---
        # if 0:
        #     arclen, az = distance(phi_geod[i], 0, phi_geod[i], delta_lon, ellipsoid)
        #     print([arclen / 1000, EW_dist_km])

        seg_str = f"{segment_dist[i]:8.1f}" if not np.isnan(segment_dist[i]) else "     n/a"

        print(
            f"{phi_geod[i]:7.5g} deg │ {phi_geoc_i:7.4f} deg │ {delta_lat:8.6f} deg │ "
            f"{delta_lat * 60:8.3f} '  │ {delta_lat * 3600:8.3f} \"  │ "
            f"{delta_h * 1000:5.2f} m │ {NS_dist_km:5.2f} km │ {seg_str} km"
        )

    print("-" * 118)
    print(f"{'':>90}1/4 circumference: sum = {np.nansum(segment_dist[1:]):.1f} km")

    phi_geod = 5.3606 # [°]
    lam = 138.7274 # [°]
    H_MSL = 3776 # [m]

    R_sat, phi_geoc, R_e, R_N, R_M, x0, z0 = geod_to_pos(phi_geod, lam, H_MSL/1000)

    R_top = R_sat.ravel()

    x, y, z = R_top[0], R_top[1], R_top[2]
    print(f"x: {x:.4f} km")
    print(f"y: {y:.4f} km")
    print(f"z: {z:.4f} km")

    # z = 592.2501 km != H_MSL = 3.7 km because z represents the vertical elevation relative to the equatorial plane in the ECEF frame of reference, which is very different from the altitude H_MSL which is the height above the ellipsoid


if __name__ == "__main__":
    main()