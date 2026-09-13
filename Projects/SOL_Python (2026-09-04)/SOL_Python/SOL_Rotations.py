# ╭───────────────────╮
# │  SOL_Orbit_2B.py  │ Two-Body unperturbed orbital propagation
# ╰───────────────────╯
#  © SpaceOrbitLAB 
"""
# (Header to include here)
# (c) SpaceOrbitLAB

(main).py

SOL_Tools/
  ├─ AstroConstants.py
  ├─ Orbit_tools.py
  ├─ Math_tools.py
  └─ Plot_tools.py

  Visual Studio Code: Color Theme: Light Modern
"""

# General Python tools
import numpy as np
import pyvista as pv  # https://docs.pyvista.org/getting-started/
                      # https://docs.pyvista.org/api/plotting/_autosummary/pyvista.plotter.add_title
import os
import sys

# SOL specific tools
from SOL_Tools.AstroConstants import *  # all variales
from SOL_Tools.Orbit_tools import *  # OrbitPropagation_2BN, OrbitPropagation_2BK, COE_to_SV
from SOL_Tools.Graph_tools import *  # DrawEarth3D, DrawAngularMomentumVector
from SOL_Tools.Math_tools import * 

from datetime import datetime, timezone
from astropy.time import Time


if __name__ ==  "__main__":

    if False:  # (debugging)
        print("Earth.mu = ", Earth.mu)
        print("Earth.ecc = ", Earth.ecc)
        print('')
        print(os.path.dirname(__file__))
    print('')

    # Initialize the 3D PyVista Scene
    pl = pv.Plotter(lighting="none", window_size = (1200, 900))  # kill the default three-light kit, we supply our own
    pl.render_window.SetPosition(100, 50)  # Position on compute screen, from top left origin
    # pl.set_background("black")  → nice, but impractical for course notes and printing

    front_light = pv.Light(light_type = "headlight", intensity = 0.9)
    pl.add_light(front_light)

    #--- Defines Epoch
    #Date_UTC = datetime.now(timezone.utc)
    # Reference epoch:
    Date_UTC = datetime(2026,9,3, 12,0,0)  # → GMST = 162.637288334314° at this date; At Noon, the Sun will be at Greenwich meridian
    if False:
        Date_UTC -= timedelta(hours=24 * -71.2778252/360)  # → brings the Sun exactly in Quebec's meridian
        print(Date_UTC)
        
    # Computes the Greenwich Mean Sidereal Time corresponding to this Epoch
    GMST = ComputeGMST(Date_UTC)  # [°]
    GMST2 = compute_gmst(Date_UTC)

    str = f'Epoch = {Date_UTC} → {Date_UTC.strftime('%Y-%m-%d %H:%M:%S')} UTC → GMST = {GMST:.4f}° → GMST = {GMST2:.4f}°'
    print(str)
    txt = pl.add_text(str, position = "upper_left", font_size = 12, font = 'courier')
    #txt = pl.add_text(str, position = "upper_left", font_size = 8, font_file='C:/Windows/Fonts/Calibri.ttf')
    #txt.prop.bold = False

    t = Time.now()  # Astropy: Computes all key time components
    #print(t.mjd, t.jd, t.iso, t.tt, t.tai)

    Sun_RA, Sun_Dec, Sun_Lon, Sun_dist_AU = sun_radec(Date_UTC)
    #TODO: Bring this into SunEphemeris()
    class Sun:
        RA = Sun_RA    # [°] Sun's Right Ascension
        Dec = Sun_Dec  # [°] Sun's Declination -𝜖 < δ < +𝜖
        Lon = Sun_RA - GMST  # [°] Sun's Longitude
        Dist_AU = Sun_dist_AU
        Dist_km = Sun_dist_AU * Earth.AU
        #RA = -71.3571   # [°] Sun's Right Ascension
        #Dec = 20  # [°] Sun's Declination -𝜖 < δ < +𝜖
    print(f'Sun RA = {Sun.RA:.4f}, Dec = {Sun.Dec:.4f}, Lon = {Sun.Lon:.4f}°, Dist = {Sun.Dist_AU:.4f} AU = {Sun.Dist_km:.0f} km\n')
    print(' GMST + Sun.Lon = ', GMST + Sun.Lon)  # ← lon ≡ RA - GMST


    ViewMode = 'ECI'    # Earth-Centered, Intertial
    ViewMode = 'ECEF'   # Earth-Centered, Earth-Fixed

    pl = DrawEarth3D(pl, ViewMode)

    # initial view, will be re-adjusted later
    pl.camera.azimuth = -35  # [°]
    pl.camera.elevation = 40   # [°]
    pl.camera.zoom(1.2)  # useful to zoom/unzoom


    # (1) Starting point at λ = 0°, φ = 0°
    xyz = Sph2Cart(0, 0, Earth.r1_km)
    if True:
        pl.add_points(np.column_stack(xyz), color = 'red', point_size = 10, render_points_as_spheres=True, lighting=False)
    else:
        point = pv.Sphere(radius = Earth.r1_km / 50, center = np.column_stack(xyz))
        pl.add_mesh(point, color = 'red', show_edges = False, ambient = 1)

    # Greenwidh Latitude (Parallel 51.5°)
    if False:
        Observer = DefineObserverLocation('Greenwich')  # TODO
        phi_Greenwich = Observer.Lat;
    else:
        phi_Greenwich = 51.4934  # [°]

    # (2) Computes point location up to latitude φ = +51.5° (Greenwich latitude)
    Lat = phi_Greenwich # [°N]
    xyz_ = RotY(-Lat) @ xyz;  # ← requires a negative sign, as latitudes are measured clockwise, while the RotY rotation
                              #   operator is defined counter-clockwise
    # == Sph2Cart(0, Lat, Earth.r1_km)
    pl.add_points(np.column_stack(xyz_), color = 'red', point_size = 10, render_points_as_spheres=True, lighting=False)

    # (3) Draws arc line up to Greenwich
    phi = np.linspace(0, Lat, 100)
    arc = Sph2Cart(np.zeros_like(phi), phi, Earth.r1_km)
    if False:
        print('arc shape =', arc.shape) 
        print('arc[start]  =', arc[:,0])
        print('arc[end]=', arc[:,-1])
    pl.add_mesh(pv.lines_from_points(arc.T), color = 'red', line_width = 4)  

    # (4a) Draw Earth Parallel at Greenwich latitude
    theta = np.linspace(0.0, 360, 100)  # [°]
    C = EXPAND_FACTOR * Sph2Cart(theta, phi_Greenwich * np.ones_like(theta), Earth.r1_km)
    pl.add_mesh(pv.lines_from_points(C.T), color = 'cyan', line_width = 2.0)
    # (4b) Using Small Circle
    lat, lon, xyz = SmallCircle(0, 90, 90 - phi_Greenwich)  # unitary SC
    xyz_ = EXPAND_FACTOR * Earth.r1_km * xyz  # scales to graph
    if False:
        print('φ = ', lat[0], lat[-1])
        print('λ = ', lon[0], lon[-1])
        print('ρ[0] = ', xyz_[:,1])
    pl.add_mesh(pv.lines_from_points(xyz_.T), color = 'green', line_width = 1.5)

    # (5) Computes point location westward to longitude λ = 71.2°W (Quebec City longitude)
    Lon = -71.2778252  # [°W]
    Lat =  46.78095    # [°N]
    xyz = Sph2Cart(0, 0, Earth.r1_km)
    xyz_ = RotZ(+Lon) @ xyz   # ← positive here, as longitude is computed
                              #   counter-clockwise, the same as RotZ operator
    pl.add_points(np.column_stack(xyz_), color = 'lime', point_size = 10, render_points_as_spheres=True, lighting=False)

    # (6) Draws arc from Equator
    phi = np.linspace(0, Lat, 100);
    arc = RotZ(+Lon) @ Sph2Cart(np.zeros_like(phi), phi, Earth.r1_km)
    pl.add_mesh(pv.lines_from_points(arc.T), color = 'lime_green', line_width = 4)
    # (7) Draws last point of arc
    pl.add_points(arc[:,-1], color = 'lime', point_size = 10, render_points_as_spheres=True, lighting=False)

    #--- Draws Sun vector pointing to Earth Center and add lighting in this direction
    if ViewMode ==  'ECEF': 

        # (8) Draw Sun vector pointing to Earth center, knowing Sun longitude λ_⊙ and Sun declination δ_⊙
        DrawSunVector(pl, Sun)

        str = f"Sun RA = {Sun.RA:.4f}°, Dec = {Sun.Dec:.4f}°, Lon = {Sun.Lon:.4f}°"
        pl.add_text(str, position = "lower_left", font_size = 8, font='courier')

        # (9) Draws Ecliptic line knowing Earth obliquity 𝜖 and GMST angle at epoch
        DrawEclipticLine(pl, GMST, Earth)


    #--- Draws the Great Circle arc between Miami & Greenwich
    city1 = 'Miami, USA'
    lat1 =  25.7636  # [°]
    lon1 = -80.1296  # [°]

    city2 = 'Greenwich'
    lat2 = phi_Greenwich  # [°]
    lon2 = 0.0  # [°]

    GC = GreatCircle2(lat1, lon1, lat2, lon2)
    GC_path = EXPAND_FACTOR * Earth.r1_km * GC.xyz
    pl.add_mesh(pv.lines_from_points(GC_path.T), color = 'blue', line_width = 2)  

    # -------------------------------------------------------------------------------
    # -------------------------------- Homework -------------------------------------
    #--------------------------------------------------------------------------------

    # Draw South/East/Z (up) plane at Laval University location λ = -71.357°, φ = 46.774
    Lon=-71.357 # [°]
    Lat=46.774 # [°]
    R=RotZ(Lon) @ RotY(-Lat)

    xyz=EXPAND_FACTOR*Sph2Cart(Lon,Lat,Earth.r1_km)
    pl.add_points(np.column_stack(xyz), color = 'red', point_size = 15, render_points_as_spheres=True, lighting=False)

    L=0.5*Earth.r1_km
    xyz_=np.column_stack(xyz)[0,:]
    S=R @ Sph2Cart(0,-90, L)
    E=R @ Sph2Cart(90, 0, L)
    Z=R @ Sph2Cart(0, 0, L)

    S_=np.column_stack(S)[0,:]
    E_=np.column_stack(E)[0,:]
    Z_=np.column_stack(Z)[0,:]

    South_line=np.array([xyz_, xyz_+S_])
    pl.add_mesh(pv.lines_from_points(South_line),color='blue',line_width=4)

    East_line=np.array([xyz_, xyz_+E_])
    pl.add_mesh(pv.lines_from_points(East_line),color='lime',line_width=4)

    Z_line=np.array([xyz_, xyz_+Z_])
    pl.add_mesh(pv.lines_from_points(Z_line),color='red',line_width=4)

    pl.add_mesh(pv.Plane(
            center=xyz_,
            direction=Z_,
            i_size=Earth.r1_km,
            j_size=Earth.r1_km),
            color='lightgray',
            opacity=0.25,
            show_edges=True)





    print('Done.')  # appears when all computations done

    #--- Shows the final graphic
    pl.show(title = '© SpaceOrbitLAB')

    print('Exit.')  # appears when closing the window

# end __main__