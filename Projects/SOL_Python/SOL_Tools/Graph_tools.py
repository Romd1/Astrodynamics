# ┌────────────────────────────┐
# │  SOL_Tools\Graph_tools.py  │ 
# └────────────────────────────┘
#  © SpaceOrbitLAB 

import numpy as np
import pyvista as pv
from pyvista import examples
import os
#from PIL import Image  # Pillow module: the standard and most popular third-party image processing library for Python
#import matplotlib.pyplot as plt

# Forces VTK to load the Matplotlib LaTeX parser backend
import vtkmodules.vtkRenderingMatplotlib  
import vtkmodules.vtkRenderingFreeType

from SOL_Tools.AstroConstants import EXPAND_FACTOR, Earth
from SOL_Tools.Math_tools import *  
from SOL_Tools.Orbit_tools import geod_to_pos


def DrawEarth3D(pl, view_mode):
# Textured 3D Earth using PyVista

    # Earth equatorial radius
    RE = Earth.r1_km  # [km]

    #--- Arrow locations and directions

    directions = np.eye(3)
    for d in directions:
        # print('d = ', d)
        axis = pv.Arrow(start = (0, 0, 0), direction = d,
                        tip_length = 0.25/5, tip_radius = 0.08/7, shaft_radius = 0.025/7,
                        tip_resolution = 32, shaft_resolution = 32,       # (could use 16)
                        scale = 1.5 * RE)
        pl.add_mesh(axis, color = 'gray', smooth_shading = True, lighting = False)

    #--- Create the Earth sphere
    if False:
        Earth_sphere = pv.Sphere(
            radius = RE,
            #theta_resolution = 72,    # sufficient for smooth rendering in ECI
            #phi_resolution = 36,
            theta_resolution = 720,   # required for nice Earth mapping 
            phi_resolution = 360,
            #start_theta = 270.001,
            #end_theta = 270.0,
        )

        #--- Manually assign spherical texture coordinates
        # This gives better control for latitude-longitude Earth maps.
        # The texture image is assumed to be:
        #   horizontal axis = longitude
        #   vertical axis   = latitude
        #   top             = North Pole

        points = Earth_sphere.points

        x = points[:, 0]
        y = points[:, 1]
        z = points[:, 2]

        r = np.sqrt(x**2 + y**2 + z**2)

        # Texture coordinates:
        # u: longitude-like coordinate
        u = 0.5 + (np.arctan2(-x, y) + np.pi/2) / (2.0 * np.pi)  # ← the +π/2 here rotates the map to have λ = 0° at the Greenwich Prime Meridian
        # v: latitude-like coordinate
        v = 0.5 + np.arcsin(z / r) / np.pi

        # Stores texture coordinates in the mesh
        Earth_sphere.active_texture_coordinates = np.column_stack((u, v))

        # print(Earth_texture.dimensions)
        # TODO: option downsampling for faster plotting 
    else:
        Earth_sphere = examples.planets.load_earth(radius = Earth.r1_km, lat_resolution = 180, lon_resolution = 360)
        # Adjusts the Earth's oblateness
        Earth_sphere.scale((1.0, 1.0, Earth.sf), inplace=True)

        # Rotates map to have Greenwich located at Prime Meridian (PRIME_MERIDIAN_OFFSET_DEG = 180)
        Earth_sphere.rotate_z(180, point = (0, 0, 0), inplace = True)

    #--- Earth texture

    # Local texture image file
    if view_mode ==  'ECI':
        # Draws a generic rotating globe to show that in ECI no ground features (i.e. continents) can be displayed, 
        # as each point of the orbit corresponds to a different Earth rotation angle
        texture_file = os.path.join(os.path.dirname(__file__), 'DATA/Images/ECI Rotating Globe.png')
        Earth_texture = pv.read_texture(texture_file)
    else:  # ECEF Globe with continents
        if True:
            texture_file = os.path.join(os.path.dirname(__file__), 'DATA/Images/World with Ice (4000x2000).jpg')
            #Earth_texture = plt.imread(texture_file)  # matplotlib
            #Earth_texture = Image.open(texture_file)  # Pillow
            Earth_texture = pv.read_texture(texture_file)   # PyVista
        else:
            Earth_texture = examples.load_globe_texture()

    
    #––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
    if view_mode ==  'ECI':
        pv.Plotter(lighting = "light kit")  
        #pv.Plotter(lighting = "three lights")  # just for ambient lighting of the ECI sphere
        AMBIENT = .5    
        DIFFUSE = 0.8
        SPECULAR = 0.1
        OPACITY = 0.9
    else:
        pv.Plotter(lighting = "none")  # will use the Sun specific angle lighting 
        AMBIENT = 0.4    # 0 = pitch-black night side, 0.5 = washed out
        DIFFUSE = 0.95   # 0.95–1 strength of the Sun-lit hemisphere
        SPECULAR = 0.0   # keep low; oceans otherwise get a plastic sheen
        OPACITY = 0.5    # use < 1 to make the surface transparent

    if False:  # Writes info 
        str = f"Ambient = {AMBIENT}, Diffuse = {DIFFUSE}, Specular = {SPECULAR}"
        pl.add_text(str, position = "lower_left", font_size = 8, font_file = 'C:/Windows/Fonts/Calibri.ttf')

    pl.add_mesh(Earth_sphere, texture = Earth_texture, smooth_shading = True, 
                ambient = AMBIENT,
                diffuse = DIFFUSE,
                specular = SPECULAR,
                opacity = OPACITY, 
                specular_power = 15)  # (what does it do??)
    
    #pl.add_mesh(Earth_sphere, show_edges = True, color = 'lightblue')
    #pl.add_mesh(Earth_sphere, style = 'points', point_size = 4, color = 'blue')
    #––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

    #--- Adds arrow identificators
    d = 1.6 * RE  # distance to origin
    if view_mode ==  'ECI':
        axes_labels = [r'$\widehat{I}$', r'$\widehat{J}$', r'$\widehat{K}$']
        # p.add_text(r'$\int_0^\infty e^{-x^2} dx = \frac{\sqrt{\pi}}{2}$', position = 'upper_left', font_size = 20)
        # TODO: Add "♈︎" next to Î
    else:
        axes_labels = [r'$\widehat{X}$', r'$\widehat{Y}$', r'$\widehat{Z}$']
        # TODO: Add "λ=0°" next to X̂
        
    pl.add_point_labels(
        d * directions,
        axes_labels,
        font_size = 18,
        show_points = False,
        text_color = 'black',
        bold = False,
        #shape = None,  # (not required) 
        #fill_shape = False, # (not required)
        shape_opacity = 0.0,  # removes the shape shadows
        #always_visible = False, # (not required)
        justification_horizontal = "center",
        justification_vertical = "center",
        # Prefer background_color over shape when using centered justification
        #background_color = None, # (not required)
    )
    # always_visible = True matters for orbit plots — without it, a label on the far side of the Earth mesh gets occluded.
    # shape = None removes the default rounded-rect background;
    # keep shape = 'rounded_rect' if you need contrast against a busy starfield.

    #--- Draws Equator line
    r = RE * EXPAND_FACTOR   
    theta = np.linspace(0.0, 360, 360)
    x = r * cosd(theta)
    y = r * sind(theta)
    z = np.zeros_like(theta)  # z ≡ 0° on Equator
    lines = np.column_stack((x, y, z))
    Plot3D(pl, lines, color = 'navy', line_width = 1)

    #--- Draws Tropics
    for phi in Earth.Obliquity * np.array((-1,1)):
        #print('phi = ', phi)
        pos = Sph2Cart(theta, phi * np.ones_like(theta), r)
        pos[2,:] *= Earth.sf
        Plot_dash_3D(pl, pos, dash_length = 2, gap_length = 1, color = 'Khaki', line_width = 1.5)

    #TODO: Draw Arctic & Antarctic Circles


    #--- Creates Prime Meridian line
    if view_mode ==  'ECEF':
        phi = np.linspace(-90, 90, 360)
        x = r * cosd(phi)
        y = np.zeros_like(phi)  # y ≡ 0° on Prime Meridian
        z = r * sind(phi)
        z *= Earth.sf
        lines = np.column_stack((x, y, z))
        Plot3D(pl, lines, color = "navy", line_width = 1)  #render_lines_as_tubes = True

    #--- Identifies the Equatorial Plane
    txt = pv.Text3D('Equatorial\nplane', depth = 0, normal = (0,0,1))
    txt.scale(Earth.r1_km / 20, inplace = True)
    txt.translate(1.0 * Earth.r1_km * np.array([1., 1., 0.]), inplace = True)
    txt.rotate_z(90, point = txt.center, inplace = True)
    pl.add_mesh(txt, color = 'gray')

    
    if False: #--- Add a soft ambient light
        amb_light = pv.Light(light_type = "scene light")
        amb_light.position = (20 * Earth.r1_km, 00 * Earth.r1_km, 0 * Earth.r1_km)  # focal_point: aim it at the Earth centre
        amb_light.intensity = 0.2
        pl.add_light(amb_light)

    #pl.render() # is the part people miss. Setting the property updates the VTK object but doesn't trigger a redraw on its own.
    
    # Camera position
    pl.camera_position = [
        (8 * RE, 0, 0),   # camera location
        (0, 0, 0),        # look-at point (Earth's center)
        (0.0, 0.0, 1.0),  # up direction = +Z/K
    ]
    # default viewing geometry
    pl.camera.azimuth = 45
    pl.camera.elevation = 15
    #pl.camera.roll = ?

    pl.camera.zoom(1.0)  # useful to zoom/unzoom

    #pl.view_xy()
    #pl.view_isometric()
    #pl.disable_parallel_projection()  ? for perspective view?

    if False:  # Draws axes and grids on the back faces
        max_pos = 9000
        pl.show_bounds(
            bounds = [-max_pos, max_pos, -max_pos, max_pos, -max_pos, max_pos],
            #grid = 'back',
            location = 'outer',
            all_edges = False,  #Removes front-facing boundary lines that obscure the view.
            #all_edges = True,  
            xtitle = 'X',
            ytitle = 'Y',
            ztitle = 'Z',
            fmt = '%.0f'  
        )

    # Use terrain-style interaction: Z remains the natural up direction
    pl.enable_terrain_style(
        mouse_wheel_zooms = True,
        shift_pans = False
    )
    pl.show_axes()  # Adds the RGB orientation marker widget

    return pl

# end DrawEarth3D()


def DrawAngularMomentumVector(pl, H_):

    if False:  # (simple line)
        line = np.column_stack(([0, H_[0]], [0, H_[1]], [0, H_[2]]))
        Plot3D(pl, line, color = 'navy', line_width = 2.5)
    else:  # Nicer arrow
        arrow_geom = pv.Arrow(
            start = (0.0, 0.0, 0.0),
            direction = H_,          # (no need to normalize)
            tip_length = 0.05,
            tip_radius = 0.01,      # ← controls tip width
            shaft_radius = 0.0036,  # ← controls shaft width
            tip_resolution = 32,    # (could use 16)
            shaft_resolution = 32,
            scale = 1.5 * Earth.r1_km)
        
        pl.add_mesh(arrow_geom, color = 'navy', ambient = 0.8, diffuse = 0.2, specular = 0.0, smooth_shading = True)

    #--- Writes identifier Angular Momentum vector with "H^"
    if False:
        pts = np.array([h_pos])
        pl.add_point_labels(pts, ['h'], font_size = 14, text_color = 'navy', bold = False, justification_horizontal = "center", justification_vertical = "center",
                        show_points = False, shape_opacity = 0.0)
    else:
        label = pv.Label(r'$\widehat{H}$', position = 1.05 * H_, size = 18)
        pl.add_actor(label)
        label.prop.italic = True
        label.prop.color = 'navy' 
        label.prop.justification_horizontal = "center"
        label.prop.justification_vertical   = "center"
        #lab.position = 2 * h_pos   # move it later

# end DrawAngularMomentumVector()    


def DrawSunVector(pl, Sun, view_mode):
    """
    Earth lit by the Sun at a given right ascension and declination.

    Two light contributions:
    1. A uniform ambient term so the continents stay readable on the night side.
    2. A directional "Sun" light placed along the RA/Dec direction, which produces the terminator.

    The ambient term is a material property of the mesh, not a light. In VTK's shader the ambient contribution is added once, independently of any light in
    the scene, so it survives into the shadowed hemisphere where no light reaches. That is exactly what we want: raising AMBIENT lifts the night side without
    washing out the day side or moving the terminator.
    """

    if view_mode == 'ECI':
        lon = Sun.RA
    else:
        lon = Sun.Lon

    # Sun-Earth vector
    # see SOL Part I – 𝜌  Radius of the orbit and 𝑅_𝜙

    #Sun_vector = Sph2Cart(lon, Sun.Dec, Sun.Dist_km)  # Sun is on the Celestial Sphere with position expressed in RA/Lon & Dec

    #--- Computes point touching Earth
    h_geod = 0  # for touching the Earth surface
    R_geod, phi_geoc, R_e, R_N, R_M, x0,z0 = geod_to_pos(Sun.Dec, lon, h_geod, Earth.r1_km, Earth.r2_km)

    Sun_vect_geod = R_geod - np.stack([0,0,z0])
    s_norm = Sun_vect_geod.ravel() / norm(Sun_vect_geod)  # normalized vector (sun_hat) {1×3}

    Sun_vec_length = Earth.r1_km

    arrow_geom = pv.Arrow(
        start = R_geod + Sun_vec_length * s_norm,      # starts above Earth's surface (R_N) and extends 
        direction = -s_norm,
        tip_length = 0.1,
        tip_radius = 0.02,     # ← controls tip width
        shaft_radius = 0.005,  # ← controls shaft width
        tip_resolution = 32,     # (could use 16)
        shaft_resolution = 32,
        scale = Sun_vec_length)   # length
    # Make the arrow mostly self-lit, since it's an annotation rather than a physical object in the scene
    pl.add_mesh(arrow_geom, color = 'cadmium_lemon', ambient = 0.8, diffuse = 0.2, specular = 0.0, smooth_shading = True) 

    # Draws the complement line going to the Earth’s center
    #vec = np.outer(R_geod * s_norm, [0, 1])  # [3×2]
    vec = np.vstack((np.array([0,0,z0]), R_geod))
    Plot3D(pl, vec, color = 'orange', line_width = 2.0)

    # draws the z-axis complement 
    vec = np.vstack((np.array([[0,0,-Earth.r2_km], [0,0,0]])))
    Plot3D(pl, vec, color = 'gray', line_width = 2)
    # draws the z₀ geodetic offset 
    vec = np.vstack((np.array([[0,0,z0], [0,0,0]])))
    Plot3D(pl, vec, color = 'red', line_width = 4)

    # TODO: computes the angle between the Earth axis of rotation and the Sun ray vector; → should be equal to Sun.Dec == 90° at equinox.
    
    if False:
        # Adds a light to reflect Sun's astrometric position
        pl.remove_all_lights()  # Kills the default 
        Sun_light = pv.Light(
            position = SUN_DISTANCE * s_norm,
            focal_point = (0.0, 0.0, 0.0),
            color = "white",
            light_type = "scene light",   # world coordinates, fixed to the scene
            intensity = 0.85,  # 1 is too high
        )
        Sun_light.positional = False  # → this makes it a directional light: rays are parallel and intensity does not fall off with distance, which is correct at 1 AU.
        pl.add_light(Sun_light)

    if False:  # (writes RA/Dec at start of arrow)
        pl.add_point_labels(
            np.array([2.1 * Earth.r1_km * s_norm]),
            [f"Sun α = {Sun.RA:.1f}° δ = {Sun.Dec:.1f}°"],
            font_size = 13,
            text_color = 'orange',
            shape = None)

    #--- Draws Terminator Line (line of sunset/sunrise)
    if view_mode == 'ECEF':  # N/A in ECI

        # 1a) rotated circle
        theta = np.linspace(0.0, 360, 400)  # [°]
        C = EXPAND_FACTOR * Sph2Cart(theta, np.zeros_like(theta), Earth.r1_km)
        theta = 90 - Sun.Dec
        C = RotZ(Sun.Lon) @ RotY(theta) @ C
        Plot_dash_3D(pl, C, dash_length = 4, gap_length = 1, color = 'indigo', line_width = 2)

        # 1b) squeezing circle doesn't work for arbitrary orientation
        C[2,:] *= Earth.sf
        Plot_dash_3D(pl, C, dash_length = 4, gap_length = 1, color = 'red', line_width = 2)

        # 2) GreatCircle
        lat_GC, lon_GC, xyz_GC = GreatCircle(lon, Sun.Dec)
        xyz_GC *= Earth.r1_km
        Plot3D(pl, xyz_GC, color="orange", line_width=3)  

        # 3) Rotated "Great Ellipse" (Great Circle on an Ellipsoid)
        obl_geod = Earth.Obliquity  # geodetic property or the Earth axis of rotation wrt to ecliptic plane around the Sun
        obl_geoc_ = atand(tand(obl_geod) * (Earth.r2_km/Earth.r1_km)**2)  # φ′ [°]
        print('obl_geod = ', obl_geod)
        print('obl_geoc_ = ', obl_geoc_)

        lat_GE, lon_GE, xyz_GE, lat_geoc = GreatEllipsoid(lon, Sun.Dec, Earth.r1_km, Earth.r2_km) 
        Plot3D(pl, xyz_GE, color="black", line_width=3)  

    #return Sun_vector
        
# end DrawSunVector()


def DrawEclipticLine(pl, GMST, Earth):
# Draws 3D Ecliptic plane line of Earth's surface

    Obliquity = Earth.Obliquity  # = 23.4358° (2026)
    obl_geoc_ = atand(tand(Obliquity) * (Earth.r2_km/Earth.r1_km)**2);  # φ′
    # print('Obliquity ε = %.4f°', Obliquity)
    # print('φ′ = %.4f°', obl_geoc_)       

    # Draws the Obliquity vector
    V = 1.2 * Earth.r1_km * RotZ(-GMST) @ RotX(Obliquity) @ np.array([[0,0,-1], [0,0,+1]]).T  # 3×2
    Plot3D(pl, V.T, color = 'blue', line_width = 1.5)
    #Plot_dash_3D(pl, C, dash_length = 4, gap_length = 1, color = 'blue', line_width = 2)

    lat, lon, S, lat_geoc = GreatEllipsoid(90 - GMST, obl_geoc_ - 90, Earth.r1_km, Earth.r2_km) 

    if False:
        Plot3D(pl, S, color = 'cadmium_lemon', line_width = 2)
    else:
        Plot_dash_3D(pl, S, dash_length = 3, gap_length = 1, color = 'cadmium_lemon', line_width = 3)

"""
==> Not really applicable to the ecliptic line, which is tied to the Earth
    # → Draws the vector normal to the ecliptic, i.e. the vector perpendicular to the Sun's ecliptic path on the Earth surface
    E = np.cross(S[:,0], S[:,1])  # 1×3 vector
    N = np.column_stack((np.zeros(3), E.T))  # Creates the vector with origin at Earth's center
    N_ = 1.5 * Earth.r1_km * N / norm(N)  # scales for display
    N_[2,:] *= Earth.sf
    Plot3D(pl, N_, color = 'cadmium_lemon', line_width = 2)

    # Note: will be similar to the Angular Momentum vector Ĥ, while this one can be computed directly by Ĥ = R̂ × V̂

    # → Computes the angle between the Sun vector and this angle
    V_ = V[:,1] - V[:,0]

    gamma = Angle(V_, S[:,0])  # ≡ acosd(np.dot(V_, S[:,0]) / (norm(V_) * norm(S[:,0])))
    print(f'  γ = {gamma:.4f}°')  # ≡ 90° by construction
    # → only true from first point at λ(0), otherwise must take into account 
"""

# end of DrawEclipticLine



#TODO: Identifies Perigee & Apogee points on the orbit

#TODO: Draws the line of apses with Ascending / Descending node symbols ☊ ☋
#TODO: Calculate & draw β☉ angle


def Plot3D(pl, xyz, **kwargs):
    """
    Plots a 3D line from an array of points using PyVista.
    
    Parameters:
    - xyz (array_like): An (N, 3) array of points defining the line.
    - plotter (pv.Plotter): The active PyVista plotter instance to add the mesh to.
    - **kwargs: Additional arguments passed directly to plotter.add_mesh.
    """

    m, n = xyz.shape
    if n != 3 and m == 3:
        xyz = xyz.T  # PyVista needs the transpose (.T), since it expects n×3:
    elif n != 3:
        raise TypeError('>>> Plot3D: argument must be a 3×N or N×3 array')
    
    # Create the PyVista polydata line mesh from the points
    mesh = pv.lines_from_points(xyz)  #TODO: find a way to pass the info closed or not
    
    # Add the mesh to the provided plotter instance
    pl.add_mesh(mesh, **kwargs)

# end Plot3D


def Plot_dash_3D(pl, points, dash_length = 3, gap_length = 1, **kwargs):
    """Generates a geometrically dashed 3D line segment array inside PyVista."""
    # PyVista does not have a native line style property for 3D geometry, but you can draw dashed lines by manually building a mesh of alternating short line segments
    # Because PyVista is built on VTK (which natively handles geometric elements rather than pixel-based line styling), PyVista does not support native 3D dashed line styles.

    m, n = points.shape
    if n != 3 and m == 3:
        points = points.T  # PyVista needs the transpose (.T), since it expects n×3:
    elif n != 3:
        raise TypeError('>>> Plot3D: argument must be a 3×N or N×3 array')
    
    num_points = len(points)

    # Track position sequence
    i = 0
    while i < num_points - 1:
        # Determine endpoints of the current visible dash segment
        end_idx = min(i + dash_length, num_points)
        dash_points = points[i:end_idx]
        
        # Build individual PolyData line segments for PyVista
        if len(dash_points) > 1:
            # Format lines array as required by PyVista PolyData: [num_pts, pt0, pt1, ...]
            lines = np.hstack(([len(dash_points)], np.arange(len(dash_points))))
            dash_mesh = pv.PolyData(dash_points, lines = lines)
            pl.add_mesh(dash_mesh, **kwargs)
            
        # Jump ahead past the dash and the invisible gap space
        i +=  dash_length + gap_length

# end Plot_dash_3D


def view(pl, az, el):
    pl.camera.azimuth = az  # [°]
    pl.camera.elevation = el   # [°]
# end view()