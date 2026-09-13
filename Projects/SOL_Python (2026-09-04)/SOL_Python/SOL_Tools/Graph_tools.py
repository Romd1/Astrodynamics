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


def DrawEarth3D(pl, view_mode):
# Textured 3D Earth using PyVista

    # Earth equatorial radius
    RE = Earth.r1_km  # [km]

    #--- Arrow locations and directions

    directions = np.eye(3)

    for d in directions:
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
        SPECULAR = 0.0  # keep low; oceans otherwise get a plastic sheen
        OPACITY = 1

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
    lines = pv.lines_from_points(np.column_stack((x, y, z)), close = True)
    pl.add_mesh(lines, color = "navy", line_width = 1)

    #--- Draws Tropics
    for phi in Earth.Obliquity * np.array((-1,1)):
        #print('phi = ', phi)
        r_ = r * cosd(phi)
        x = r_ * cosd(theta)
        y = r_ * sind(theta)
        z = r * sind(phi) + np.zeros_like(theta)  # z ≡ ±𝜖° on Tropics of Cancer and Capricorn
        #lines = pv.lines_from_points(np.column_stack((x, y, z)), close = True)
        #pl.add_mesh(lines, color = "Khaki", line_width = 0.5)  
        add_3d_dashed_line(pl, x, y, z, dash_length = 2, gap_length = 1, color = 'Khaki', line_width = 1.5)

    #TODO: Draw Arctic & Antarctic Circles


    #--- Creates Prime Meridian line
    if view_mode ==  'ECEF':
        phi = np.linspace(-90, 90, 360)
        x = r * cosd(phi)
        y = np.zeros_like(phi)  # y ≡ 0° on Prime Meridian
        z = r * sind(phi)
        lines = pv.lines_from_points(np.column_stack((x, y, z)), close = False)
        pl.add_mesh(lines, color = "navy", line_width = 1)  #render_lines_as_tubes = True

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
        line = pv.lines_from_points(np.column_stack(([0, H_[0]], [0, H_[1]], [0, H_[2]])))
        pl.add_mesh(line, color = 'navy', line_width = 2.5)
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


def DrawSunVector(pl, Sun):
    """
    Earth lit by the Sun at a given right ascension and declination.

    Two light contributions:
    1. A uniform ambient term so the continents stay readable on the night side.
    2. A directional "Sun" light placed along the RA/Dec direction, which produces the terminator.

    The ambient term is a material property of the mesh, not a light. In VTK's shader the ambient contribution is added once, independently of any light in
    the scene, so it survives into the shadowed hemisphere where no light reaches. That is exactly what we want: raising AMBIENT lifts the night side without
    washing out the day side or moving the terminator.
    """

    SUN_DISTANCE = 1.0e9  # [km] Only the direction matters for a directional light, but a large value keeps the geometry unambiguous.

    Sun_Earth_ECI = Sph2Cart(Sun.Lon, Sun.Dec, Earth.r1_km, True).ravel()  # point on Earth Surface
    s_norm = Sun_Earth_ECI / norm(Sun_Earth_ECI)  # sun_hat

    if False:
        line = pv.lines_from_points(np.column_stack(    Earth.r1_km * np.outer(s_norm, [1, 2])     ))
        pl.add_mesh(line, color = 'orange', line_width = 2.5)
    else:
        arrow_geom = pv.Arrow(
            start = 2 * Earth.r1_km * s_norm,  # starts at Earth's surface
            direction = -s_norm,
            tip_length = 0.1,
            tip_radius = 0.02,     # ← controls tip width
            shaft_radius = 0.005,  # ← controls shaft width
            tip_resolution = 32,     # (could use 16)
            shaft_resolution = 32,
            scale = Earth.r1_km)   # to Earth's surface
        # Make the arrow mostly self-lit, since it's an annotation rather than a physical object in the scene
        pl.add_mesh(arrow_geom, color = 'cadmium_lemon', ambient = 0.8, diffuse = 0.2, specular = 0.0, smooth_shading = True) 

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

    #--- Draws Terminator Line
    theta = np.linspace(0.0, 360, 100)  # [°]
    C = EXPAND_FACTOR * Sph2Cart(theta, np.zeros_like(theta), Earth.r1_km)
    C = RotZ(Sun.Lon) @ RotY(90 - Sun.Dec) @ C
    pl.add_mesh(pv.lines_from_points(C.T), color = 'dim_gray', line_width = 2.0)
        
# end DrawSunVector()


def DrawEclipticLine(pl, GMST, Earth):
# Draws 3D Ecliptic plane line of Earth's surface

    # Creates a ring 
    theta = np.linspace(0.0, 360, 101)
    # Note: Matlab linspace(0, 360, 101) is inclusive at both ends, so the first and last points coincide — 101 points describing 100 segments. 
    # That's what you want for a closed ring you'll draw as a polyline, and it matches MATLAB exactly. 
    # If you ever want 101 distinct points instead, use np.linspace(0, 360, 101)[:-1] or np.arange(0, 360, 3.6).
    ring = np.array([cosd(theta), sind(theta), np.zeros_like(theta)])  # along the equator Z = 0°

    # Rotates the ring about X axis
    Obliquity = Earth.Obliquity  # = 23.436° (2024)
    S = RotX(Obliquity) @ ring  # positive, as inclination  [3x101]

    # Rotates disk about Z axis: over the primary direction or the First Point of Aries (♈︎)
    S = RotZ(-GMST) @ S  # ← negative here, because from ECI to ECEF

    # ECEF because GMST if from I to X, and ECI is the frame of reference of the current graph

    S_ = EXPAND_FACTOR * Earth.r1_km * S  # scales for graphic
    if True:
        # PyVista needs the transpose (.T), since it expects n×3:
        pl.add_mesh(pv.lines_from_points(S_.T), color = 'cadmium_lemon', line_width = 2)
    else:
        x = S[0,:]
        y = S[1,:]
        z = S[2,:]
        add_3d_dashed_line(pl, x, y, z, dash_length = 2, gap_length = 1, color = 'cadmium_lemon', line_width = 1)

    # → Draws the vector normal to the ecliptic, i.e. the vector perpendicular to the Sun's ecliptic path on the Earth surface
    E = np.cross(S[:,1], S[:,2])  # 1×3 vector
    N = np.column_stack((np.zeros(3), E.T))  # Creates the vector with origin at Earth's center
    N_ = 1.5 * Earth.r1_km * N / norm(N)  # scales for display
    pl.add_mesh(pv.lines_from_points(N_.T), color = 'cadmium_lemon', line_width = 2)
   
    # Note: will be similar to the Angular Momentum vector Ĥ, while this one can be computed directly by Ĥ = R̂ × V̂

    # → Computes the angle between the Sun vector and this angle
    if False:
        gamma = np.rad2deg(np.arccos(np.dot(E, S[:,1]) / (np.linalg.norm(E) * np.linalg.norm(S[:,1]))));
    else:
        gamma = Angle(E, S[:,1])

    print(f'  γ = {gamma:.4f}°')  # ≡ 90° by construction

# end of DrawEclipticLine


#TODO: Identifies Perigee & Apogee points on the orbit

#TODO: Draws the line of apses with Ascending / Descending node symbols ☊ ☋
#TODO: Calculate & draw β☉ angle

def add_3d_dashed_line(plotter, x, y, z, dash_length = 5, gap_length = 3, **kwargs):
    """Generates a geometrically dashed 3D line segment array inside PyVista."""
    # PyVista does not have a native line style property for 3D geometry, but you can draw dashed lines by manually building a mesh of alternating short line segments
    # Because PyVista is built on VTK (which natively handles geometric elements rather than pixel-based line styling), PyVista does not support native 3D dashed line styles.

    points = np.column_stack((x, y, z))
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
            plotter.add_mesh(dash_mesh, **kwargs)
            
        # Jump ahead past the dash and the invisible gap space
        i +=  dash_length + gap_length