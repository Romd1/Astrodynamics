import numpy as np
from SOL_Tools.Math_tools import * 
from SOL_Tools.Orbit_tools import * 

from tracking import *

R_geoc = geod_to_pos(45,0 , 1000)

pos = np.asarray(R_geoc[0], dtype=float).reshape(3) # ECEF position [km]
Lat_geoc = R_geoc[1]            # Geocentric latitude [deg]
R_e = R_geoc[2]                 # Geocentric radius [km]
z0 = R_geoc[6]                  # Auxiliary z-coordinate [km]

print(ECEF_to_AER(np.vstack((0,0,-1)), 0, 0))