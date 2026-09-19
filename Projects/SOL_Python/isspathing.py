import numpy as np
from SOL_Tools.Tracking_Obj import * # Definit les classes et methodes propres à un observateur et une cible


R_ECEF_sat = np.array([1457.6991, -4505.404, 4858.8987]) # Position de l'Iss à l'époque ou elle est la plus proche

Obs = Observer(46.781,-71.278,0.090) # Observateur
Iss = Target(R_ECEF_sat) # Iss à la position de tracking

Azim,Elev,Range,SEZ = Obs.aer(Iss) # Get Azimut, Elevation, Range, and the SEZ position
Sub_Lat_geod, Sub_Lon, Sub_ECEF = Iss.subsatellite_point() # Get the sub-satellite point in Geod angles and ECEF coord

# Formatage
Sub_Lat_geod = Sub_Lat_geod.item()
Sub_Lon = Sub_Lon.item()
Sub_ECEF = np.asarray(Sub_ECEF).flatten()

# Affichage
print("\n=== AER ===")
print(f"Azimuth             : {Azim}°")
print(f"Elevation           : {Elev}°")
print(f"Range               : {Range} km")
print(f"SEZ                 : {SEZ} km")

print("=== Subsatellite point ===")
print(f"Latitude géodésique : {Sub_Lat_geod:.6f}°")
print(f"Longitude           : {Sub_Lon:.6f}°")
print(f"ECEF                : [{Sub_ECEF[0]:.6f}, "
      f"{Sub_ECEF[1]:.6f}, {Sub_ECEF[2]:.6f}] km")