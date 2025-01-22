from dustmaps.bayestar import BayestarQuery
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
from dustmaps.config import config

# Set the local data storage path for dustmaps
config['data_dir'] = '/Users/baotong/dustmaps/'  
# Load the Bayestar extinction model (local data)
bayestar = BayestarQuery(version='bayestar2019')
coords = SkyCoord(ra=266.41683, dec=-29.00781, unit=(u.deg, u.deg), frame='icrs')
# Query extinction distribution in "samples" mode
ebv_samples = bayestar(coords, mode='samples')  
distances = bayestar.distances  # Retrieve distance nodes (in kpc)
print("Distance (kpc) vs E(B-V):")
for d, e in zip(distances, ebv_samples[0]):
    print(f"Distance: {d:.2f} kpc, E(B-V): {e:.4f}")
# Interpolate to find E(B-V) at a specific target distance
target_distance = 2.2  # Target distance in kpc
distances_pc = np.array(distances)  # Convert distance nodes to array (already in kpc)
ebv_target = np.interp(target_distance, distances_pc, ebv_samples[0])
print(f"\nTarget Distance: {target_distance} kpc, E(B-V): {ebv_target:.4f}")

R_G, R_BP, R_RP = 2.74, 3.14, 2.05  # Extinction coefficients
A_G = R_G * ebv_target
A_BP = R_BP * ebv_target
A_RP = R_RP * ebv_target
print(f"A_G = {A_G:.2f}, A_BP = {A_BP:.2f}, A_RP = {A_RP:.2f}")

m_observed_G = 15.0  # Example observed G-band magnitude
m_corrected_G = m_observed_G - A_G
print(f"Corrected G magnitude: {m_corrected_G:.2f}")
