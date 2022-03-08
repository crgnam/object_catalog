import spiceypy as spice
import numpy as np
from datetime import datetime

from celnav_py import Catalog

# Load SPICE kernels:
spice.furnsh("meta_kernel.tm")

# Define the observer position and time:
observer_position = np.array([3e8, 0, 0])
now = datetime.now() # current date and time
et = spice.str2et(now.strftime("%Y-%m-%d %H:%M:%S"))
magnitude_maximum = 10
sun_angle_minimum = np.deg2rad(15)

# Create/load the catalog:
my_catalog = Catalog('catalog.pickle')

# Obtain the states of all objects in the catalog:
asteroid_positions = my_catalog.get_asteroid_states(et)

# Obtain the position of all planets:
planet_positions, spk_ids, planet_H = my_catalog.get_planet_states(et)

# Determine the potentially visible objects:
visible_asteroids = my_catalog.get_visible(observer_position, asteroid_positions, my_catalog.H, my_catalog.G, magnitude_maximum, sun_angle_minimum)
visible_planets   = my_catalog.get_visible(observer_position, planet_positions, planet_H, np.nan, magnitude_maximum, sun_angle_minimum)

# Calculate the orbital elements of the planets/moons:
planet_spk = spk_ids[visible_planets]
_,mu_sun = spice.bodvrd( 'SUN', 'GM', 1 )
planet_elements = np.zeros((planet_spk.size, 1+1+6+3))
for idx, spk_id in enumerate(planet_spk):

    # Obtain the gravity constants if a planet:
    if np.mod(spk_id,100) == 99:
        _,radii = spice.bodvrd( str(spk_id), 'RADII', 3 )
        R = np.mean(radii)
        _,J2 = spice.bodvrd( str(spk_id), 'J2', 1 )
        _,mu = spice.bodvrd( str(spk_id), 'GM', 1 )
        mu = mu[0]
        J2 = J2[0]
        parent_object = 'SUN'
    else:
        mu = np.nan
        R = np.nan
        J2 = np.nan
        parent_object = str(spk_id - np.mod(spk_id,100) + 99)

    # Calculate the object's orbital elements:
    state,_ = spice.spkezr(str(spk_id), et, 'J2000','NONE', parent_object)
    if parent_object == 'SUN':
        spice_elems = spice.oscelt(state,et,mu_sun)
    else:
        _,mu_parent = spice.bodvrd( parent_object, 'GM', 1 )
        spice_elems = spice.oscelt(state,et,mu_parent)

    # Convert spice elems to kepler elements:
    rp = spice_elems[0]
    e = spice_elems[1]
    i = spice_elems[2]
    node = spice_elems[3]
    peri = spice_elems[4]
    M = spice_elems[5]
    a = rp/(1-e)

    # Store the data into a matrix:
    planet_elements[idx,:] = np.array([spk_id, et, a,e,i,peri,node,M, mu,R,J2])

# Format elements for asteroids:
asteroid_elements = np.nan*np.zeros((np.sum(visible_asteroids), 1+1+6+3))
asteroid_elements[:,0] = my_catalog.spkid[visible_asteroids]
asteroid_elements[:,1] = my_catalog.epoch[visible_asteroids]
asteroid_elements[:,2] = my_catalog.a[visible_asteroids]
asteroid_elements[:,3] = my_catalog.e[visible_asteroids]
asteroid_elements[:,4] = my_catalog.i[visible_asteroids]
asteroid_elements[:,5] = my_catalog.peri[visible_asteroids]
asteroid_elements[:,6] = my_catalog.node[visible_asteroids]
asteroid_elements[:,7] = my_catalog.M[visible_asteroids]

# Write the Catalog .csv:
output_elements = np.vstack((planet_elements, asteroid_elements))
np.savetxt("catalog.csv", output_elements, delimiter=",")

# Write results to a csv (for plotting in MATLAB):
asteroid_data = np.hstack((asteroid_positions.T, visible_asteroids.reshape((-1,1)) ))
planet_data   = np.hstack((planet_positions.T, visible_planets.reshape((-1,1)) ))
output_data   = np.vstack((planet_data, asteroid_data))

np.savetxt("test_data.csv", output_data, delimiter=",")
np.savetxt("planet_positions.csv", planet_positions.T)