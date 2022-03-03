import spiceypy as spice
import numpy as np
from datetime import datetime

from celnav_py import Catalog

# Load SPICE kernels:
spice.furnsh("meta_kernel.tm")

# Define the observer position and time:
observer_position = np.array([3e8, 0, 0])
now = datetime.now() # current date and time
epoch = spice.str2et(now.strftime("%Y-%m-%d %H:%M:%S"))
magnitude_maximum = 12
sun_angle_minimum = np.deg2rad(30)

# Create/load the catalog:
my_catalog = Catalog('catalog.pickle')

# Obtain the states of all objects in the catalog:
object_positions = my_catalog.get_states(et = epoch)

# Determine the potentially visible objects:
visible_objects = my_catalog.get_visible(observer_position, object_positions, my_catalog.H, my_catalog.G, magnitude_maximum, sun_angle_minimum)

# Write results to a csv (for plotting in MATLAB):
output_data = np.hstack((object_positions.T, visible_objects.reshape((-1,1)) ))
np.savetxt("test_data.csv", output_data, delimiter=",")