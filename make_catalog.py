import spiceypy as spice

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

from celnav_py import Catalog

def set_axes_equal(ax):
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

spice.furnsh("meta_kernel.tm")

my_catalog = Catalog()

now = datetime.now() # current date and time
epoch = spice.str2et(now.strftime("%Y-%m-%d %H:%M:%S"))
et_array = np.arange(epoch, epoch + 5*365*86400, int(10*86400))
r = my_catalog.get_states(et = et_array)
print(r)

# Plot the results:
ax = plt.axes(projection='3d')

# Data for a three-dimensional line
zline = np.linspace(0, 15, 1000)
xline = np.sin(zline)
yline = np.cos(zline)
ax.plot3D(r[0,:], r[1,:], r[2,:], 'gray')
set_axes_equal(ax)
plt.show()