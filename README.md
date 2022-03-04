# CelNav Catalog
This is a project meant for generating a complete catalog of objects from JPL's HORIZONS system, and identifying potentially visible objects which can then be given to an onboard system for celestial navigation.

## Generate Catalog with Python
1. **Download the relevant SPICE kernels:** `cd generic_kernels/; ./get_kernels.sh`

## Building C++ Integrator
```
mkdir build
cd build
cmake ..
make
```
## Run the C++ Integrator
- `./celnav <catalog.csv> <ephemeris-time>`

where `<catalog.csv>` is the path to a catalog .csv file created by the python `Catalog` class, and `<ephermis-time>` is the ephemeris time to evaluate the catalog objects at.