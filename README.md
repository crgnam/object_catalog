# CelNav Catalog
This is a project meant for generating a complete catalog of objects from JPL's HORIZONS system, and identifying potentially visible objects which can then be given to an onboard system for celestial navigation.

## Generate Catalog with Python
*Coming Soon*

## Building C++ Integrator
```
mkdir build
cd build
cmake ..
make
```
## Run the C++ Integrator
- `./propagate <catalog.csv> <ephemeris-time>`

where `<catalog.csv>` is the path to a catalog .csv file created by python, and `<ephermis-time>` is the ephemeris time to evaluate the catalog objects at.