#!/bin/bash

# Leap second kernel:
wget -nc https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls

# Generic barycenter ephemeris:
wget -nc https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp

# Planetary satellites:
wget -nc https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/mar097.bsp
wget -nc https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/jup344.bsp
wget -nc https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/jup365.bsp
wget -nc https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/sat393.bsp
wget -nc https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/ura111.bsp
wget -nc https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/nep086.bsp
wget -nc https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/plu058.bsp

# Planetary gravity constants:
wget -nc https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de431.tpc

# Planetary radii:
wget -nc https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc