ARTS2 radar droplet reconstruction
=================================

``radar_single_scattering.py`` generates the former ``droplet_50um.xml``
particle using ``ParticleHabit.sphere`` and the native ARTS Mie solver.
The original diameter (50 micrometres), frequency (94 GHz), temperatures
(268.15, 273.15, 278.15 K) and 181 scattering angles are preserved.

Water permittivity follows Atmlab ``eps_water_liebe93.m``, including its
relaxation-frequency coefficient 146 (not 146.4 in ``mie/epswater93.m``).
The formula is included directly in the test; no Atmlab installation is needed.
Compared with the original data using identical refractive indices, relative
errors were below 6.7e-13 for extinction/absorption and 1.7e-12 for the phase
matrix normalized by its peak. The small F34 component agreed within 1.2e-10
relative to its own peak.

The habit now has physical sphere mass, with density 1000 kg/m3, and explicit
forward/backscatter matrices. Particle abundance remains number density.
Original radar reference values and acceptance tolerances are unchanged.
