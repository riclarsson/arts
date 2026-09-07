T-matrix interface
==================

The Python interface is ``pyarts3.arts.tmatrix``.  ``available()`` reports
whether the backend was built; ``extended_precision()`` identifies the selected
variant.  ``fixed`` computes a particle and evaluates its amplitude and phase
matrices at one illumination/scattering geometry.  ``random`` computes a
randomly oriented size distribution, defaulting to an effectively monodisperse
particle.  These are particle calculations, not workspace radiative-transfer
methods or an automatic scattering-habit generator.

Units and polarization
----------------------

Angles are in degrees.  For fixed particles and the default and power-law
size distributions, use a common length unit for radius and wavelength.
Metres give SI outputs.  ``radius_ratio=1`` selects volume-equivalent radius;
other positive values select surface-area-equivalent radius for these shapes.
``shape=-1`` means spheroid, with horizontal/rotational axis ratio;
``shape=-2`` means cylinder, with diameter/length ratio.

The fixed result contains a complex 2 by 2 Jones amplitude matrix (length),
a ``Muelmat`` phase matrix (length squared), and orientation-averaged
extinction and scattering cross sections (length squared).  The Jones matrix
is not a ``Specmat``, which represents a complex 4 by 4 Mueller matrix.
The amplitude-to-Mueller conversion preserves the ARTS2 Stokes convention.

Random results contain ``MuelmatVector`` phase matrices in the scattering-plane
basis at equally spaced scattering angles from 0 to 180 degrees.  They are
dimensionless, with the first element integrating to :math:`4\pi` over solid
angle.  Multiply by ``scattering / (4*pi)`` to obtain differential scattering
cross sections.  The nonzero elements are F11, F22, F33, F44, symmetric F12,
and antisymmetric F34, following the original solver convention.

Results remain valid after subsequent calls.

See :doc:`concept.tmatrix` for the physical definitions and normalization.
Build configuration and implementation notes are in :doc:`dev.tmatrix`.
