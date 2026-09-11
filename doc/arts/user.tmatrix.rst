T-matrix interface
==================

The Python interface is ``pyarts3.arts.tmatrix``.  ``available()`` reports
whether the backend was built; ``extended_precision()`` identifies the selected
variant.  ``fixed`` computes a particle and evaluates its amplitude and phase
matrices at one illumination/scattering geometry.  ``fixed_batch`` computes
one particle and evaluates a matrix of geometries, returning a list of fixed
results in row order.  Each row contains incident zenith, scattered zenith,
incident azimuth, scattered azimuth, alpha and beta, all in degrees.
``random`` computes a
randomly oriented size distribution, defaulting to an effectively monodisperse
particle.  These functions calculate individual optical results.  Use
``ParticleHabit.tmatrix`` below to prepare a habit for scattering species.

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

Using particles as scattering species
-------------------------------------

``ParticleHabit.tmatrix`` generates a totally randomly oriented habit in memory.
Supply temperature and frequency grids, volume-equivalent diameters in metres,
a complex refractive-index matrix with axes (temperature, frequency), and
material density in kg/m³.  The refractive index uses a nonnegative imaginary
part.  Shape and aspect-ratio conventions match the direct T-matrix interface.

For example, a single-size ice population can be set up as follows::

    import numpy as np
    import pyarts3 as pa

    A = pa.arts
    density = 917.0
    habit = A.ParticleHabit.tmatrix(
        t_grid=[250.0, 280.0],
        f_grid=[229e9, 231e9],
        diameters=[200e-6],
        refractive_index=np.full((2, 2), 1.78 + 0.01j),
        density=density,
        aspect_ratio=1.5,
        angles=181,
    )
    number = A.ScatteringSpeciesProperty(
        "ice", A.ParticulateProperty.NumberDensity)
    psd = A.MonodispersePSD(number, 250.0, 280.0)
    ws = pa.Workspace()
    ws.scat_species = [A.ScatteringHabit(
        habit, psd, density * np.pi / 6, 3.0)]

The constant refractive index above is illustrative; supply the material model
appropriate to the calculation.  Set the ``number`` property in the atmosphere
to the particle number density in m⁻³.  The explicit mass-size relation uses
volume-equivalent diameter: mass equals ``density * pi / 6 * diameter**3``.
``MonodispersePSD`` requires exactly one diameter; habits containing several
sizes can be combined with the other ARTS PSDs.

The factory computes each particle size separately.  ``ScatteringHabit`` then
applies the atmospheric PSD and interpolates the stored optical properties;
it does not rerun T-matrix during bulk-property evaluation.  Choose the
sampling grids to resolve changes in the optical properties.  ``angles`` sets
an equally spaced scattering-angle grid from 0 to 180 degrees.

The generated data already have physical cross-section units, including
forward and backscatter matrices.  No additional phase normalization is needed
before passing the habit to ``ScatteringHabit``.  The factory describes totally
random orientations; it does not generate aligned or azimuthally random
particle populations.
