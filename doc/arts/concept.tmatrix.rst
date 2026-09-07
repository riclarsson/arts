T-matrix particle scattering
============================

The optional Mishchenko T-matrix backend computes scattering by homogeneous
spheroids and finite circular cylinders.  Enable it with
``-DENABLE_TMATRIX=ON``.  A Fortran compiler is required only when enabled.
``-DENABLE_TMATRIX_QUAD=ON`` selects the retained extended-precision solver;
inputs and outputs remain double precision, and both solvers retain some
single-precision internal storage.  The two backends are alternative builds.

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

Each call owns its output and holds one common mutex throughout computation
and amplitude evaluation.  Results survive subsequent calls.  Concurrent
calls are safe but do not run the Fortran solver concurrently.

Port and reference tests
------------------------

The complete ARTS2 ``3rdparty/tmatrix`` directory was retained, including its
original programs, parameter files, README, license and four reference files.
Array limits and numerical algorithms were preserved.  The build file was
adapted for optional Fortran compilation, position-independent code, compiler
symbol mangling, and GNU equivalents of the old quad-precision intrinsics.
The reference programs remain optional targets ``tmatrix_ampld`` and
``tmatrix_tmd``.

The following defects in the ARTS-specific sources were corrected:

* Initialize the Fortran error string on entry, and consistently use default
  Fortran integers for the quiet flag in callers and callees.
* Return immediately when the random-orientation stored expansion exceeds
  NPN4; the old wrapper reported failure but continued writing the arrays.
* Return distribution-averaged CEXT and CSCA, matching the returned albedo and
  phase matrix, rather than cross sections from the last size quadrature point.
* In the quad fixed wrapper, honor the integer shape argument, call its own
  ACJB routine, and stop redirecting standard output to a file named ``test``.
* Declare all six output arrays in the quad random wrapper.  Complete the
  ARTS LAPACK symbol prefix on its ``tmzswap`` implementation.

``tests/core/tmatrix/reference.py`` exercises the selected backend through
nanobind and the C++ interface, using the unchanged ``.ref`` files.  It compares
complex amplitudes, Mueller elements, size-distribution cross sections,
albedo, asymmetry and effective size statistics.  Printed precision determines
the comparison tolerances.  Timing and printed expansion-coefficient tables
are not part of the exposed result or these comparisons.

The fixed double reference uses volume-equivalent radius; the fixed quad
reference uses surface-area-equivalent radius.  The random double reference
uses refractive index 1.53+0.008i; the quad reference uses 1.33+0.001i.
These differences are inputs, not relaxed numerical comparisons.

The random reference's two distributions are evaluated separately, preserving
its seven and four size quadrature points.  The final extra CEXT/CSCA line in
``tmatrix_tmd.ref`` records the old return-value bug; the test compares the
unchanged distribution-average summaries printed earlier in that same file.
No reference numbers were replaced.  Radius normalization at the C++ boundary
for the power-law and gamma distributions avoids the original POWER routine's
absolute root-bracket scale when using SI-sized particles.

Tests also cover unit scaling, invalid inputs, convergence failure, ownership,
and interleaved Python-thread calls.  The original backend retains its
convergence envelope and some internal fatal error paths; the port is not a
claim that arbitrary particle parameters converge.  Input angle validation
prevents the directly exposed AMPL angle-error STOP paths.
