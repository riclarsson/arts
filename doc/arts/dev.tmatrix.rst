T-matrix implementation and validation
======================================

The optional Mishchenko T-matrix backend computes scattering by homogeneous
spheroids and finite circular cylinders.  Enable it with
``-DENABLE_TMATRIX=ON``.  A Fortran compiler is required only when enabled.
``-DENABLE_TMATRIX_QUAD=ON`` selects the retained extended-precision solver;
inputs and outputs remain double precision, and both solvers retain some
single-precision internal storage.  The two backends are alternative builds.

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

See :doc:`user.tmatrix` for the interface and :doc:`concept.tmatrix` for
physical definitions and normalization.

Native particle-habit integration
--------------------------------

``ParticleHabit::tmatrix`` in ``particle_habit_tmatrix.cc`` generates native
TRO gridded single-scattering data, without legacy-data adapters or files.
The scattering library links to the optional-backend adapter; when the backend
is disabled the factory reports its unavailability through the same exception
path as the direct interface.

For each diameter, temperature, and frequency, the factory uses the direct
random-orientation interface with one size quadrature point at the specified
volume-equivalent radius.  It does not perform the atmospheric PSD integration.
The Mueller matrices are multiplied by ``Csca/(4*pi)`` and packed in native
TRO order F11, F12, F22, F33, F34, F44.  Extinction is Cext and absorption is
Cext minus Csca.  Both endpoint scattering matrices are extracted from the
populated phase data.

Particle mass is density times the equivalent-sphere volume.  For spheroids,
maximum dimension is the larger axial diameter.  For cylinders it is the
largest point-to-point distance, including both cylinder length and diameter.
Temperature/frequency grids and all material inputs are validated before
starting the solver calls.  The existing T-matrix mutex serializes Fortran
access, and bulk evaluation subsequently uses the ordinary habit machinery.

``tests/core/tmatrix/habit.py`` checks the direct-to-native normalization and
component ordering, forward/backscatter extraction, metadata, and the analytic
small-sphere Rayleigh limit.  It also exercises temperature/frequency
interpolation and repeated number-density scaling through
``ArrayOfScatteringSpecies``, verifying that these evaluations leave the
original particle data intact.
