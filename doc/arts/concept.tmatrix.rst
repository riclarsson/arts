T-matrix particle scattering
============================

The T-matrix method relates an incident electromagnetic field to the field
scattered by a particle.  Both fields are expanded in vector spherical waves;
for expansion-coefficient vectors :math:`\boldsymbol a` and
:math:`\boldsymbol b`, the linear relation is

.. math::

   \boldsymbol b = \boldsymbol T\boldsymbol a.

For a given particle and wavelength, the same T-matrix can be used for different
illumination and observation directions.  The particle model considered here
is a homogeneous spheroid or finite circular cylinder.

Particle size and orientation
-----------------------------

Volume-equivalent radius is the radius of a sphere with the same volume as the
particle.  Surface-area-equivalent radius is the radius of a sphere with the
same surface area.  These definitions agree for a sphere, but generally differ
for nonspherical particles.

A fixed-orientation calculation specifies the particle orientation as well as
the incident and scattered directions.  A random-orientation calculation
averages over particle orientations.  Averaging over a particle size
distribution is an additional operation.

Amplitude and phase matrices
----------------------------

The complex 2 by 2 Jones amplitude matrix acts on electric-field polarization
components.  The corresponding real 4 by 4 Mueller phase matrix acts on the
Stokes vector.  Its elements are quadratic combinations of the complex
amplitudes.  With an amplitude measured in length units, the phase matrix has
units of length squared.

For randomly oriented particles of the shapes considered here, the normalized
phase matrix in the scattering-plane basis has the form

.. math::

   \boldsymbol F(\Theta) =
   \begin{pmatrix}
   F_{11} & F_{12} & 0 & 0 \\
   F_{12} & F_{22} & 0 & 0 \\
   0 & 0 & F_{33} & F_{34} \\
   0 & 0 & -F_{34} & F_{44}
   \end{pmatrix},\qquad
   \int_{4\pi} F_{11}\,d\Omega = 4\pi.

These matrix elements are dimensionless.  Multiplication by
:math:`C_{\mathrm{sca}}/(4\pi)` gives the differential scattering cross-section
matrix, where :math:`C_{\mathrm{sca}}` is the scattering cross section.
For a size distribution, the cross sections and normalized phase matrix must
refer to the same distribution average.

See :doc:`user.tmatrix` for ARTS parameter names, units, and result types.
