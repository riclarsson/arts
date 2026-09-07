ARTS2 Monte Carlo source inputs
==============================

The original pressure, latitude and longitude grids are literal NumPy arrays
in ``mc_general_arts2.py``. Their XML values and printed precision are preserved
exactly, including the distinct near-zero longitude node. No grid XML files
are retained. The atmospheric preparation previously
stored here as ``TestMonteCarloDataPrepare.*`` is now performed in
``mc_general_arts2.py``; no generated fields or binary sidecars are retained.

Temperature, altitude and O2/N2/H2O abundances use the shared ARTS3 AFGL
tropical profiles in ``planets/Earth/afgl/tropical/``, interpolated in log
pressure onto the original MC pressure grid. Horizontal expansion is performed
at runtime. The top of atmosphere remains the altitude corresponding to the
original grid's minimum pressure (1 Pa), not the full 120-km AFGL extent.
No local copies of the tropical profiles are retained.

The original raw particle population is constructed directly in the test:
a zero array of shape (9, 10, 11), with indices ``[2:7, 2:8, 2:9]`` set to
34712.0922774 particles/m³. The original pressure/latitude/longitude coordinates
are given alongside it. This exactly reproduces the former ``pnd_field_raw.xml``
without storing the array. It is interpolated onto the original cloudbox
subgrid in log pressure, latitude and longitude, then embedded in a zero-valued
full atmospheric field.

The standard AFGL abundances replace the slightly different original MC
source abundances. Original MC I/Q reference values, four-standard-error
acceptance, seed and photon counts are unchanged. Temperature, altitude and
particle-density reconstruction were verified against the removed snapshots
to floating-point precision before this abundance substitution.

The large azimuthally random ice XML is replaced by an in-memory T-matrix
calculation in the test, which now requires ENABLE_TMATRIX.  The oblate
100-micrometre volume-equivalent-radius spheroid, aspect ratio 1.5, original
REFICE material values, frequency/temperature grids and angular grids are
preserved.  Full-precision SI wavelengths and quadrature replace historical
rounding.  Against the old XML, maximum differences normalized by each array's
peak were 0.180 ppm (phase), 0.089 ppm (extinction), and 0.744 ppm (absorption)
with the extended-precision backend.  No MC acceptance numbers were changed.
