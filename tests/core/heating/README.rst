Heating-rate regression port
============================

The three arrays in ``arts2.py``'s ``REFERENCE_HEATING`` dictionary retain
the numeric values from the original ``*REFERENCE.xml`` files from ARTS 2
commit ``f90b319698617ad34c03dcae5238513d61cc4641`` (the local ARTS 2.6.18
audit tree). No reference values or original tolerances were changed.

``arts2.py`` reruns radiative transfer with ARTS3's ``cppdisort`` (and, for
the tropical case, ``cppvdisort``), then uses the maintained
``pyarts3.recipe.heating_rates`` integration and heating conversion functions.
It exercises both orders of angular and frequency integration. CMake's Python
test discovery registers this script and ``diagnostics.py``.

Fixed inputs
------------

The Garand regression loads the ARTS3 XML-data profile with::

    field = arts.AtmField()
    field.extendxml("planets/Earth/Garand/garand_profile.0.xml")
    field["O2"] = 0.2095
    field["N2"] = 0.7808

Altitude and pressure levels come from this field, and temperature is
interpolated from it at each ray sample. CO2 is already present in the profile.
The ARTS XML data path must provide this converted profile (the XML reader
also accepts its compressed ``.xml.gz`` form).

Only the absorption samples remain frozen pending the lookup-table port:
The ``EXTINCTION`` dictionary in ``arts2.py`` holds ``k0`` through ``k2``
for the three Garand ray subdivisions and ``extinction`` for the tropical
pressure levels. These are explicit NumPy array literals, with no separate
input fixture files. Decimal literals preserve every coefficient exactly. These are total scalar
gas extinction [1/m], summed over species, with axes (frequency, altitude).
They retain the original spectroscopy and H2O dilution treatment; the O2/N2
constants above do not regenerate absorption in this intermediate port.

All other intermediate properties are computed in the test:

* The frequency grid is the original ten-point linear grid, 3e11 to 3e13 Hz.
* The zenith quadrature and weights are double Gauss-Legendre with three
  nodes per hemisphere.
* Garand ray samples are subdivided to the original 10 km maximum path step;
  temperatures are evaluated from the loaded atmospheric field.
* The tropical pressure grid has 81 logarithmically spaced levels from
  101300 to 1 Pa. Temperature and altitude are interpolated in log pressure
  from ``planets/Earth/afgl/tropical/p.xml`` and ``t.xml``. The altitude
  coordinate supplies the original altitude profile. These profiles extend
  to 120 km; only the original test pressure range is used. Their additional
  upper levels do not change the values within that range. This atmosphere
  is distinct from Garand.
* The spherical radius is the original ARTS2 Earth setting, 6378100 m,
  specified in the test rather than stored in a file. It is retained rather
  than substituting a different modern Earth-radius convention.

No computed grids, quadrature, thermodynamic snapshots, NPZ files or binary
sidecars are stored. Reference heating rates are embedded in the test, with
exact numeric equality verified against the original XML values.
``tools/export_arts2_heating.py`` produces an ``extinction.py`` snippet with
the literal extinction and reference dictionaries for updating the test.
Reference values are read from the existing ARTS2 XML files, not recomputed. It
executes the old reference checks before extraction. This remains an RT and
heating regression until the absorption lookup handling is ported too.

Explicit adaptations
--------------------

* The Garand controlfile used first-order clear-sky emission: arithmetic mean
  of endpoint Planck radiances within each path segment. DISORT is configured
  with that constant layer source, with the original ray subdivision. ARTS3
  Planck evaluation is used; its constants are not patched to old values.
* The tropical Python test used CDISORT's linear-in-optical-depth Planck source.
  Its ``DTAUCPR > 1e-4`` rule is preserved: thinner layers retain the top
  Planck value and zero slope. Omitting that rule changes upper-atmosphere
  heating beyond the original tolerance. Both ARTS3 solvers receive the same
  source coefficients and non-scattering atmosphere.
* ARTS2 stored downward irradiance as negative; the recipe returns positive
  upward and downward hemispheric fluxes. The adapter uses ``up - down``.
  DISORT propagation cosines are mapped to ARTS viewing zenith angles.
* Original reference fluxes are at every atmospheric level, including TOA.
  Tests explicitly request all optical-depth boundaries from ``solver.flux``.
  Workspace ``DisortFlux`` instead omits TOA and evaluates at each layer's
  lower boundary. Its DFDT is **not** a heating rate in K/s.
* The old interior pressure stencil was a secant over neighboring levels,
  not a general nonuniform three-point derivative. Its lower boundary mixed
  upward flux at level 0 with downward flux at level 2. Its varying-gravity
  loop omitted the last level, leaving heating zero at TOA. Those behaviors
  exist only in the named ``legacy_heating`` test adapter. The public
  ``from_flux`` uses a nonuniform derivative and evaluates both boundaries.
* Original heat capacities remain 1006 and 1003.5 J/(kg K). The simple case
  uses 9.80665 m/s². The varying case retains the old equatorial gravity
  expression, including its extra 0.033895 m/s² centrifugal subtraction,
  and the original radius/altitude scaling. ARTS3's Earth gravity operator
  would be a physical-model change and is not substituted in reference checks.
* Tolerances remain absolute: 1e-9 K/s for Garand, 1e-6 K/s for both tropical
  results, and 1e-14 K/s for exchanging integration order. No relative
  tolerance is added. The unchanged zero at TOA is explicitly tested.

``diagnostics.py`` independently checks the supported formulas with a quadratic
flux on nonuniform pressure levels (both orderings and both boundaries),
isotropic hemispheres, frequency integration, analytic beam heating and
thermal cooling in both solvers, workspace DFDT sampling locations, spectral
extinction weighting, and invalid inputs. This separates preserved historical
behavior from the public scientific interface.

The port also corrects scalar DISORT's delta-M DFDT beam term. Total flux
retains the scaled beam after the direct/diffuse redistribution, so its
derivative must use that beam in mean intensity too. Previously DFDT used
the physical direct beam and disagreed with finite differences of net flux.
``cpp.fast.disort-test-delta-m-plus-settings`` checks all three scalar flux
interfaces against derivatives in physical optical depth. This correction
does not affect the non-scattering original heating-rate reference cases.
