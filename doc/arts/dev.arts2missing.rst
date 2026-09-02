.. _dev-arts2missing:

ARTS 2 Feature-Port Status
==========================

This is a living audit of scientific and numerical functionality that existed
in ARTS 2 and is not yet fully accounted for in ARTS 3.  It is not an API
compatibility promise and it is not a list of ARTS 2 workspace methods to copy.
The companion :doc:`dev.arts2notintended` page records deliberate redesigns
and features that must not be ported literally.

The main rule for this audit is to track capabilities, not controlfile syntax.
An ARTS 2 helper is not missing merely because it has no ARTS 3 workspace
method.  Plain Python, NumPy, a Python recipe, or a composition of smaller ARTS
3 methods is an acceptable and often preferred replacement.  A feature remains
interesting when the underlying physical result, numerical algorithm, data
model, or practical workflow cannot be produced through those interfaces.

Status terms
------------

**Missing**
  No ARTS 3 implementation or usable replacement was found.  A port or an
  explicit decision not to port is required.

**Coverage gap**
  The implementation or its building blocks exist, but no active regression
  demonstrates parity with the corresponding ARTS 2 result.

**Python decision**
  The ARTS 2 operation is suitable for Python.  Decide whether a maintained
  recipe/helper is needed; do not automatically recreate a workspace method.

**Deferred**
  The configured test inventories did not exercise the feature.  A later
  method-, data-, and optional-build audit must classify it.

Initial audit snapshot
----------------------

This first snapshot was made on 2026-09-02 from ARTS 3 branch ``ai-port`` at
``566654276`` and an ARTS 2.6.18 source/build tree.

* The configured ARTS 3 build registers 190 CTests: 30 Python examples, 128
  Python regression tests, and 32 C++ tests.
* The configured ARTS 2 build registers 201 CTests.  Many are paired
  ``ctlfile.*`` and ``converted.*`` executions of the same case and must not be
  counted as two capabilities.
* ``cmake --build build --target check -j 12`` selected 110 ARTS 2 checks.
  It excludes converted duplicates and the ``slow``, ``xmldata``, ``nocheck``,
  and ``planettoolbox`` groups.  108 passed.
  ``pyarts.fast.artscomponents.lookup.TestLookup_h2o`` and the aggregate
  ``pytest`` check segfaulted in the old Python extension; those two
  environment-sensitive failures are not evidence of ARTS 3 feature gaps.
* The runs used the ARTS XML and catalogue data trees identified by CMake.
  Both builds were configured as ``Release``; the ARTS 2 build had Fortran and
  NetCDF disabled.  Test counts are configuration-dependent, especially for
  LGPL, SHTNS, FASTEM, RT4, T-matrix, and NetCDF options.

The classification below combines that CTest comparison with a narrow source
audit where a test name alone could not distinguish a missing feature from an
ARTS 3 replacement.  It does not yet compare the complete method inventories.

Confirmed implementation gaps
-----------------------------

DOIT scattering solver
~~~~~~~~~~~~~~~~~~~~~~

**Status: Missing.**  ARTS 2 registers ``TestDOIT``,
``TestDOITaccelerated``, ``TestDOITprecalcInit``,
``TestDOITsensorInsideCloudbox``, and the XML-data-dependent
``TestDOITpressureoptimization``.  The selected fast DOIT checks passed in the
ARTS 2 audit.  No DOIT implementation, workspace method, example, or test was
found in ARTS 3.

ARTS 3 DISORT and Monte Carlo calculations are useful alternatives, but they
are not a deterministic polarized DOIT implementation.  A DOIT port must use
the ARTS 3 atmospheric-field and scattering-species models.  It must not
restore the formal ARTS 2 cloudbox, global spatial grids, or global
``stokes_dim`` workspace variable merely to preserve the old setup sequence.

FASTEM surface model
~~~~~~~~~~~~~~~~~~~~

**Status: Missing; ARTS 3 surface-boundary design required.**  FASTEM is an
intended ARTS 3 capability, but the current surface interface cannot represent
the complete model.  FASTEM can provide separate emissivity and reflectivity
and includes an atmospheric-transmittance correction, whereas the existing
closed-surface reflectance path applies Kirchhoff consistency and has no
explicit atmospheric dependency.

Design a surface-boundary contract that represents emission, reflection, and
the required atmospheric coupling explicitly.  The implementation must not
silently discard FASTEM outputs or weaken the consistency guarantees of the
existing reflectance agenda.  The surviving ``ENABLE_FASTEM`` configure and
link hooks are legacy integration residue; they do not provide this interface
or an ARTS 3 implementation.

Microwave gas refractivity
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Status: Missing.**  The selected ARTS 2 ``TestRefractPlanets`` check
exercises both the Earth-specific and general microwave gas refractivity
models (``GasMicrowavesEarth`` and ``GasMicrowavesGeneral``).  ARTS 3 has the
refractive ray-path machinery and a tested Harvey-98 visible/near-infrared
water/steam refractive-index model, but no implementation of these microwave
dry/moist-gas formulae was found.

Port the supported formulae as propagation-property inputs to the common ARTS
3 path model; do not restore atmosphere-dimensional modes or the ARTS 2
``refr_index_air_agenda`` setup interface.  Add Earth and non-Earth reference
cases, including a geometric-path comparison.

Heating-rate diagnostic
~~~~~~~~~~~~~~~~~~~~~~~

**Status: Missing packaged functionality; Python decision.**  ARTS 2
``TestHeatingRates`` and ``Test_HeatingRate`` calculate a radiance/irradiance
field and then derive heating rate from flux divergence.  ARTS 3 provides and
tests spectral radiance fields, spectral flux profiles, DISORT flux fields,
frequency integration, and the ``AtmosphericFlux`` and
``SpectralAtmosphericFlux`` Python recipes.  It does not provide a heating-rate
method or recipe.

The missing final calculation is well suited to Python.  Prefer a documented,
unit-aware recipe with an ARTS 2 reference regression unless a C++ kernel is
needed for performance or for use inside an agenda.

Confirmed gaps with unresolved port intent
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Status: Missing; decision required.**  Source inspection found several ARTS
2 capabilities that were excluded from the configured baseline by optional
dependencies or were otherwise not registered there.  Their absence from the
110-test ``check`` selection is not evidence that they were ported:

* RT4 and its cross-solver/hybrid comparison are absent.  DISORT and the ARTS
  3 Monte Carlo solvers are not numerical parity for RT4.
* The Mishchenko T-matrix single-scattering-data generator is absent.  Reading
  or using existing single-scattering data is a different capability.
* The NetCDF/libRadtran ``WriteMolTau`` exporter has no packaged equivalent;
  this is likely best implemented as a Python/xarray exporter if still needed.

Decide whether each item remains supported science or a maintained workflow.
Accepted items need ARTS 3-native interfaces and reference tests; rejected
items should move to :doc:`dev.arts2notintended` with the reason and date.

Active coverage gaps
--------------------

These are not all missing implementations.  They are cases where the first
audit cannot yet make a defensible parity claim.

Line-by-line derivatives
~~~~~~~~~~~~~~~~~~~~~~~~

**Status: Coverage gap, high priority.**
``tests/core/lbl/voigt_derivatives_perturbed.py`` is registered as a CTest, but
its substantive body is currently below ``if 0`` and does not run.  The
Jacobian comparisons in ``voigt_lte_pure.py`` and ``voigt_lte_mem.py`` are also
commented out.  ARTS 2 has active shape-derivative and VMR-derivative checks.

Re-enable or replace these tests with small finite-difference cases that cover
temperature, species amount, wind, magnetic field, line center/strength, and
the supported line-shape parameters.  A test that imports successfully while
executing no assertions must not count as port evidence.

CIA derivatives
~~~~~~~~~~~~~~~

**Status: Coverage gap.**  ARTS 2 ``TestCIADerivs`` checks temperature, VMR,
and wind sensitivities.  ARTS 3 implements CIA Jacobian output and tests CIA
values and file I/O, but no equivalent finite-difference derivative regression
was found.  Add a compact test for every CIA target that the implementation
claims to support.

Transmission including particulate extinction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Status: Building blocks present; coverage and dependency gap.**  ARTS 2
``TestTransmissionWithScat`` covers a refracted three-dimensional, full-Stokes
transmission path through particles.  ARTS 3 can compose
``spectral_propmat_scat_pathFromPath``,
``spectral_propmat_pathAddScattering``,
``spectral_tramat_pathFromPath``, and
``spectral_radCumulativeTransmission``.  No current test composes those pieces
into a passive-transmission regression.  Add a geometric-path particle test
now, without reintroducing cloudbox state.  Exact parity with the ARTS 2 case
also depends on the missing microwave gas refractivity model recorded above.

Absorption-only particle radiative transfer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Status: Building blocks present, coverage gap.**  The slow ARTS 2
``TestAbsParticle`` check validates end-to-end thermal radiative transfer for
absorption-only particles against a reference.  ARTS 3 has the relevant
particle extinction/absorption-vector and thermal source-correction machinery,
but no active equivalent numerical assertion was found.  Add a compact thermal
reference case; a plotting example or a zero-absorption gas-scattering test is
not sufficient parity evidence.

Hydrostatic equilibrium
~~~~~~~~~~~~~~~~~~~~~~~

**Status: Replacement present, coverage gap.**  ARTS 2 ``TestHSE`` adjusts an
altitude field as a function of pressure with ``z_fieldFromHSE``.  ARTS 3 has
independent altitude-indexed atmospheric quantities and provides the
design-native reverse representation, ``atm_fieldHydrostaticPressure``, which
constructs pressure as a function of altitude.  The old direction is not a
literal port target.  No active ARTS 3 regression was found for the replacement
method; add one covering its supported options and boundary behavior.
ARTS 3 implements HSE-coupled path-length temperature derivatives through the
``hse_derivative`` option, but no current test exercises it.  Add an
analytic-versus-perturbed case covering the full-Stokes/Zeeman transmission
interaction represented by ARTS 2 ``TestTjacStokes4_transmission``.

Field regridding
~~~~~~~~~~~~~~~~

**Status: Replacement present, coverage gap.**  The old
``GriddedFieldPRegrid``, ``AtmFieldPRegrid``, and master-grid refinement family
does not fit the ARTS 3 field model.  ARTS 3 exposes ``AtmData.regrid`` and
``AtmField.regrid`` to Python and lets atmospheric quantities interpolate
independently.  Generic interpolation is tested, but no direct atmospheric
field regrid test was found.  Add binding-level tests rather than restoring the
pressure-grid workspace-method family.

Line cutoffs and line mixing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Status: Present, focused coverage incomplete.**  ARTS 3 line-by-line code
implements per-band cutoffs, and current camera/lookup tests use a fixed cutoff.
The ARTS 2 Python cutoff check sweeps cutoff values against reference radiances;
there is no equally focused ARTS 3 regression.  Add a small sweep, including
cutoff Jacobians where supported.

ARTS 3 ECS O2 and CO2 tests exercise the new line-mixing representation, so the
old ``TestMolOxyAdaptation`` API is not a port target.  Strengthen the ECS tests
with numerical reference or perturbation assertions before declaring complete
scientific parity.

Predefined instruments and end-to-end radiances
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Status: Seven definitions present with incomplete end-to-end coverage; other
presets missing pending decision.**  ARTS 3 has compiled definitions and setup
tests for AMSU-A, AMSU-B, ICI, HIRS, AVHRR, MVIRI, and SEVIRI.
``tests/core/sensor/predefined_sensor.py`` checks channel counts, normalization,
and AMSU-A horizontal/vertical Stokes-projection weights.  A separate custom
AMSU-A-like Builder example checks simulated channel radiances, but does not
validate the compiled predefined AMSU-A response end to end.

The ARTS 2 ICI, HIRS, AVHRR, MVIRI, SEVIRI, and MetMM tests run complete
forward/batch calculations.  Add compact end-to-end reference cases for the
remaining accepted definitions.  In particular, record whether AVHRR's two
zero-response solar rows are an intentional limitation.  Do not port the
MetMM controlfile framework or its dedicated grid/response methods; use the
predefined sensor API or ``pyarts3.arts.sensor.Builder``.  The ARTS 2 MetMM
presets ATMS, DEIMOS, HATPRO, ISMAR up/down, MARSS, MHS, MWHS-2, and SAPHIR are
not present in the ARTS 3 predefined-sensor enumeration.  Decide which data
definitions to port individually rather than restoring the framework wholesale.

Raw calibration and time-series corrections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Status: Python decision.**  ARTS 2 registers ``raw/calib.py`` and
``raw/corr.py`` for hot/cold calibration, timestamp sorting, a simple
tropospheric correction, and time averaging.  No packaged ARTS 3 equivalents
were found.  These array/time-series operations are natural Python operations.
Determine whether users need maintained recipes with the exact ARTS 2
semantics; otherwise document a Python replacement and do not add workspace
methods.

Additional scientific parity checks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Status: Coverage gaps or Python decisions.**  The initial test comparison
also identified the following narrower cases.  They should become compact
regressions or receive an explicit explanation of why existing coverage is
sufficient:

* CIA has value and I/O tests but lacks both the derivative check above and an
  end-to-end CIA radiance reference.
* Faraday rotation is implemented, but no direct ARTS 2-equivalent propagation
  or radiance regression was found.
* The Mars catalogue test covers end-to-end Earth-versus-Mars isotopologue-ratio
  changes and their effect on clear-sky radiance.  Wider ARTS 2 planetary
  isotopologue cases remain unmatched.
* Cross-section-fit point evaluation is tested; an end-to-end XFIT radiance
  reference is not.
* The predefined-sensor test already checks horizontal/vertical Stokes
  projection.  Add scan-angle polarization rotation and an end-to-end polarized
  radiance case.
* Geographic surface-type dispatch and scaled-species propagation helpers
  appear expressible by composing current field/propmat operations in Python.
  Decide whether stable recipes are needed instead of adding compatibility
  workspace methods.

Deferred beyond this first pass
-------------------------------

The ARTS 2 tree contains further features that were not registered in the
configured CTest inventory or were hidden behind disabled optional builds.
The confirmed absences found in the narrow source audit are classified above.
Covariance-only controlfiles and some antenna and single-polarization cases
still need enough inspection to distinguish missing science from Python
orchestration, replacement APIs, external-library responsibilities, and
deliberate removals.

A later audit should compare:

#. the complete ARTS 2 and ARTS 3 public method/group inventories;
#. optional build configurations and their tests;
#. catalogue/XML data formats and model data availability;
#. examples and user workflows that were never registered as tests; and
#. numerical reference outputs, not merely method-name similarity.

Updating this page
------------------

This page is only a backlog of work that remains.  Remove an entry as soon as
its implementation and required coverage are complete; do not retain completed
features as historical status records.  For each remaining entry, record the
ARTS 2 capability and evidence, the applicable ARTS 3 design constraint, and
the remaining action.  A port is complete only when the result is usable
through a supported C++ or Python interface and an active regression checks
meaningful physics or numerics.  If Python orchestration is sufficient, add a
test and documentation or a recipe; do not recreate the ARTS 2 scripting
language.
