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

ECS line-by-line derivatives
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Status: Missing; design deferred.**  LTE Voigt, mirrored LTE Voigt, and
line-NLTE derivatives now have active finite-difference coverage for their
supported atmospheric, spectroscopic-line, and line-shape targets.  ECS
derivatives remain intentionally outside that implementation: an ECS target
must also account for the ECS parameter set, and its derivative representation
should be designed before adding a partial interface.  Do not treat the
non-ECS derivative support as evidence that ECS derivatives are complete.

Confirmed gaps with unresolved port intent
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Status: Missing; decision required.**  Source inspection found several ARTS
2 capabilities that were excluded from the configured baseline by optional
dependencies or were otherwise not registered there.  Their absence from the
110-test ``check`` selection is not evidence that they were ported:

* RT4 and its cross-solver/hybrid comparison are absent.  ARTS 3 already has
  polarized VDISORT and Monte Carlo solvers; the remaining question is whether
  the independent RT4 algorithm and its legacy comparisons justify a port.
* General azimuthal orientation averaging and a native ARO T-matrix habit
  factory remain unconnected.  The vertically aligned, axisymmetric particle
  in ``tests/core/scat/mc_general_arts2.tmatrix.py`` is already generated with
  ``tmatrix.fixed_batch`` and consumed through the legacy ARO adapter.  That
  case needs no orientation averaging and is not an outstanding port.  Scope
  any further work to orientation distributions and native integration beyond
  this case (see :doc:`dev.tmatrix`).
* The NetCDF/libRadtran ``WriteMolTau`` exporter has no packaged equivalent;
  this is likely best implemented as a Python/xarray exporter if still needed.

Decide whether each item remains supported science or a maintained workflow.
Accepted items need ARTS 3-native interfaces and reference tests; rejected
items should move to :doc:`dev.arts2notintended` with the reason and date.

Active coverage gaps
--------------------

These are not all missing implementations.  They are cases where the first
audit cannot yet make a defensible parity claim.

CIA derivatives
~~~~~~~~~~~~~~~

**Status: Coverage gap.**  ARTS 2 ``TestCIADerivs`` checks temperature, VMR,
and wind sensitivities.  ARTS 3 implements CIA Jacobian output and tests CIA
values and file I/O, but no equivalent finite-difference derivative regression
was found.  Add a compact test for every CIA target that the implementation
claims to support.

Transmission including particulate extinction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Status: Building blocks present; coverage gap.**  ARTS 2
``TestTransmissionWithScat`` covers a refracted three-dimensional, full-Stokes
transmission path through particles.  ARTS 3 can compose
``spectral_propmat_scat_pathFromPath``,
``spectral_propmat_pathAddScattering``,
``spectral_tramat_pathFromPath``, and
``spectral_radCumulativeTransmission``.  No current test composes those pieces
into a passive-transmission regression.  Add a geometric-path particle test
now, without reintroducing cloudbox state.

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

**Status: HSE temperature coupling missing; path-length derivative validation
missing.**  The pressure-construction replacement is implemented and tested in
``tests/core/atm/hse.py``.  ARTS 2 ``z_fieldFromHSE`` adjusts altitude at fixed
pressure levels; ARTS 3 ``atm_fieldHydrostaticPressure`` constructs pressure at
fixed altitudes.  Restoring the old altitude-field representation is not a
port target.

The remaining derivative work must distinguish those two constraints:

* ``hse_derivative`` enables an existing endpoint path-length correction in
  transmission and emission calculations when a temperature Jacobian target
  is present.  It is disabled by default and has no active regression.  This
  correction describes layer expansion at fixed pressure levels; it is not
  the temperature derivative of ``atm_fieldHydrostaticPressure``.  Define its
  supported perturbation model and check it against finite differences, first
  for an isothermal column and then for nonuniform temperatures.  The
  full-Stokes/Zeeman interaction in ARTS 2
  ``TestTjacStokes4_transmission`` also remains unvalidated.
* The altitude-indexed ``Atm::HydrostaticPressure`` callable stores pressure
  and gradient values computed at construction.  Temperature updates do not
  rebuild it, and no temperature-to-HSE-pressure Jacobian is supplied.
  Supporting that coupling requires an explicit dependency and propagation
  of the induced pressure changes through the radiative-transfer calculation.
  Validate it by perturbing temperature and rebuilding HSE pressure at fixed
  altitudes.  Detecting the callable's type alone must not automatically
  enable the existing path-length correction.

The pressure-construction test uses a fixed specific gas constant.
Composition-derived gas constants and nonuniform horizontal fields still lack
focused coverage; these are validation gaps, not missing HSE solvers.

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
The LTE and NLTE finite-difference derivative tests also exercise pressure
derivatives with a line cutoff.  What remains unmatched is the ARTS 2 Python
cutoff sweep against reference radiances.  Add a focused radiance sweep; do
not classify the already tested pressure-derivative branches as missing.

ARTS 3 ECS O2 and CO2 tests execute and plot the new line-mixing representation, so the
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
