.. _dev-arts2notintended:

ARTS 2 Features Not Intended for ARTS 3
=======================================

This page records ARTS 2 interfaces and implementation patterns that should
not be restored merely to make a test or method list look complete.  The
scientific behavior behind an item can still be required; unresolved
capabilities belong in :doc:`dev.arts2missing` and must be implemented using
the ARTS 3 design.

ARTS 2 scripting language and executable
----------------------------------------

The ARTS 2 controlfile language is not an ARTS 3 compatibility target.  This
includes its parser, the ``arts`` controlfile executable, syntax-only checks,
controlfile-to-Python conversion, and command-line introspection tied to that
runtime.  The ARTS 2 CTests named ``cmdline.*``, ``converted.*``,
``TestWSMCalls``, and the controlfile-conversion example therefore do not
identify missing scientific functionality.

Use the Python interface for orchestration:

* Python ``for`` statements replace controlfile ``ForLoop`` helpers.
* Python assignment, copying, lists, and callables replace generic ``Copy``,
  ``AgendaCopy``, and ``ArrayOfAgenda`` convenience patterns where no native
  ARTS object semantics are required.
* Python loops and array assembly replace ``ybatch``, ``DOBatch``, and
  ``YCalcAppend`` control-language workflows.  If a large workload needs a
  parallel C++ kernel, track that kernel as a performance capability rather
  than restoring the language construct.
* Reusable high-level setups belong in documented ``pyarts3.recipe`` callables.

Python itself is fully supported and encouraged.  The restriction in the ARTS
3 design is narrower: do not dynamically add Python implementations that
pretend to be native workspace methods.  Prototype or extend native execution
through callback operators, agendas, and ordinary Python callables with clear
ownership and performance expectations.

Always-3D geometry, not atmosphere modes
----------------------------------------

ARTS 3 has a spheroidal reference ellipsoid and an always-three-dimensional
geometric view.  Positions use altitude above the ellipsoid, geodetic latitude,
and longitude; directions use local zenith and azimuth angles.

Do not port the following ARTS 2 mechanisms literally:

* ``AtmosphereSet1D`` and ``AtmosphereSet2D`` as global model modes;
* general/master ``z_grid``, ``lat_grid``, and ``lon_grid`` workspace
  variables that determine every field's shape; or
* separate 1-D, 2-D, and 3-D propagation-path implementations and tests.

Local grids can still be method inputs, and individual atmospheric quantities
can be constant or low-dimensional.  Tests should exercise those quantities
through the common three-dimensional coordinate and path representation.

No global Stokes or Jacobian switches
-------------------------------------

The ARTS 2 ``stokes_dim`` workspace variable is not intended to return.  ARTS
3 data types represent full Stokes quantities.  Sensors select polarization
through sparse Stokes-projection weights, and a solver may expose a local
computational option when that is useful; neither case requires a global
dimension switch.

Likewise, do not restore a global ``jacobian_do`` flag.  ARTS 3 represents
requested derivatives with ``jac_targets`` and their setup/finalization.  The
early design notes used the tentative name ``jacobian_setup``; the current API
name is authoritative.  An empty target set or a method's explicit inputs
should determine whether derivative work is needed.  A limited ``do_jac``
switch survives in the OEM/inversion wrapper and can clear ``jac_targets``;
it is not a deep, global switch for arbitrary forward-model methods and should
not grow into one.

No formal cloudbox
------------------

The ARTS 2 ``cloudbox_on``, ``cloudbox_limits``, checked flags, cloudbox field,
and cloudbox interpolation/setup machinery are not part of the ARTS 3 model.
Scattering solvers should consume explicit atmospheric fields, scattering
species, paths/domains, and solver settings.

This decision does not remove the need for scattering solvers.  In particular,
DOIT remains a missing numerical capability even though its old cloudbox API
must not be ported.  ``TestCloudboxAuto`` and
``Testcloudbox_fieldInterp2Azimuth`` are therefore retired interface tests;
their useful interpolation or boundary behavior must be tested through the new
solver's actual inputs.

Composite fields and fewer workspace variables
----------------------------------------------

ARTS 3 intentionally favors fewer workspace variables, more explicit generic
inputs/outputs, and stronger high-level types.  Dimensions and invariants
should be fixed by the type when practical.  Do not split ``AtmField``,
``SurfaceField``, ``SubsurfaceField``, model state, or scattering state back
into the old global grid plus tensor collection merely to reproduce an ARTS 2
method signature.

Small named methods are still encouraged when they clarify the data flow.
Prefer strict, composable methods and static meta-methods over optional branches
that silently change user-facing behavior.

Sensor representation
---------------------

The ARTS 2 collection of global sensor position, line-of-sight, frequency-grid,
response, polarization, and normalization workspace variables has been
replaced by ``measurement_sensor``.  It is a list of observation elements;
each element references frequency and position/line-of-sight grids (which can
be shared) and stores sparse Stokes-vector projection weights needed to produce
a measurement value.  Metadata is associated with observation ranges rather
than a parallel set of loosely coupled workspace variables.

Port channel definitions and verified measurement behavior, not the old sensor
workspace-variable plumbing.  Prefer ``measurement_sensorFromPredefined`` or
``pyarts3.arts.sensor.Builder`` to instrument-specific controlfile frameworks
such as MetMM.

Scattering data representation
------------------------------

The ARTS 2 scattering species, particle-number-density fields, scattering data
arrays, and agenda setup were marked for a complete redesign.  ARTS 3 uses
typed scattering-species properties, particle habits and PSDs, gas or
parametric scatterers, and an array of scattering species that evaluates bulk
properties at an atmospheric point.

Port a PSD, habit, optical-property model, or solver result against this API.
Do not restore the old ``ScatSpeciesInit``/``pnd_field``/``scat_data`` object
layout or assume that a passing ARTS 2 setup sequence defines the new
interface.

Failure and validation behavior
-------------------------------

ARTS 2 relied heavily on eager ``*_checkedCalc`` flags and large up-front
consistency passes.  ARTS 3 methods and agendas should work or fail naturally,
validate values when they are used, and check dimensions shared by multiple
inputs.  Only ``std::runtime_error`` should escape C++ ARTS code, with the
binding presenting an appropriate Python exception and preserving the error
chain.

Do not port silent atmosphere-class switches, class-ignore flags, or
solver-side zeroing.  If a calculation needs a quantity set to zero, the caller
should construct that field explicitly.  A solver must not mutate or zero
unrelated atmospheric quantities; it may ignore quantities outside its
contract without treating their mere presence as an error.

ARTS 2 checks retired by these decisions
----------------------------------------

The following test families should be mapped to ARTS 3 behavior, not copied:

* ``TestPpath1D`` and ``TestPpath2D``: cover their geometry through the common
  3-D path implementation.
* ``TestForloop``, ``TestAgendaCopy``, ``TestArrayOfAgenda``,
  ``TestDOBatch``, and ``TestYCalcAppend``: use Python orchestration unless an
  underlying numerical or parallel primitive is absent.
* ``TestCloudboxAuto`` and cloudbox interpolation tests: test explicit solver
  inputs; do not restore cloudbox state.
* ``TestHSE`` and the old pressure-grid regridding family: test
  ``atm_fieldHydrostaticPressure`` and independent field interpolation rather
  than global pressure/altitude grid mutation.
* instrument controlfiles: retain accepted response data and forward-model
  references through observation elements/builders, not instrument-specific
  workspace-variable frameworks.
* ``jacobianAdjustAndTransform``: transformations belong to individual typed
  targets and model-state mappings, not one late global adjustment pass.
* ``atmfields_checkedCalc``, ``atmgeom_checkedCalc``,
  ``cloudbox_checkedCalc``, and ``sensor_checkedCalc``: replace checked-state
  flags with local contracts and use-time errors.

Legacy choices that remain
--------------------------

Some older mechanisms remain deliberately or transitionally.  Their presence
is not permission to expand them indiscriminately.

XML file format
~~~~~~~~~~~~~~~

The ARTS 3 design explicitly keeps XML "for now" and prefers storing data in
its natural form.  XML I/O and selected old-format readers therefore remain a
supported legacy choice.  New data should still use typed structures and a
documented schema; compatibility parsing should not dictate the in-memory
design.

``f_grid``
~~~~~~~~~~

The ``f_grid`` working variable and its name survived.  It is a local grid for
the methods that consume it, not a global/master grid system.  Renaming it to
``frequency_grid`` was left as an open question, not a settled migration task.
Do not perform a partial rename without a project-wide decision.

Workspace methods, variables, agendas, and arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These concepts remain, but with narrower roles.  Operators are for focused
pure work; agendas represent impure or configurable ideas; static meta-methods
compose smaller methods.  ``Array<>`` storage remains appropriate when each
element acts independently and should not be forced into a matpack tensor.

Loader-level missing-data options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The current ARTS 3 implementation exposes explicit ingestion policies such as
``missing_is_zero``, ``ignore_missing``, and ``name_missing``.  They operate at
the data-loading boundary and do not contradict the rule that a required model
quantity must be present when used.  Do not use these policies as precedent
for solver-side silent defaults or broaden them beyond explicit ingestion
contracts without a design decision.

OEM model-state update cost
~~~~~~~~~~~~~~~~~~~~~~~~~~~

ARTS 3 OEM can update the full model state on the fly, even though copying or
updating the composite atmosphere, surface, and line data can be more expensive
than the narrower ARTS 2 tensor updates.  The design notes explicitly accept
that trade-off.  Optimize it without changing the model-state semantics unless
measurements show that a separate mode is necessary.

Residual control-language tooling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Editor grammar assets remain under ``tools/arts-control-lang``.  They are
legacy tooling, not a commitment to ship or test the ARTS 2 parser/runtime.
Remove or repurpose them separately; do not infer language compatibility from
their presence.

Stale optional-feature hooks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The current tree still contains ``ENABLE_RT4`` and ``ENABLE_TMATRIX``
configure definitions, conditional link clauses for the two libraries, and a
T-matrix enum value in
``src/core/scattering/optproperties.h``.  No corresponding ARTS 3
implementations were found in this audit.  These hooks are stale build/API
residue, not evidence that the capabilities are already supported.  Remove
them when safe, or reconnect them only as part of an ARTS 3-native port.

Open questions, not design decisions
------------------------------------

The early ARTS 3 notes leave several questions open.  Do not place either side
of these questions in the missing-feature backlog without a new decision:

* whether ``f_grid`` should eventually be renamed;
* how a solver-local Stokes-size optimization should be expressed;
* whether propagation paths should branch and how crossings are requested;
* whether every internal radiative-transfer path should always use the full
  Planck function;
* whether Kokkos should replace OpenMP;
* how sensor normalization should be explained; and
* how Jacobian target and covariance information should be unified.

Port-review checklist
---------------------

Before porting an ARTS 2 test or method, ask:

#. Is a physical, numerical, data, or performance capability actually absent?
#. Can ordinary Python or a maintained recipe express it clearly?
#. Which ARTS 3 field, sensor, scattering, path, or Jacobian type owns it?
#. Does the proposed port accidentally restore global grids, cloudbox state,
   ``stokes_dim``, checked flags, or controlfile-language conveniences?
#. Can a compact active regression compare a meaningful result without a
   runtime dependency on the temporary ARTS 2 source tree?

If only the old spelling or orchestration is missing, do not port it.  If the
result itself cannot be produced, add it to :doc:`dev.arts2missing` with the
evidence and an ARTS 3-native acceptance test.
