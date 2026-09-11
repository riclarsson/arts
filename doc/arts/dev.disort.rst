DISORT implementation and validation
=====================================

The scalar and vector solvers are implemented in ``src/core/disort-cpp``.
See :doc:`concept.disort` for the equations and notation used below.

Implementation array layouts
============================

The similarly named arrays in the two cores do not always store the same
factorization:

.. list-table::
   :header-rows: 1
   :widths: 24 35 41

   * - Mathematical object
     - CPP-DISORT
     - VDISORT
   * - Layer transport matrix :math:`\boldsymbol A`
     - Implicit in ``D_pos``, ``D_neg``, ``apb``, and ``amb``; only the
       :math:`N\times N` reduced product ``sqr`` is diagonalized
     - Local :math:`4N_q\times4N_q` real ``A_real`` is copied to a complex
       matrix and diagonalized
   * - Propagation constants :math:`\boldsymbol K`
     - ``K_collect[m, layer, eigen]`` with negative half followed by its
       positive partners
     - ``K_collect[alpha, m, layer, eigen]`` after complex lexicographic sort
   * - Eigenvectors :math:`\boldsymbol G`
     - ``G_collect[m, layer, state, eigen]``
     - ``G_collect[alpha, m, layer, state, eigen]``
   * - Modal constants :math:`\boldsymbol c_\ell`
     - The band-solve ``RHS`` is overwritten by the constants and copied to
       ``C_collect[m, layer, eigen]``
     - The complex band-solve ``rhs`` is copied into ``GC_collect``
   * - ``GC_collect``
     - Caches :math:`G_{qe}c_e` for ordinary exponential evaluation; the
       conservative pair is evaluated from ``C_collect`` and its coupled
       two-column basis
     - Despite the name, stores only :math:`c_e`; multiplication by
       :math:`G_{qe}` occurs during field reconstruction
   * - Beam particular solution :math:`\boldsymbol B`
     - ``B_collect[m, layer, stream]``
     - ``B_collect[alpha, m, layer, stream]`` of Stokes vectors
   * - Polynomial particular solution
     - Coefficients are cached in
       ``source_collect[layer, stream, power]``; ``SRC0``, ``SRC1``, and
       ``SRCB`` retain boundary values used by the ordinary fast assembly
     - Coefficients are Stokes vectors in
       ``source_collect[alpha, m, layer, stream, power]``
   * - Layer-bottom total field
     - ``um[layer, m, stream]``
     - ``um[layer, alpha, m, stream]`` of Stokes vectors

In particular, ``um`` is a cached total quadrature field, not an eigenvector or
a modal-coefficient array.

Numerical behavior and limitations
**********************************

The following details are intentional and should be considered when changing
or comparing the cores:

* **Conservative scattering.**  Exact :math:`\omega=1` gives a defective zero
  pair in the zeroth cosine mode.  Both cores keep the input albedo unchanged
  and represent the physical energy pair by the constant/linear centered basis
  above.  The same basis is selected through the near-conservative interval
  :math:`[1-10^{-8},1]`; ordinary modes retain the fast anchored-exponential
  path.  VDISORT stabilizes the one physically required energy-conservation
  pair.  A contrived Mueller operator with additional independently conserved
  polarization quantities can contain more zero pairs and is not covered by
  this single-pair representation.  Generalizing this is intentionally deferred
  until a physically occurring scattering model requiring it is identified.
* **Complex VDISORT modes.**  A real physical solution can use complex
  eigenpairs.  The reconstructed imaginary part is required to cancel to a
  relative tolerance of about :math:`2\,10^{-8}` before the real part is
  returned.  Failure indicates an ill-conditioned or incorrectly assembled
  system, not a physical complex radiance.
* **VDISORT spectral split.**  Anchoring requires half the non-neutral modes to
  propagate toward each boundary.  After sorting, the implementation checks
  this sign split.  A relative neutral band of :math:`10^{-10}` times the
  eigenspectrum scale permits mathematically zero or nearly imaginary modes to
  occupy either half; a non-neutral sign imbalance or a mode on the wrong side
  is rejected before boundary assembly.

* **No absorption cutoff.**  An absorption-optical-depth shortcut is not
  implemented.  The complete physical solution is evaluated instead.
* **No special-boundary shortcut.**  ``IBCND=1`` albedo/transmission is not a
  separate algorithm here.  The same quantities are obtained from the regular
  boundary-value solution.


* **Update cost.**  Constructing or updating either solver performs the
  eigendecompositions and global boundary solves.  Evaluation at cached
  quadrature points is much cheaper.  VDISORT intentionally contains no
  OpenMP parallel region in these core algorithms; ARTS normally parallelizes
  independent frequencies outside the solver.
* **Thread use.**  A fully constructed solver may be shared for read-only
  evaluation when each caller owns its scratch objects.  Updating the solver
  or sharing scratch storage concurrently is unsafe.

Validation scope
****************

The scalar core is checked against published numerical reference values for
the supported plane-parallel test problems, including arbitrary angles, fluxes,
DFDT, delta-M corrections, physical BRDFs, conservative scattering, multilayer
continuity, and delta-M-plus.  The pseudo-spherical test case is intentionally
unsupported.

The strict scalar embedding of VDISORT is checked against the same supported
reference problems, with :math:`Q=U=V=0` asserted.  These tests establish the
normalization and scalar limit of the vector equations.  Analytic polarized
two-stream tests additionally cover :math:`I/Q` coupling, complex :math:`U/V`
eigenpairs, a polarized direct beam, polarized absorption, vector internal
sources, and a polarized reflecting boundary.  These focused tests do not
replace comparison with an independent general-purpose polarized reference,
particularly for reference-plane rotations, many-stream Mueller problems, and
genuinely polarized IMS/TMS corrections.

Caller conventions are described in :doc:`user.disort`.
