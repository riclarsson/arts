Radiative-transfer layer options
================================

Select the within-layer approximation with
:attr:`~pyarts3.workspace.Workspace.rte_option`.

.. list-table:: Radiative-transfer layer options
   :header-rows: 1
   :widths: 18 30 52

   * - Option
     - Within-layer model
     - Step
   * - ``constant``
     - Endpoint-average propagation matrix and source
     - Ordinary matrix exponential with a layer-average source.
   * - ``lintau``
     - Constant propagation matrix and linear source
     - Uses the linear-source operator :math:`\Lambda` derived in :doc:`concept.radiative_transfer`.
   * - ``linprop``
     - Linear propagation and linear source
     - Uses the exact scalar linear-propagation source integral.  Polarized
       transfer uses endpoint-average transmission with a commutator-free
       augmented-source correction.
   * - ``magop``
     - Linear propagation matrix and layer-average source
     - Uses the second-order Magnus exponent and the ordinary source step.
   * - ``magop_linsrc``
     - Linear propagation matrix and linear source
     - Adds Magnus-ordered transmission to the augmented linear-source
       operator by retaining the first propagation-matrix commutator.

All five options propagate analytical derivatives of their transmittance
and, where applicable, source operators.  ``magop`` and ``magop_linsrc``
are most useful when the polarized propagation matrices at the layer
endpoints do not commute.  For scalar transfer, ``linprop`` evaluates the
linear-propagation source integral exactly.  For polarized transfer,
``linprop`` retains the ordinary endpoint-average transmission, whereas
``magop_linsrc`` also includes the Magnus ordering correction.

See :doc:`concept.radiative_transfer` for the equations behind these options.
