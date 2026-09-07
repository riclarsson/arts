Radiative-transfer implementation
=================================

The notation below follows the linear-source derivation in
:doc:`concept.radiative_transfer`.

.. admonition:: Implementation note
   :class: tip

   It is very important to implement the expression for :math:`\Lambda_0`
   in a numerically stable way.  The expression above is not
   stable for small :math:`K_0 r` (i.e., large :math:`T_0`) *as written*.
   The key instability stems from the subtraction of two nearly equal
   terms in :math:`1 - T_0`.  The IEEE floating point standard
   provides a function :code:`expm1(x)`, which computes :math:`e^x - 1`
   in a numerically stable way for small :math:`x`.

   Likewise, the matrix expansion of :math:`1 - T_0` might be unstable.

   So we use a special solution implementing our own version of the reduced
   Cayley-Hamilton theorem to compute :math:`\Lambda_0` in a numerically
   stable way
   for matrices that conform to the propagation matrix notation in ARTS.
   This makes use of the inversion of :math:`K` to remove components from
   the expansion of the matrix exponential that would otherwise
   cause numerical instability.

