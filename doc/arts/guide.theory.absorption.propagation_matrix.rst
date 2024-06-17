Propagation Matrix
==================

The propagation matrix conceptually describes how :artsdoc:`Spectral Radiance`
propagates through a system. The propagation matrix is a square matrix
with strict symmetries, and it has the form

.. math::

   \mathbf{K} = \left[ \begin{array}{rrrr}
        A & B & C & D \\
        B & A & U & V \\
        C &-U & A & W \\
        D &-V &-W & A
    \end{array} \right],

where
:math:`A` describes the total power reduction,
:math:`B` describes the difference in power reduction between horizontal and vertical linear polarizations,
:math:`C` describes the difference in power reduction between plus 45 and minus 45 linear polarizations,
:math:`D` describes the difference in power reduction between right and left circular polarizations,
:math:`U` describes the phase delay between right and left circular polarizations,
:math:`V` describes the phase delay between plus 45 and minus 45 linear polarizations, and
:math:`W` describes the phase delay between horizontal and vertical linear polarizations.
The unit of all of these is m :math:`^{-1}`.
