Absorption
**********

Absorption is the physical process that reduces 
the :artsdoc:`Spectral Radiance` of a beam of light as it passes through a medium.
It can be described be Beer's law:

.. math::
  \vec{I}(z) = \vec{I}(0) \exp(-\mathbf{K} z),

where :math:`\vec{I}(z)` is the :artsdoc:`Spectral Radiance` of the light at some distance
:math:`z`, and :math:`\mathbf{K}` is the :artsdoc:`Propagation Matrix` of the medium.

.. include:: guide.theory.absorption.propagation_matrix.rst
.. include:: guide.theory.absorption.spectral_radiance.rst
.. include:: guide.theory.absorption.line-by-line.rst
