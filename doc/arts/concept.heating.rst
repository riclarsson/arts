Radiative heating rates
=======================

A radiative heating rate is a gas temperature tendency.  Positive values mean
warming.  Heating depends on the divergence of energy flux, not on
angle-integrated radiance alone: the flux integral contains the projected-area
factor given by the direction cosine.

Let :math:`F=F_\uparrow-F_\downarrow` be the upward net flux, integrated over
frequency.  In a plane-parallel atmosphere, with altitude :math:`z` increasing
upwards, density :math:`\rho`, and mass-specific heat capacity :math:`c_p`,

.. math::

   \frac{dT}{dt} = -\frac{1}{\rho c_p}\frac{dF}{dz}.

Under hydrostatic balance, :math:`dp/dz=-\rho g`, this becomes

.. math::

   \frac{dT}{dt} = \frac{g}{c_p}\frac{dF}{dp}.

These vertical relations omit horizontal flux divergence and spherical-area
divergence.  Finite differences on pressure levels approximate the derivative;
the pressure stencil and boundary treatment affect that approximation.

Optical-depth derivatives
-------------------------

Optical depth increases downwards, so that
:math:`d\tau_\nu/dz=-k_{\mathrm{ext},\nu}`.  For the spectral net flux,

.. math::

   D_\nu = \frac{dF_\nu}{d\tau_\nu},\qquad
   \frac{dT}{dt} = \frac{1}{\rho c_p}
       \int k_{\mathrm{ext},\nu}D_\nu\,d\nu.

Total extinction is required here: the absorption fraction is already
represented in the flux derivative.  Frequency-dependent extinction must be
applied before frequency integration, and it must correspond to the
optical-depth coordinate used for the derivative.  A local optical-depth
derivative and a finite difference of flux on pressure levels need not give
identical numerical results.

See :doc:`concept.disort` for the discrete-ordinate flux expressions and
:doc:`user.heating` for units, sampling locations, and recipe examples.
