.. _Sec Empirical microwave surface emissivity:

Empirical microwave surface emissivity
#######################################

ARTS provides the TESSEM2 and TELSEM2 empirical microwave emissivity models
through the ``Tessem`` and ``Telsem`` options of
:meth:`~pyarts3.Workspace.spectral_surf_refl_agendaSet`.  Load their model data
with :class:`~pyarts3.arts.TessemNN` or
:class:`~pyarts3.arts.TelsemAtlas`, or with the corresponding workspace ASCII
reader methods.

TESSEM2
*******

TESSEM reads the following local values from ``surf_field``:

- ``"t"``: surface skin temperature [K];
- ``"wind speed"``: 10-m wind speed [m/s]; and
- ``"salinity"``: salinity [kg/kg].

Each value can use any normal :class:`~pyarts3.arts.SurfaceData`
representation, including a constant, a latitude/longitude field, or a
function.  Surface Jacobian targets for all three inputs are evaluated by
finite differences using the target perturbation.

TELSEM2
*******

TELSEM obtains its spatial variation from its monthly atlas.  It uses
``surf_field`` only for the local surface geometry and evaluates the atlas for
the ray point and local surface normal.

Surface-boundary behavior
*************************

Both agenda options return polarized specular reflectance.  Thermal emission
is supplied by the closed-surface agenda, which applies Kirchhoff consistency.

FASTEM remains an intended future capability, but it cannot be added as
another option of this reflectance agenda.  Its atmospheric correction and
potentially independent emissivity and reflectivity require a different
surface-boundary contract, which still needs to be designed.
