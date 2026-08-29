"""Seeded regression tests for the ARTS3-native MCRadar interface."""

import numpy as np
import pyarts3 as pyarts

A = pyarts.arts

ws = pyarts.Workspace()
ws.surf_fieldPlanet(option="Earth")
ws.surf_field[A.SurfaceKey("t")] = 280.0
ws.atm_fieldRead(
    toa=30e3,
    basename="planets/Earth/afgl/tropical/",
    missing_is_zero=1,
)
ws.abs_speciesSet(species=[])
ws.spectral_propmat_agendaAuto()
ws.ray_path_observer_agendaSetGeometric(
    max_step_option="step",
    add_crossings=1,
    remove_non_atm=1,
)

# A horizontally homogeneous, isotropic cloud.  Extinction is deliberately
# small enough that both the cloud top and bottom remain visible.
extinction = A.ScatteringSpeciesProperty(
    "cloud", A.ParticulateProperty.Extinction
)
ssa = A.ScatteringSpeciesProperty(
    "cloud", A.ParticulateProperty.SingleScatteringAlbedo
)
grids = ws.atm_field["t"].data.grids
altitude = np.asarray(grids[0])
cloud = (altitude >= 5e3) & (altitude <= 15e3)
ws.atm_field[extinction] = A.GriddedField3(
    data=np.where(cloud, 1e-4, 0.0)[:, None, None], grids=grids
)
ws.atm_field[ssa] = A.GriddedField3(
    data=np.where(cloud, 0.8, 0.0)[:, None, None], grids=grids
)
ws.scat_species = [A.HenyeyGreensteinScatterer(extinction, ssa, 0.0)]

ws.mc_antenna.set_gaussian_fwhm(0.6, 0.6)
kwargs = dict(
    frequency=94e9,
    sensor_pos=[20e3, 0.0, 0.0],
    sensor_los=[180.0, 0.0],
    mc_y_tx=[1.0, 1.0, 0.0, 0.0],
    range_bins=np.arange(0.0, 20.001e3, 1e3),
    mc_seed=1729,
    mc_max_iter=64,
    mc_max_scatorder=1,
    unit="1",
)

ws.MCRadar(**kwargs)
signal = np.array(ws.radar_signal)
error = np.array(ws.radar_error)

assert signal.shape == (20, 4)
assert error.shape == signal.shape
assert np.all(np.isfinite(signal))
assert np.all(np.isfinite(error))
assert np.all(error >= 0.0)
assert np.any(signal[:, 0] > 0.0)

# The same seed must give bitwise-identical Monte Carlo output.
ws.MCRadar(**kwargs)
assert np.array_equal(signal, np.asarray(ws.radar_signal))
assert np.array_equal(error, np.asarray(ws.radar_error))

# Ze is a linear conversion of the native return, including its uncertainty.
ws.MCRadar(**(kwargs | {"unit": "Ze"}))
ze = np.array(ws.radar_signal)
ze_error = np.array(ws.radar_error)
mask = signal != 0.0
factor = ze[mask] / signal[mask]
assert np.allclose(factor, factor[0], rtol=1e-13)
error_mask = error != 0.0
assert np.allclose(ze_error[error_mask] / error[error_mask], factor[0])

# Invalid higher scattering orders are rejected explicitly until the collision
# sampler is added, rather than silently returning a single-scattering result.
try:
    ws.MCRadar(**(kwargs | {"mc_max_scatorder": 2}))
except RuntimeError:
    pass
else:
    raise AssertionError("MCRadar accepted an unsupported scattering order")
