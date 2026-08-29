"""Monte Carlo integration with sensor weighting and time-tagged metadata."""

import numpy as np
import pyarts3 as pyarts

A = pyarts.arts
ws = pyarts.Workspace()

ws.surf_fieldPlanet(option="Earth")
ws.surf_field[A.SurfaceKey("t")] = 280.0
ws.atm_fieldRead(toa=30e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1)
ws.abs_speciesSet(species=[])
ws.scat_species = []
ws.spectral_propmat_agendaAuto()
ws.ray_path_observer_agendaSetGeometric(
    max_step_option="step", add_crossings=1, remove_non_atm=1
)
ws.spectral_rad_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_rad_surface_agendaSet(option="Blackbody")


@pyarts.arts_agenda(ws=ws, fix=False)
def spectral_rad_observer_agenda(ws):
    ws.spectral_radMonteCarlo(mc_seed=9182, mc_min_iter=8, mc_max_iter=8)


freq = A.AscendingGrid([220e9, 230e9])
pos = [20e3, 0.0, 0.0]
ws.measurement_sensorInit()
ws.measurement_sensorAddSimple(freq_grid=freq, pos=pos, los=[180.0, 0.0])
ws.measurement_sensor_meta[-1].time = A.Time("2001-01-01 00:00:00")
ws.measurement_sensorAddSimple(freq_grid=freq, pos=pos, los=[0.0, 0.0])
ws.measurement_sensor_meta[-1].time = A.Time("2001-01-01 00:01:00")

ws.measurement_vecFromSensor()

expected_surface = [A.physics.planck(f, 280.0) for f in freq]
expected_space = [A.physics.planck(f, 2.725) for f in freq]
assert np.allclose(np.asarray(ws.measurement_vec), expected_surface + expected_space)
assert len(ws.measurement_sensor_meta) == 2
assert all(meta.count == 2 for meta in ws.measurement_sensor_meta)
assert ws.measurement_sensor_meta[0].time != ws.measurement_sensor_meta[1].time

ws.measurement_sensor_metaFromMeasurementVec()
assert np.allclose(np.asarray(ws.measurement_sensor_meta[0].data.data), expected_surface)
assert np.allclose(np.asarray(ws.measurement_sensor_meta[1].data.data), expected_space)
