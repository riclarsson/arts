import pyarts
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy as copy


alt = 420e3
aas = [0]
zas = [180]
lats = [67]
lons = [20]

noise = 0.1
NFREQ = 1001

N = len(lats)

assert len(aas) == N and len(zas) == N and len(lats) == N and len(lons) == N

ws = pyarts.workspace.Workspace()

# %% Sampled frequency range

line_f0 = 118750348044.712
ws.frequency_grid = np.linspace(-50e6, 50e6, 1001) + line_f0

# %% Species and line absorption

ws.absorption_speciesSet(species=["O2-66"])
ws.ReadCatalogData()
ws.absorption_bandsSelectFrequency(fmin=40e9, fmax=120e9, by_line=1)
ws.absorption_bandsSetZeeman(species="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet

ws.surface_fieldSetPlanetEllipsoid(option="Earth")
ws.surface_field[pyarts.arts.SurfaceKey("t")] = 295.0
ws.atmospheric_fieldRead(
    toa=120e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atmospheric_fieldIGRF(time="2000-03-11 14:39:37")

# %% Checks and settings

ws.spectral_radiance_unit = "Tb"
ws.spectral_radiance_observer_agendaSet(option="EmissionUnits")
ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")
ws.ray_path_observer_agendaSet(option="Geometric")

# %% Extract single point magnetic field

uk = pyarts.arts.AtmKey.mag_u
vk = pyarts.arts.AtmKey.mag_v
wk = pyarts.arts.AtmKey.mag_w

ws.atmospheric_fieldRegrid(parameter=uk, alt = [80e3], lat = lats, lon = lons)
ws.atmospheric_fieldRegrid(parameter=vk, alt = [80e3], lat = lats, lon = lons)
ws.atmospheric_fieldRegrid(parameter=wk, alt = [80e3], lat = lats, lon = lons)

u = copy(ws.atmospheric_field[uk].data.data)
v = copy(ws.atmospheric_field[vk].data.data)
w = copy(ws.atmospheric_field[wk].data.data)

# %% Jacobian

ws.RetrievalInit()
ws.RetrievalAddMagneticField(component="U", matrix=np.diag(np.ones((len(u))) * 1e-6**2))
ws.RetrievalAddMagneticField(component="V", matrix=np.diag(np.ones((len(v))) * 1e-6**2))
ws.RetrievalAddMagneticField(component="W", matrix=np.diag(np.ones((len(w))) * 1e-6**2))
ws.RetrievalFinalizeDiagonal()


@pyarts.workspace.arts_agenda(ws=ws, fix=True)
def inversion_iterate_agenda(ws):
    ws.UpdateModelStates()
    ws.measurement_vectorFromSensor()
    ws.measurement_vector_fittedFromMeasurement()

# %% Core calculations

for i in range(N):
  pos = [119e3, lats[i], lons[i]]
  los = [zas[i], aas[i]]
  ws.ray_pathGeometric(pos=pos, los=los, max_step=1000.0)
  ws.spectral_radianceClearskyEmission()
  ws.spectral_radianceApplyUnitFromSpectralRadiance()
  ws.measurement_sensorSimple(pos=pos, los=los)
  ws.measurement_vectorFromSensor()

  ws.measurement_vector_fitted = []
  ws.model_state_vector = []
  ws.measurement_jacobian = [[]]

  ws.atmospheric_field[uk].data += 10e-6
  ws.model_state_vector_aprioriFromData()

  ws.measurement_vector_error_covariance_matrixConstant(value=noise**2)
  ws.measurement_vector += np.random.normal(0, noise, NFREQ)

  #ws.OEM(method="gn")