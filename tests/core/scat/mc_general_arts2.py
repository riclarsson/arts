"""ARTS2 passive-MC I/Q acceptance cases using the original ARO ice data."""

from pathlib import Path

import numpy as np
import pyarts3 as pyarts
from scipy.interpolate import RegularGridInterpolator


A = pyarts.arts
DATA = Path(__file__).with_name("data") / "mc_general_arts2"


def load(name):
    return np.asarray(pyarts.xml.load(str(DATA / name)))


p = load("p_grid.xml")
lat = load("lat_grid.xml")
lon = load("lon_grid.xml")
# Reconstruct AtmFieldsCalcExpand1D from source profiles, not its serialized
# 3D outputs. ARTS2 interpolates these quantities linearly in log pressure.
raw_p = A.GriddedField3.fromxml("planets/Earth/afgl/tropical/p.xml")
raw_t = A.GriddedField3.fromxml("planets/Earth/afgl/tropical/t.xml")
np.testing.assert_array_equal(raw_p.grids[0], raw_t.grids[0])
log_p = np.log(np.asarray(raw_p.data).ravel()[::-1])
z = np.interp(np.log(p), log_p, np.asarray(raw_p.grids[0])[::-1])
t_profile = np.interp(np.log(p), log_p, np.asarray(raw_t.data).ravel()[::-1])
shape = (len(p), len(lat), len(lon))
t = np.broadcast_to(t_profile[:, None, None], shape)
vmr = []
for species in ("O2", "N2", "H2O"):
    raw = A.GriddedField3.fromxml(f"planets/Earth/afgl/tropical/{species}.xml")
    np.testing.assert_array_equal(raw_p.grids[0], raw.grids[0])
    profile = np.interp(
        np.log(p),
        log_p,
        np.asarray(raw.data).ravel()[::-1],
    )
    vmr.append(np.broadcast_to(profile[:, None, None], shape))

# Original cloudbox indices: [13, 19, 23, 39, 23, 40], inclusive.
# Interpolate the raw particle population in (log pressure, latitude, longitude).
# Original pnd_field_raw.xml: a constant population inside a zero-padded box.
pnd_pressure = [
    110000.0,
    21618.7922264,
    21617.7922264,
    20390.6782501,
    19233.2202635,
    18141.4643086,
    17111.6808705,
    17110.6808705,
    1e-5,
]
pnd_latitude = [-90.0, -1.9, -1.8, -1.08, -0.36, 0.36, 1.08, 1.8, 1.9, 90.0]
pnd_longitude = [
    -180.0,
    -1.9,
    -1.8,
    -1.2,
    -0.6,
    -2.22044604925e-16,
    0.6,
    1.2,
    1.8,
    1.9,
    180.0,
]
raw_pnd = np.zeros((9, 10, 11))
raw_pnd[2:7, 2:8, 2:9] = 34712.0922774  # particles / m³
pnd_interpolate = RegularGridInterpolator(
    (np.log(pnd_pressure), pnd_latitude, pnd_longitude), raw_pnd
)
cloud_points = np.stack(
    np.meshgrid(np.log(p[13:20]), lat[23:40], lon[23:41], indexing="ij"), axis=-1
)
pnd = pnd_interpolate(cloud_points)

ws = pyarts.Workspace()
# Keep the original MC pressure-grid top, not the full AFGL 120-km extent.
ws.atm_field.top_of_atmosphere = float(z[-1])
grids = [z, lat, lon]
names = ["altitude", "latitude", "longitude"]


def field(data, field_grids=grids):
    return A.GeodeticField3(data=data, grids=field_grids, grid_names=names)


ws.atm_field["t"] = field(t)
ws.atm_field["p"] = field(np.broadcast_to(p[:, None, None], t.shape))
for species, data in zip(("O2", "N2", "H2O"), vmr, strict=True):
    ws.atm_field[species] = field(data)

number_density = A.ScatteringSpeciesProperty(
    "arts2_ice", A.ParticulateProperty.NumberDensity
)
pnd_full = np.zeros_like(t)
pnd_full[13:20, 23:40, 23:41] = pnd
ws.atm_field[number_density] = field(pnd_full)

legacy = pyarts.xml.load(str(DATA / "azi-random_f229-231T214-225r100NP-1ar1_5ice.xml"))
meta = A.ScatteringMetaData()
# This is a one-particle population; size metadata only identifies its sole bin.
meta.mass = 1.0
meta.diameter_volume_equ = 1.0
meta.diameter_max = 1.0
habit = A.ParticleHabit.from_legacy_aro([legacy], [meta])
psd = A.MonodispersePSD(number_density, 0.0, 400.0)
ws.scat_species = [A.ScatteringHabit(habit, psd, 1.0, 1.0)]

ws.abs_speciesSet(species=["O2-PWR98", "N2-SelfContStandardType", "H2O-PWR98"])
# These predefined models have no fitted catalog payload; register their model
# names locally so the test remains independent of tmp/ and ARTS_DATA_PATH.
ws.abs_predef_dataReadSpeciesSplitCatalog(
    basename=str(DATA) + "/", name_missing=1, ignore_missing=1
)
ws.spectral_propmat_agendaAuto()
ws.max_stepsize = 3e3
ws.ray_path_observer_agendaSetGeometric(
    max_step_option="step", add_crossings=1, remove_non_atm=1
)
ws.spectral_rad_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_rad_surface_agendaSet(option="Blackbody")
ws.surf_fieldPlanet(option="Earth")
ws.surf_field[A.SurfaceKey("h")] = 500.0
ws.surf_field[A.SurfaceKey("t")] = float(np.interp(500.0, z, t[:, 0, 0]))

common = dict(
    frequency=230e9,
    sensor_pos=[95000.1, 7.61968838781, 0.0],
    sensor_los=[99.7841941981, 180.0],
    mc_seed=1729,
    mc_min_iter=128,
    mc_max_iter=128,
    mc_max_scatorder=30,
)


def check(reference):
    ws.MCGeneral(**common)
    rj = (
        np.asarray(ws.mc_spectral_rad)
        * A.constants.c**2
        / (2.0 * A.constants.k * common["frequency"] ** 2)
    )
    error = (
        np.asarray(ws.mc_error)
        * A.constants.c**2
        / (2.0 * A.constants.k * common["frequency"] ** 2)
    )
    accepted = np.abs(rj[:2] - reference) <= 4.0 * error[:2]
    assert np.all(accepted), (
        f"ARTS2 I/Q={np.asarray(reference)} K; "
        f"ARTS3 I/Q={rj[:2]} +/- {error[:2]} K; "
        f"ARTS3 U/V={rj[2:]} +/- {error[2:]} K"
    )


ws.mc_antenna.set_pencil_beam()
check((198.7, 7.9))

ws.mc_antenna.set_gaussian_fwhm(0.1137, 0.239)
check((198.6, 7.6))
