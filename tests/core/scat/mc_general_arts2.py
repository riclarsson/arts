"""ARTS2 passive-MC I/Q acceptance cases using T-matrix-generated ARO ice data."""

from pathlib import Path

import numpy as np
import pyarts3 as pyarts
from scipy.interpolate import RegularGridInterpolator


A = pyarts.arts
DATA = Path(__file__).with_name("data") / "mc_general_arts2"


def load(name):
    return np.asarray(pyarts.xml.load(str(DATA / name)))


# Warren/Wiscombe/Gao REFICE, the material model used by historical arts_scat.
# Source: https://radiativetransfer.org/svn/rt/PyARTS/trunk/src/REFICE.f
# Evaluated with REAL inputs wavelength = c/f * 1e6 [micrometres], T [K],
# and COMPLEX (single precision) output, matching its original interface.
# Rows: frequency; columns: temperature. These are NOT fitted optical data.
REFRACTIVE_INDEX = np.array([
    [1.7800652980804443 + 0.00246472773142159j,
     1.7805060148239136 + 0.003070596605539322j],
    [1.7800800800323486 + 0.0024795890785753727j,
     1.7805144786834717 + 0.0030892009381204844j],
])
FREQUENCIES = [229e9, 231e9]
TEMPERATURES = [214.0, 225.0]
ANGLES = np.arange(0.0, 181.0, 10.0)
C = 299792458.0


def generate(legacy_numerics=False, quadrature=10):
    # Match the historical integration in theta and relative azimuth, both
    # 0..180 degrees, with sin(theta) and a factor 2 for the other azimuth half.
    if legacy_numerics and quadrature == 10:
        # Original arts_math.py printed its quadrature to eight decimal places.
        x = np.array([.14887434, .43339539, .67940957, .86506337, .97390653])
        w = np.array([.29552422, .26926672, .21908636, .14945135, .06667134])
        x, w = np.r_[x, -x], np.r_[w, w]
    else:
        x, w = np.polynomial.legendre.leggauss(quadrature)
    theta = 90 * (1 + x)
    weights = 2 * (np.pi / 2)**2 * w[:, None] * w[None, :] * np.sin(np.deg2rad(theta))[:, None]

    # Matrix rows: incident zenith, scattered zenith, incident azimuth,
    # scattered azimuth, alpha, beta. Oblate symmetry axis stays vertical.
    phase_geometries = [[i, s, 0, a, 0, 0] for s in ANGLES for a in ANGLES for i in ANGLES]
    forward_geometries = [[z, z, 0, 0, 0, 0] for z in ANGLES]
    integral_geometries = [[i, s, 0, a, 0, 0] for s in ANGLES for i in theta for a in theta]
    geometries = np.array(phase_geometries + forward_geometries + integral_geometries)
    nphase = len(phase_geometries)
    nz = len(ANGLES)
    phase = np.empty((2, 2, nz, nz, nz, 1, 16))
    extinction = np.empty((2, 2, nz, 1, 3))
    absorption = np.empty((2, 2, nz, 1, 2))
    orders = []
    for fi, frequency in enumerate(FREQUENCIES):
        wavelength = C / frequency
        radius = 100e-6
        unit = 1.0
        if legacy_numerics:
            # arts_scat stored wavelength in float32 micrometres before calling
            # the double-precision solver, then converted amplitudes back to SI.
            wavelength = float(np.float32(wavelength * 1e6))
            radius, unit = 100., 1e-6
        for ti, _ in enumerate(TEMPERATURES):
            m = REFRACTIVE_INDEX[fi, ti]
            results = A.tmatrix.fixed_batch(radius, wavelength, 1.5, m.real, m.imag,
                                             geometries, accuracy=.001, shape=-1)
            orders.append(int(results[0].order))
            phase[fi, ti] = (np.array([np.asarray(r.phase) for r in results[:nphase]]) * unit**2).reshape(nz, nz, nz, 1, 16)
            # Optical theorem: directional extinction, not the orientation-
            # averaged scalar cross section returned by the T-matrix solve.
            for zi, r in enumerate(results[nphase:nphase + nz]):
                s = np.asarray(r.amplitude) * unit
                lam = wavelength * unit
                k = A.Propmat([lam * (s[0, 0] + s[1, 1]).imag,
                               lam * (s[0, 0] - s[1, 1]).imag,
                               0., 0., 0., 0., lam * (s[1, 1] - s[0, 0]).real])
                extinction[fi, ti, zi, 0] = np.asarray(k)[[0, 1, 6]]
            # Integrate column 0 at fixed outgoing direction, as in arts_scat;
            # subtract scattering out of thermal equilibrium from K[:,0].
            z = np.array([np.asarray(r.phase) for r in results[nphase + nz:]]) * unit**2
            z = z.reshape(nz, quadrature, quadrature, 4, 4)
            scattering = np.einsum('ij,zijc->zc', weights, z[..., :, 0])
            for zi in range(nz):
                a = A.Stokvec([extinction[fi, ti, zi, 0, 0] - scattering[zi, 0],
                               extinction[fi, ti, zi, 0, 1] - scattering[zi, 1], 0., 0.])
                absorption[fi, ti, zi, 0] = np.asarray(a)[:2]
    ssd = A.SingleScatteringData()
    ssd.ptype = A.PType.PTYPE_AZIMUTH_RND
    ssd.description = (
        "ARTS3 T-matrix MC particle: equal-volume radius 100 um, oblate spheroid "
        "aspect ratio 1.5, alpha=beta=0, accuracy=0.001, original REFICE material. "
        f"Extended precision: {A.tmatrix.extended_precision()}; "
        f"legacy numerics: {legacy_numerics}; absorption quadrature: {quadrature}. "
        "Full angular grid evaluated directly."
    )
    ssd.f_grid, ssd.T_grid = FREQUENCIES, TEMPERATURES
    ssd.za_grid, ssd.aa_grid = ANGLES, ANGLES
    ssd.pha_mat_data, ssd.ext_mat_data, ssd.abs_vec_data = phase, extinction, absorption
    return ssd, orders



def main():
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

    legacy, _ = generate()
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


if __name__ == "__main__":
    main()
