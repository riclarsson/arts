"""Hydrostatic pressure against analytic ideal-gas columns; no external data."""

import numpy as np
import pyarts3 as pa

A = pa.arts
G, R, P0, Z0 = 9.81, 287.0, 90000.0, 1200.0


def column(altitudes, temperatures, option, fixed=False):
    ws = pa.Workspace()
    ws.atm_field.top_of_atmosphere = float(altitudes[-1])
    ws.gravity_operator = A.NumericTernaryOperator(lambda z, lat, lon: G)
    ws.atm_field['t'] = A.GeodeticField3(
        data=np.asarray(temperatures)[:, None, None],
        grids=[altitudes, [0.0], [0.0]],
        grid_names=['Altitude', 'Latitude', 'Longitude'],
    )
    p0 = P0
    if fixed:
        p0 = A.GeodeticField2(data=[[P0]], grids=[[0.0], [0.0]],
                             grid_names=['Latitude', 'Longitude'])
    ws.atm_fieldHydrostaticPressure(
        alt_grid=altitudes, p0=p0, fixed_specific_gas_constant=R,
        fixed_atmospheric_temperature=float(temperatures[0]) if fixed else -1,
        hydrostatic_option=option,
    )
    return ws.atm_field['p']


# A nonzero reference altitude catches confusion between absolute altitude and
# height above p0. Include off-grid points and extrapolation at both ends.
z = np.linspace(Z0, Z0 + 4000, 41)
query = np.array([Z0 - 100, Z0, Z0 + 175, z[-1], z[-1] + 100])
for temperature in (100.0, 273.15, 1000.0):
    expected = P0 * np.exp(-G * (query - Z0) / (R * temperature))
    for fixed in (False, True):
        pressure = column(z, np.full(z.size, temperature), 'HypsometricEquation', fixed)
        actual = np.array([pressure(h, 0.0, 0.0) for h in query])
        np.testing.assert_allclose(actual, expected, rtol=2e-14)
        assert np.all(np.diff(actual) < 0)

# Positive temperature increasing with altitude (an inversion). Integrating
# dp/p = -g dz/(R T) gives the independent closed form below. Both available
# discretizations should approach it under vertical refinement.
lapse = 0.01
height = 4000.0
expected = P0 * (1 + lapse * height / 220.0)**(-G / (R * lapse))
for option in ('HypsometricEquation', 'HydrostaticEquation'):
    errors = []
    for layers in (40, 80, 160):
        altitudes = np.linspace(Z0, Z0 + height, layers + 1)
        temperature = 220.0 + lapse * (altitudes - Z0)
        pressure = column(altitudes, temperature, option)
        np.testing.assert_allclose(pressure(Z0, 0.0, 0.0), P0, rtol=1e-14)
        actual = pressure(altitudes[-1], 0.0, 0.0)
        assert 0 < actual < P0
        errors.append(abs(actual / expected - 1))
    assert errors[-1] < 0.002, (option, errors)
    assert 0.4 < errors[1] / errors[0] < 0.6, (option, errors)
    assert 0.4 < errors[2] / errors[1] < 0.6, (option, errors)
