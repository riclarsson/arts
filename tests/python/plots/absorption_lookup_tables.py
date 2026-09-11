import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pyarts


def make_table(scale, water=False):
    table = pyarts.arts.AbsorptionLookupTable()
    table.f_grid = pyarts.arts.AscendingGrid([1e9, 2e9, 3e9, 4e9])
    table.log_p_grid = pyarts.arts.DescendingGrid(np.log([1e5, 5e4, 1e4]))
    table.t_pert = pyarts.arts.AscendingGrid([-10.0, 10.0])
    table.w_pert = pyarts.arts.AscendingGrid([0.5, 1.5] if water else [])
    table.t_atmref = [290.0, 270.0, 240.0]
    table.water_atmref = [0.01, 0.005, 0.001]

    shape = (2, 2 if water else 1, 3, 4)
    table.xsec = scale * (1.0 + np.arange(np.prod(shape))).reshape(shape)
    return table


tables = pyarts.arts.AbsorptionLookupTables()
tables[pyarts.arts.SpeciesEnum.H2O] = make_table(1e-27, water=True)
tables[pyarts.arts.SpeciesEnum.O2] = make_table(2e-28)

# Generic dispatch, species selection, pressure selection, and perturbation
# dimensions are covered without creating any output files.
fig, axes = pyarts.plot(
    tables,
    mode="cross_section",
    pressures=[1e5, 1e4],
    temperature_perturbation=1,
    water_perturbation=1,
    cols=2,
)
assert len(fig.axes) == 2
assert all(len(axis.lines) == 2 for axis in fig.axes)
assert {axis.get_title() for axis in fig.axes} == {"H2O", "O2"}
plt.close(fig)

fig, axis = pyarts.plots.AbsorptionLookupTables.plot(
    tables,
    mode="cross_section",
    species="O2",
)
assert len(fig.axes) == 1
assert axis.get_title() == "O2"
plt.close(fig)

# H2O uses its stored reference VMR. O2 exercises the explicit profile input.
fig, axis = pyarts.plot(
    tables,
    mode="opacity",
    vmr_profiles={"O2": 0.21},
    altitude=[0.0, 5e3, 16e3],
)
assert axis.get_yscale() == "log"
assert len(axis.lines) == 4  # H2O, O2, total, and opacity=1.
assert {line.get_label() for line in axis.lines[:3]} == {"H2O", "O2", "Total"}
assert all(np.all(np.isfinite(line.get_ydata())) for line in axis.lines)
plt.close(fig)
