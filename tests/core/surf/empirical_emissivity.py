from pathlib import Path

import numpy as np
import pyarts3 as pyarts


data = Path(__file__).parent
horizontal = pyarts.arts.TessemNN.from_ascii(str(data / "tessem_sav_net_H.txt"))
vertical = pyarts.arts.TessemNN.from_ascii(str(data / "tessem_sav_net_V.txt"))
model_input = [1.0e10, 0.0, 0.0, 273.14999, 0.003]

# Reference values from ARTS 2 TestTessem.
eh = np.asarray(horizontal(model_input))
ev = np.asarray(vertical(model_input))
np.testing.assert_allclose(eh, [0.395911], atol=1e-6, rtol=0)
np.testing.assert_allclose(ev, [0.374513], atol=1e-6, rtol=0)

# Exercise the ARTS 3 surface-field integration, including agenda setup.
ws = pyarts.Workspace()
ws.freq_grid = [1.0e10]
ws.surf_fieldEarth()
ws.surf_field["t"] = model_input[3]
ws.surf_field["wind speed"] = model_input[2]
ws.surf_field["salinity"] = model_input[4]
ws.ray_point.pos = [0.0, 0.0, 0.0]
ws.ray_point.los = [180.0, 0.0]
ws.tessem_neth = horizontal
ws.tessem_netv = vertical
ws.spectral_surf_refl_agendaSet(option="Tessem")
ws.spectral_surf_refl_agendaExecute()

r = np.asarray(ws.spectral_surf_refl[0])
rv = 1.0 - ev[0]
rh = 1.0 - eh[0]
np.testing.assert_allclose(r[0, 0], 0.5 * (rv + rh), atol=1e-12)
np.testing.assert_allclose(r[0, 1], 0.5 * (rv - rh), atol=1e-12)

# A minimal atlas checks TELSEM ASCII parsing, grid inversion, and the 53-degree
# base-channel behavior without requiring the externally distributed atlas.
atlas = pyarts.arts.TelsemAtlas.from_ascii(str(data / "telsem_atlas.txt"), month=1)
lat, lon = atlas.coordinates(330034)
assert atlas.cell_number(lat, lon) == 330034
np.testing.assert_allclose(
    atlas.emissivity(lat, lon, 53.0, 19.35e9), [0.90, 0.80], atol=1e-12
)

ws.telsem_atlas = atlas
ws.ray_point.pos = [0.0, lat, lon]
ws.ray_point.los = [127.0, 0.0]
ws.freq_grid = [19.35e9]
ws.spectral_surf_refl_agendaSet(option="Telsem")
ws.spectral_surf_refl_agendaExecute()
