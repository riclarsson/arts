import pyarts3 as pyarts

ws = pyarts.Workspace()

@pyarts.arts_agenda
def create_agenda(ws: pyarts.Workspace):
  ws.abs_bandsFromModelState()

@pyarts.in_parallel(ws=ws)
def create_agenda_parallel(ws: pyarts.Workspace):
  ws.abs_speciesSet(species=["H2O", "CO2"])
  ws.atm_pointInit()

assert ws.abs_species == ["H2O", "CO2"]
