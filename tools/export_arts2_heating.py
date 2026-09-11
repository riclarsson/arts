#!/usr/bin/env python3
"""Export fixed RT inputs with ARTS2's Python extension (not pyarts3).

Usage: python export_arts2_heating.py /path/to/arts2 /path/to/output
Set PYTHONPATH to the ARTS2 build/python and ARTS_DATA_PATH to its
controlfiles and XML data directories. Output extinction.py provides the literal
EXTINCTION and REFERENCE_HEATING blocks for arts2.py. Reference results are
read from the existing ARTS2 reference files, never regenerated.
"""

import argparse
import os
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
from pyarts.workspace import Workspace
from pyarts.workspace.workspace import Include

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("arts2", type=Path)
parser.add_argument("output", type=Path)
args = parser.parse_args()
root, output = args.arts2.resolve(), args.output.resolve()
output.mkdir(parents=True, exist_ok=True)
literals = ["import numpy as np", "", "# fmt: off", "EXTINCTION = {"]

references = ["REFERENCE_HEATING = {"]

for name, folder in [
    ("TestHeatingRates", "controlfiles"),
    ("Test_HeatingRate", "controlfiles-python"),
]:
    source = root / folder / "artscomponents/heatingrates"
    os.chdir(source)
    if name == "TestHeatingRates":
        ws = Workspace()
        Include(ws, name + ".arts")
    else:
        # Do not allow the original script's automatic LUT recalculation path:
        # it writes into the old source tree and would change the inputs.
        script = (source / (name + ".py")).read_text()
        script = script.replace(
            "recalc = True", 'raise RuntimeError("Original LUT required")'
        )
        namespace = {}
        exec(compile(script, str(source / (name + ".py")), "exec"), namespace)
        ws = namespace["ws"]
    ws.propmat_clearsky_fieldCalc()
    p, z, t = (
        np.array(getattr(ws, key).value).ravel()
        for key in ("p_grid", "z_field", "t_field")
    )
    data = {
        "extinction": np.array(ws.propmat_clearsky_field.value).sum(axis=0).squeeze()
    }
    if name == "TestHeatingRates":
        vmr = np.array(ws.vmr_field.value).squeeze()
        for i, za in enumerate(np.array(ws.za_grid.value)[:3]):
            maximum_dz = 1e4 * abs(np.cos(np.deg2rad(za)))
            nz = np.concatenate(
                [
                    np.linspace(lo, hi, int(np.ceil((hi - lo) / maximum_dz)) + 1)[:-1]
                    for lo, hi in zip(z[:-1], z[1:])
                ]
                + [[z[-1]]]
            )
            ws.p_grid = np.exp(np.interp(nz, z, np.log(p)))
            ws.z_field = nz[:, None, None]
            ws.t_field = np.interp(nz, z, t)[:, None, None]
            ws.vmr_field = np.array([np.interp(nz, z, v) for v in vmr])[
                :, :, None, None
            ]
            ws.propmat_clearsky_fieldCalc()
            data.update(
                {
                    f"k{i}": np.array(ws.propmat_clearsky_field.value)
                    .sum(axis=0)
                    .squeeze(),
                }
            )
    literals.append(f'    "{name}": {{')
    keys = ("k0", "k1", "k2") if name == "TestHeatingRates" else ("extinction",)
    for key in keys:
        literals.append(f'        "{key}": np.array([')
        for row in data[key]:
            literals.append(
                "            [" + ", ".join(repr(float(v)) for v in row) + "],"
            )
        literals.append("        ]),")
    literals.append("    },")
    for reference in sorted(source.glob(name + ".heating*REFERENCE.xml")):
        key = reference.name.removesuffix("REFERENCE.xml")
        tokens = ET.parse(reference).find(".//Tensor3").text.split()
        references.append(f'    "{key}": np.array([')
        references.extend("        " + value + "," for value in tokens)
        references.append("    ]),")

references.append("}")
literals.extend(["}", "", *references, "# fmt: on", ""])
(output / "extinction.py").write_text("\n".join(literals))
