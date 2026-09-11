"""Code-defined predefined instruments through a workspace method."""

import numpy as np
import pyarts3 as pyarts


def stokes_sum(obsel):
    return np.asarray(obsel.weight_matrix.dense).sum(axis=(0, 1))


ws = pyarts.Workspace()
position = [817e3, 0.0, 0.0]
line_of_sight = [135.0, 0.0]

for style, channel_count in (
    (pyarts.arts.PredefinedSensor.AMSUA, 15),
    (pyarts.arts.PredefinedSensor.AMSUB, 5),
    (pyarts.arts.PredefinedSensor.ICI, 11),
):
    ws.measurement_sensorFromPredefined(
        pos=position,
        los=line_of_sight,
        sensor=style,
    )
    assert len(ws.measurement_sensor) == channel_count
    assert len(ws.measurement_sensor_meta) == 1
    assert ws.measurement_sensor_meta[0].count == channel_count
    assert all(np.isclose(stokes_sum(obsel)[0], 1.0) for obsel in ws.measurement_sensor)

for style, channel_count in (("HIRS", 19), ("MVIRI", 3), ("SEVIRI", 12)):
    ws.measurement_sensorFromPredefined(
        pos=position,
        los=line_of_sight,
        sensor=style,
    )
    assert len(ws.measurement_sensor) == channel_count
    assert all(np.isclose(stokes_sum(obsel)[0], 1.0) for obsel in ws.measurement_sensor)

# The copied NOAA-19 AVHRR fast response only defines its thermal channel;
# its two solar channels are retained as zero-response rows.
ws.measurement_sensorFromPredefined(
    pos=position,
    los=line_of_sight,
    sensor=pyarts.arts.PredefinedSensor.AVHRR,
)
assert len(ws.measurement_sensor) == 3
assert np.allclose([stokes_sum(obsel)[0] for obsel in ws.measurement_sensor], [0.0, 0.0, 1.0])

#Channel numbers are one - based and order - preserving.AMSU - A channel 5 is
#horizontal and channel 1 is vertical in the instrument definition.
ws.measurement_sensorFromPredefined(
    pos=position,
    los=line_of_sight,
    sensor="AMSU-A",
    n=1,
    channels=[5, 1],
)
assert len(ws.measurement_sensor) == 2
assert np.allclose(stokes_sum(ws.measurement_sensor[0]), [1.0, -1.0, 0.0, 0.0])
assert np.allclose(stokes_sum(ws.measurement_sensor[1]), [1.0, 1.0, 0.0, 0.0])

#Frequency sampling is controlled directly by an integer.
ws.measurement_sensorFromPredefined(
    pos=position,
    los=line_of_sight,
    sensor="ICI",
    n=81,
    channels=[1],
)
assert len(ws.measurement_sensor) == 1
assert len(ws.measurement_sensor[0].f_grid) > 100
