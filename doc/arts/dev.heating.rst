Heating-rate reference tests
============================

The active tests in ``tests/core/heating`` embed all three original heating
reference arrays and retain their absolute tolerances. Their README records the
frozen optical inputs, ray subdivisions, source choices, signs, and gravity
model. Historical pressure-stencil and boundary quirks are isolated in a test
adapter; they are not defaults of the public recipe. The regression reruns
ARTS3 radiative transfer and does not regenerate old absorption lookup tables.

The public interface is described in :doc:`user.heating`; the physical
relations are in :doc:`concept.heating`.
