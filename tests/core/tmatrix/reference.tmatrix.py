"""Original ARTS2 reference values, exercised through ARTS3's Python API."""

from pathlib import Path
import re
from concurrent.futures import ThreadPoolExecutor
from threading import Barrier

import numpy as np
import pyarts3 as pa

TM = pa.arts.tmatrix
QUAD = TM.extended_precision()
REFERENCES = Path(__file__).resolve().parents[3] / "3rdparty" / "tmatrix"


def number(text):
    return float(text.replace("D", "e"))


def fixed(radius=10.0, wavelength=2 * np.pi, **kw):
    kw.setdefault("radius_ratio", 0.1 if QUAD else 1.0)
    args = dict(
        radius=radius,
        wavelength=wavelength,
        aspect_ratio=0.5,
        refractive_real=1.5,
        refractive_imag=0.02,
        theta_incident=56.0,
        theta_scattered=65.0,
        phi_incident=114.0,
        phi_scattered=128.0,
        alpha=145.0,
        beta=52.0,
    )
    args.update(kw)
    return TM.fixed(**args)


def random(radius=1.0, quadrature=5):
    return TM.random(
        radius=radius,
        wavelength=0.5,
        aspect_ratio=2.0,
        refractive_real=1.33 if QUAD else 1.53,
        refractive_imag=0.001 if QUAD else 0.008,
        radius_ratio=0.5,
        distribution=3,
        b=0.1,
        gamma=0.5,
        size_quadrature=quadrature,
        lower_radius_ratio=0.89031,
        upper_radius_ratio=1.56538,
    )


def small_random():
    # COMMON-block ownership/thread checks do not need the expensive reference
    # size distribution. Keep a nonspherical, absorbing particle in these checks.
    return TM.random(radius=0.1, wavelength=1.0, aspect_ratio=2.0,
                     refractive_real=1.5, refractive_imag=0.01)


def test_fixed():
    text = (
        REFERENCES / ("tmatrix_amplq.ref" if QUAD else "tmatrix_ampld.ref")
    ).read_text()
    result = fixed()
    assert isinstance(result.phase, pa.arts.Muelmat)
    if not QUAD:
        assert result.order == int(re.search(r"NMAX=\s*(\d+)", text)[1])
    values = re.findall(r"S\d\d=\s*([^ ]+) \+ i\*\s*([^\n]+)", text)
    actual = np.asarray(result.amplitude).ravel()
    for value, (real, imag) in zip(actual, values, strict=True):
        # One unit in the last printed decimal place, independently for each
        # real/imaginary component. Never recompute or replace reference values.
        for got, reference in [(value.real, real), (value.imag, imag)]:
            exponent = int(reference.split("D")[1])
            decimals = len(reference.split("D")[0].split(".")[1])
            np.testing.assert_allclose(
                got, number(reference), rtol=0, atol=10.0 ** (exponent - decimals)
            )
    phase = np.loadtxt(text.split("PHASE MATRIX\n")[1].splitlines()[:4])
    np.testing.assert_allclose(result.phase, phase, rtol=0, atol=1e-4)
    # The returned matrices own their data; a new COMMON-block calculation
    # must not alter a previously returned result.
    saved = np.array(result.amplitude)
    scaled = fixed(radius=10e-6, wavelength=2e-6 * np.pi)
    np.testing.assert_allclose(
        np.asarray(scaled.amplitude) / 1e-6, saved, rtol=2e-6, atol=1e-8
    )
    np.testing.assert_allclose(
        np.asarray(scaled.phase) / 1e-12, result.phase, rtol=2e-6, atol=1e-6
    )
    small_random()
    np.testing.assert_array_equal(result.amplitude, saved)
    # The old quad wrapper ignored the shape argument and always used spheroids.
    cylinder = fixed(shape=-2)
    assert not np.allclose(cylinder.amplitude, saved)
    np.testing.assert_array_equal(result.amplitude, saved)
    return saved


def test_random():
    text = (REFERENCES / ("tmatrix_tmq.ref" if QUAD else "tmatrix_tmd.ref")).read_text()
    summary = re.findall(
        r"CEXT=\s*(\S+)\s+CSCA=\s*(\S+)\s+W=\s*(\S+)\s+<COS>=\s*(\S+)", text
    )
    tables = re.findall(r"      <\s+F11[^\n]+\n((?:[^\n]+\n){19})", text)
    assert len(tables) == 2
    # The original single call scans radii 1 and .5, with 7 and 4 size
    # quadrature points. The ARTS3 interface computes one distribution/call.
    for i, (radius, quadrature) in enumerate([(1.0, 5), (0.5, 2)]):
        result = random(radius, quadrature)
        expected = np.loadtxt(tables[i].splitlines())[:, 1:]
        assert isinstance(result.phase, pa.arts.MuelmatVector)
        phase = np.asarray(result.phase)
        actual = phase[:, [0, 1, 2, 3, 0, 2], [0, 1, 2, 3, 1, 3]]
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-4)
        np.testing.assert_array_equal(phase[:, 1, 0], phase[:, 0, 1])
        np.testing.assert_array_equal(phase[:, 3, 2], -phase[:, 2, 3])
        for name, reference in zip(
            ["extinction", "scattering", "albedo", "asymmetry"], summary[i], strict=True
        ):
            np.testing.assert_allclose(
                getattr(result, name),
                number(reference),
                rtol=0,
                atol=1e-5 if name in ["extinction", "scattering"] else 1e-6,
            )
        np.testing.assert_allclose(
            result.scattering / result.extinction, result.albedo, rtol=1e-14
        )
        np.testing.assert_allclose(result.effective_radius, radius, rtol=0, atol=5e-5)
        np.testing.assert_allclose(result.effective_variance, 0.1, rtol=0, atol=5e-5)
    # The power-law root solver used an absolute lower bound of 1e-5.
    # The ARTS3 boundary must permit the same particle specified in metres.
    # Reuse the second reference result, including its four size nodes.
    base = result
    si = TM.random(
        radius=0.5e-6,
        wavelength=0.5e-6,
        aspect_ratio=2.0,
        refractive_real=1.33 if QUAD else 1.53,
        refractive_imag=0.001 if QUAD else 0.008,
        radius_ratio=0.5,
        distribution=3,
        b=0.1,
        gamma=0.5,
        size_quadrature=2,
    )
    np.testing.assert_allclose(si.phase, base.phase, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(si.scattering / 1e-12, base.scattering, rtol=1e-12)
    # The third CEXT/CSCA footer in the ARTS2 reference is the old wrapper bug:
    # last-particle cross sections paired with distribution-averaged albedo.
    # Compare the correct, unchanged distribution summaries above instead.
    assert len(summary) == (2 if QUAD else 3)


def test_errors_and_threads():
    for kwargs in [
        {"radius": -1},
        {"accuracy": float("nan")},
        {"shape": 0},
        {"beta": 181},
        {"phi_incident": -1},
    ]:
        try:
            fixed(**kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError(f"invalid input accepted: {kwargs}")
    try:
        fixed(radius=1e4)
    except RuntimeError as e:
        assert "CONVERGENCE" in str(e)
    else:
        raise AssertionError("array limit did not produce a Python exception")

    # Exercise real concurrent entry with cheap particles. Repeating the large
    # reference cases here only serializes expensive solves on the solver mutex.
    expected = np.array(fixed(radius=0.1, wavelength=1.0).amplitude)
    expected_random = np.array(small_random().phase)
    barrier = Barrier(4)

    def evaluate(i):
        barrier.wait(timeout=30)
        random_result = small_random() if i % 2 else None
        result = fixed(radius=0.1, wavelength=1.0)
        if random_result is not None:
            np.testing.assert_array_equal(random_result.phase, expected_random)
        return np.array(result.amplitude)

    with ThreadPoolExecutor(max_workers=4) as pool:
        for actual in pool.map(evaluate, range(8)):
            np.testing.assert_array_equal(actual, expected)


def test_fixed_batch():
    args = (100e-6, 299792458 / 230e9, 1.5, 1.78, .003)
    geometries = np.array([[0., 90., 0., 30., 0., 0.],
                           [56., 65., 114., 128., 145., 52.],
                           [90., 90., 0., 0., 0., 0.]])
    batch = TM.fixed_batch(*args, geometries)
    assert len(batch) == len(geometries)
    saved = [np.array(r.phase) for r in batch]
    # New solver calls must not overwrite earlier rows or retained batches.
    for i, geometry in enumerate(geometries):
        scalar = TM.fixed(*args, *geometry)
        np.testing.assert_array_equal(batch[i].amplitude, scalar.amplitude)
        np.testing.assert_array_equal(batch[i].phase, scalar.phase)
        np.testing.assert_array_equal(batch[i].phase, saved[i])
    assert TM.fixed_batch(*args, np.empty((0, 6))) == []
    for invalid in [np.zeros((2, 5)),
                    np.array([[0., 0., 0., 0., 0., 181.]]),
                    np.array([[0., 0., float('nan'), 0., 0., 0.]])]:
        try:
            TM.fixed_batch(*args, invalid)
        except ValueError:
            pass
        else:
            raise AssertionError('invalid batch geometry accepted')


assert TM.available()
test_fixed()
test_random()
test_errors_and_threads()

test_fixed_batch()
