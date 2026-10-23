"""Numerical and end-to-end tests of the CPU reference, using independent mpmath."""
import csv
import io
import math
import os
from pathlib import Path
import subprocess
import sys
import unittest

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
SCANNER = Path(os.environ.get("RH_SCANNER", ROOT / "build/rhscan-cpu")).resolve()
sys.path.insert(0, str(ROOT / "tools"))
from verify import verify

mp.mp.dps = 50


def scan(start, end, step, executable=SCANNER):
    result = subprocess.run(
        [str(executable), str(start), str(end), str(step)],
        capture_output=True, text=True, check=True,
        env={**os.environ, "OMP_NUM_THREADS": "2"},
    )
    return list(csv.DictReader(io.StringIO(result.stdout)))


def values(points):
    result = subprocess.run(
        [str(ROOT / "build/probe"), *map(str, points)],
        capture_output=True, text=True, check=True,
    )
    return list(map(float, result.stdout.splitlines()))


class ScannerTests(unittest.TestCase):
    def test_independent_values(self):
        points = [200, 201.5, 500, 1000, 10000, 1000000, 1e8, 1e9]
        for t, actual in zip(points, values(points)):
            with self.subTest(t=t):
                # Empirical regression tolerance, not an analytic error bound.
                tolerance = 0.2 * t ** -0.75 + 2e-5
                self.assertLess(abs(actual - float(mp.siegelz(t))), tolerance)

    def test_removable_singularities(self):
        points = [2 * math.pi * (n + p + delta) ** 2
                  for n in [6, 40] for p in [0.25, 0.75]
                  for delta in [-1e-10, 0, 1e-10]]
        for t, actual in zip(points, values(points)):
            self.assertTrue(math.isfinite(actual))
            self.assertLess(abs(actual - float(mp.siegelz(t))), 0.2 * t ** -0.75)

    def test_sum_cutoff_transitions(self):
        points = [2 * math.pi * n * n + delta
                  for n in [6, 20] for delta in [-1e-8, 0, 1e-8]]
        for t, actual in zip(points, values(points)):
            self.assertLess(abs(actual - float(mp.siegelz(t))), 0.2 * t ** -0.75)

    def test_known_region(self):
        rows = scan(200, 240, 0.05)
        # Independent count used only to validate this particular small region.
        self.assertEqual(len(rows), int(mp.nzeros(240) - mp.nzeros(200)))
        for row in rows:
            a, b = mp.mpf(row["t_left"]), mp.mpf(row["t_right"])
            self.assertLess(mp.siegelz(a) * mp.siegelz(b), 0)
        root = verify(rows[0]["t_left"], rows[0]["t_right"])
        self.assertLess(abs(mp.siegelz(root)), mp.mpf("1e-32"))
        if SCANNER.name == "rhscan-cuda":
            cpu_rows = scan(200, 240, 0.05, ROOT / "build/rhscan-cpu")
            self.assertEqual(len(rows), len(cpu_rows))
            for gpu, cpu in zip(rows, cpu_rows):
                self.assertEqual(gpu["t_left"], cpu["t_left"])
                self.assertEqual(gpu["t_right"], cpu["t_right"])
                for key in ["z_approx_left", "z_approx_right"]:
                    self.assertAlmostEqual(float(gpu[key]), float(cpu[key]), delta=1e-10)

    def test_cross_batch_brackets(self):
        # Locate an approximate crossing, then put it across the batch boundary.
        root = float(mp.zetazero(80).imag)
        step = 0.000001
        short = scan(root - 0.02, root + 0.02, step)
        self.assertEqual(len(short), 1)
        left = float(short[0]["t_left"])
        start = left - 65535 * step
        rows = scan(start, root + 0.02, step)
        self.assertEqual(len(rows), 1)
        self.assertAlmostEqual(float(rows[0]["t_left"]), left, delta=step / 10)

    def test_invalid_inputs(self):
        for args in [(199, 240, .05), (240, 200, .05), (200, 240, 0),
                     (200, 240, "nan"), (200, "inf", .05), (200, 240, "1e-300"),
                     (200, 240, "junk"), (200, 1e10, 1)]:
            result = subprocess.run([str(SCANNER), *map(str, args)],
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0, args)

    def test_final_partial_interval(self):
        rows = scan(200, 201.3, 0.8)
        self.assertEqual(len(rows), 1)
        self.assertEqual(float(rows[0]["t_right"]), 201.3)
        self.assertEqual(float(rows[0]["t_left"]), 200.8)

    def test_verifier_rejects_false_candidate(self):
        with self.assertRaises(ValueError):
            verify("200", "200.01")


if __name__ == "__main__":
    subprocess.run(["make", "cpu", "build/probe"], cwd=ROOT, check=True)
    unittest.main(verbosity=2)
