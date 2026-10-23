#!/usr/bin/env python3
"""Independent numerical bracket checks; not a rigorous zero count or RH test."""
import argparse
import csv
import sys
import mpmath as mp


def verify(left, right):
    a, b = mp.mpf(left), mp.mpf(right)
    if not mp.isfinite(a) or not mp.isfinite(b) or not 200 <= a < b <= 10**9:
        raise ValueError("invalid bracket")
    fa, fb = mp.siegelz(a), mp.siegelz(b)
    if fa == 0:
        return a
    if fb == 0:
        return b
    if mp.sign(fa) == mp.sign(fb):
        raise ValueError("high-precision endpoints have the same sign")
    for _ in range(160):
        mid = (a + b) / 2
        fm = mp.siegelz(mid)
        if fm == 0:
            return mid
        if mp.sign(fm) == mp.sign(fa):
            a, fa = mid, fm
        else:
            b = mid
        if b - a < mp.mpf("1e-35"):
            return (a + b) / 2
    raise ValueError("bisection did not converge")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidates")
    args = parser.parse_args()
    mp.mp.dps = 50
    failed = checked = 0
    writer = csv.writer(sys.stdout)
    writer.writerow(["t_left", "t_right", "numerical_root", "status"])
    with open(args.candidates, newline="") as stream:
        reader = csv.DictReader(stream)
        if not {"t_left", "t_right"}.issubset(reader.fieldnames or []):
            parser.error("missing candidate CSV columns")
        for row in reader:
            left, right = row["t_left"], row["t_right"]
            try:
                root = verify(left, right)
                writer.writerow([left, right, mp.nstr(root, 40), "numerically_checked"])
                checked += 1
            except ValueError as error:
                writer.writerow([left, right, "", str(error)])
                failed += 1
    print(f"Numerically checked: {checked}; rejected: {failed}. No completeness claim.", file=sys.stderr)
    return int(failed > 0)


if __name__ == "__main__":
    sys.exit(main())
