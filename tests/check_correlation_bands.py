"""Finite correlation-band audit for the centered Type II window; not a proof."""
import math

import numpy as np

from check_fourier_reduction import mangoldt
from check_global_vaughan import convolution, mobius


def floor_cuberoot(n):
    q = int(round(n ** (1.0 / 3.0)))
    while (q + 1) ** 3 <= n:
        q += 1
    while q**3 > n:
        q -= 1
    return q


def centered_coefficients(n):
    h = n // 2
    u = v = floor_cuberoot(n)
    limit = 2 * n
    lam = mangoldt(limit)
    mu = mobius(limit)
    one = np.zeros(limit)
    one[1:] = 1.0
    log = np.zeros(limit)
    log[1:] = np.log(np.arange(1, limit))

    mu_gt_v = mu.copy()
    mu_gt_v[:v + 1] = 0.0
    lam_gt_u = lam.copy()
    lam_gt_u[:u + 1] = 0.0
    b_u = convolution(one, lam_gt_u)
    type_ii = convolution(mu_gt_v, b_u)

    a = float(sum(mu[k] / k for k in range(1, v + 1)))
    c = float(sum(mu[k] * math.log(k) / k for k in range(1, v + 1)))
    ell = float(sum(lam[d] / d for d in range(1, u + 1)))
    center = 1.0 + c + a * ell

    centered = type_ii - center + a * log
    centered[0] = 0.0
    return centered, h, u, v


def window_values(c, n, h):
    prefix = np.concatenate(([0.0], np.cumsum(c)))
    x = np.arange(n, 2 * n)
    return prefix[x + 1] - prefix[x - h + 1]


def pair_overlap(n, h, left, shift):
    """Number of target windows containing both left and left+shift."""
    return max(0, min(2 * n - 1, left + h - 1) - max(n, left + shift) + 1)


def check(n, do_bands=True):
    coeff, h, u, v = centered_coefficients(n)
    w = window_values(coeff, n, h)
    energy = float(w @ w)

    diagonal = 0.0
    for left in range(1, 2 * n):
        k = max(0, min(2 * n - 1, left + h - 1) - max(n, left) + 1)
        diagonal += coeff[left] ** 2 * k

    shifts = {}
    abs_shifts = {}
    for shift in range(1, h):
        signed = 0.0
        absolute = 0.0
        for left in range(1, 2 * n - shift):
            k = pair_overlap(n, h, left, shift)
            if k:
                product = coeff[left] * coeff[left + shift]
                signed += product * k
                absolute += abs(product) * k
        shifts[shift] = 2.0 * signed
        abs_shifts[shift] = 2.0 * absolute

    off_total = float(sum(shifts.values()))
    identity_error = abs(energy - diagonal - off_total)
    assert identity_error <= 1e-8 * max(1.0, energy)

    if not do_bands:
        return off_total

    m = math.isqrt(n)
    short_end = min(h - 1, math.floor(math.log(n)))
    last_start = max(1, h - m)
    bands = {
        "short": range(1, short_end + 1),
        "middle": range(short_end + 1, last_start),
        "last": range(last_start, h),
    }
    band_values = {}
    band_absolute = {}
    for name, band in bands.items():
        indices = [shift for shift in band if shift in shifts]
        band_values[name] = float(sum(shifts[shift] for shift in indices))
        band_absolute[name] = float(sum(abs_shifts[shift] for shift in indices))

    positive_count = sum(value > 0 for value in shifts.values())
    total_abs = float(sum(abs_shifts.values()))
    last_abs = band_absolute["last"]
    c_energy = float(coeff[1:2 * n] @ coeff[1:2 * n])
    endpoint_bound = m * (m + 1) * c_energy
    assert last_abs <= endpoint_bound + 1e-8 * max(1.0, endpoint_bound)

    return {
        "N": n,
        "U": u,
        "V": v,
        "H": h,
        "E": energy,
        "D_diagonal": float(diagonal),
        "O": off_total,
        "identity_error": identity_error,
        "positive_shift_count": positive_count,
        "shift_count": len(shifts),
        "O_abs": total_abs,
        "bands": band_values,
        "bands_abs": band_absolute,
        "M": m,
        "endpoint_abs": last_abs,
        "endpoint_bound": endpoint_bound,
    }


if __name__ == "__main__":
    print("N,U,V,H,E,D,O,identity_error,positive_shifts,total_shifts,|O|sum,")
    print("short,short_abs,middle,middle_abs,last,last_abs,M,endpoint_abs,endpoint_bound")
    for n in (16, 64, 216, 512, 1024, 4096):
        row = check(n)
        fields = [
            row[k] for k in ("N", "U", "V", "H", "E", "D_diagonal", "O",
                             "identity_error", "positive_shift_count", "shift_count", "O_abs")
        ]
        for band in ("short", "middle", "last"):
            fields.extend((row["bands"][band], row["bands_abs"][band]))
        fields.extend((row["M"], row["endpoint_abs"], row["endpoint_bound"]))
        print(",".join(str(x) for x in fields))

    positives = [n for n in range(8, 257, 2) if check(n, do_bands=False) > 0]
    print("even N in 8..256 with positive total O:", positives)
    print("Finite numerical checks only; no analytic correlation estimate is proved.")
