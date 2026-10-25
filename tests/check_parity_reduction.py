"""Finite parity/window identities for dyadic N; not an RH proof."""
import math

import numpy as np

from check_fourier_reduction import mangoldt


def window(a, n, h):
    prefix = np.concatenate(([0.0], np.cumsum(a)))
    x = np.arange(n, 2 * n)
    return prefix[x + 1] - prefix[x - h + 1]


def overlap(n, h, left, right):
    lo, hi = min(left, right), max(left, right)
    return max(0, min(2 * n - 1, lo + h - 1) - max(n, hi) + 1)


def check(n, compare_odd_loop=False):
    h = n // 2
    top = 2 * n - 1
    lam = mangoldt(top + 1)
    a = np.zeros(top + 1)
    r = np.zeros(top + 1)
    c = np.zeros(top + 1)
    odd = np.arange(top + 1) % 2 == 1
    a[1:] = lam[1:] - 1.0
    r[1:] = lam[1:] - 2.0 * odd[1:]
    c[1:] = (lam[1:] - 2.0) * odd[1:]

    wa = window(a, n, h)
    wr = window(r, n, h)
    wc = window(c, n, h)
    assert np.allclose(wa, wr, atol=1e-10, rtol=1e-12)

    target = np.zeros(n)
    target[:h] = math.log(2.0)
    reduction = wa - wc
    reduction_error = reduction - target
    assert np.allclose(reduction_error, 0.0, atol=1e-10, rtol=1e-12)
    error2 = float(reduction @ reduction)
    expected_error2 = h * math.log(2.0) ** 2
    assert abs(error2 - expected_error2) < 1e-9

    energy = float(wc @ wc)
    diagonal = 0.0
    for m in range(1, top + 1):
        diagonal += c[m] ** 2 * overlap(n, h, m, m)
    even_offdiag = energy - diagonal

    # c is supported only on odd indices, so all odd shifts vanish exactly.
    odd_shift_abs_c = 0.0
    for shift in range(1, h, 2):
        for left in range(1, top - shift + 1):
            if c[left] and c[left + shift]:
                odd_shift_abs_c += 2 * abs(c[left] * c[left + shift]) * overlap(
                    n, h, left, left + shift
                )
    assert odd_shift_abs_c < 1e-9

    log2 = math.log(2.0)
    exceptional_abs = 0.0
    exceptional_signed = 0.0
    for m in range(1, top + 1, 2):
        k = overlap(n, h, n, m)
        exceptional_abs += 2 * log2 * abs(lam[m] - 2.0) * k
        exceptional_signed += 2 * log2 * (lam[m] - 2.0) * k
    odd_abs_bound = 2 * h * log2 * sum(abs(lam[m] - 2.0) for m in range(1, top + 1, 2))
    assert exceptional_abs <= odd_abs_bound + 1e-9 * max(1.0, odd_abs_bound)

    loop_signed = None
    loop_abs_pairs = None
    if compare_odd_loop:
        loop_signed = 0.0
        loop_abs_pairs = 0.0
        for shift in range(1, h, 2):
            for left in range(1, top - shift + 1):
                k = overlap(n, h, left, left + shift)
                product = r[left] * r[left + shift] * k
                loop_signed += 2 * product
                loop_abs_pairs += 2 * abs(r[left] * r[left + shift]) * k
        assert abs(loop_signed - exceptional_signed) < 1e-8
        assert abs(loop_abs_pairs - exceptional_abs) < 1e-8

    return {
        "N": n,
        "H": h,
        "E": energy,
        "diagonal": diagonal,
        "even_offdiag": even_offdiag,
        "odd_abs_r": exceptional_abs,
        "odd_abs_bound": odd_abs_bound,
        "odd_signed_r": exceptional_signed,
        "error2": error2,
        "expected_error2": expected_error2,
        "odd_loop_signed": loop_signed,
        "odd_loop_abs_pairs": loop_abs_pairs,
    }


if __name__ == "__main__":
    print("N,H,E,diagonal,even_offdiag,oddabs_r,odd_abs_bound,odd_signed_r,error2,expected_error2")
    for n in (16, 64, 256, 1024, 4096):
        row = check(n, compare_odd_loop=(n in (16, 64)))
        print(",".join(str(row[k]) for k in (
            "N", "H", "E", "diagonal", "even_offdiag", "odd_abs_r",
            "odd_abs_bound", "odd_signed_r", "error2", "expected_error2"
        )))
        if row["odd_loop_signed"] is not None:
            print(f"  odd-h loop check: signed={row['odd_loop_signed']}, "
                  f"absolute pair sum={row['odd_loop_abs_pairs']}")
    print("Finite numerical checks only; no RH statement is proved.")
