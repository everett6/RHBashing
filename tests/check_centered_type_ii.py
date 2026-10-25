"""Numerical audit of a centered Type II window decomposition; not an RH test."""
import math
import numpy as np

from check_fourier_reduction import mangoldt
from check_global_vaughan import convolution, mobius


def window_sums(a, n):
    """Return sum_{r=0}^{H-1} a[x-r] for x=N,...,2N-1."""
    h = n // 2
    prefix = np.concatenate(([0.0], np.cumsum(a)))
    x = np.arange(n, 2 * n)
    return prefix[x + 1] - prefix[x - h + 1]


def check(n):
    h = n // 2
    u = v = math.isqrt(math.isqrt(n))
    limit = 2 * n

    lam = mangoldt(limit)
    mu = mobius(limit)
    one = np.zeros(limit)
    one[1:] = 1.0
    log = np.zeros(limit)
    log[1:] = np.log(np.arange(1, limit))

    # Verify the exact truncated Vaughan identity at every coefficient < 2N.
    small_mu = mu.copy()
    small_mu[v + 1:] = 0.0
    small_lam = lam.copy()
    small_lam[u + 1:] = 0.0
    exact_b = convolution(one, lam - small_lam)
    term1 = convolution(small_mu, log)
    term3 = convolution(convolution(small_mu, one), small_lam)
    term4 = convolution(mu - small_mu, exact_b)
    reconstruction = term1 + small_lam - term3 + term4
    vaughan_error = float(np.max(np.abs(lam - reconstruction)))
    assert np.allclose(lam, reconstruction, atol=1e-10, rtol=1e-10)

    # Requested Type II B = mu_{>V} * (1 * Lambda_{>U}).
    mu_gt_v = mu.copy()
    mu_gt_v[:v + 1] = 0.0
    lam_gt_u = lam.copy()
    lam_gt_u[:u + 1] = 0.0
    b_u = convolution(one, lam_gt_u)
    type_ii = convolution(mu_gt_v, b_u)

    a_const = float(sum(mu[k] / k for k in range(1, v + 1)))
    c_const = float(sum(mu[k] * math.log(k) / k for k in range(1, v + 1)))
    l_const = float(sum(lam[d] / d for d in range(1, u + 1)))

    x = np.arange(n, 2 * n, dtype=float)
    lower = x - h
    j = x * np.log(x) - x - lower * np.log(lower) + lower
    t_ii = window_sums(type_ii, n)
    r = t_ii + a_const * j - h * (c_const + a_const * l_const + 1.0)
    s = window_sums(lam, n) - h
    diff = s - r
    e1 = window_sums(term1, n) - (a_const * j - h * c_const)
    e3 = window_sums(term3, n) - h * a_const * l_const
    assert np.allclose(diff, e1 - e3, atol=1e-8, rtol=1e-9)

    psi_u = float(np.sum(lam[1:u + 1]))
    proposed_bound = 2.0 * v * math.log(2 * n) + v * psi_u
    max_error = float(np.max(np.abs(diff)))
    l2_error = float(np.linalg.norm(diff))
    assert max_error <= proposed_bound + 1e-9

    return {
        "N": n,
        "U": u,
        "V": v,
        "vaughan_error": vaughan_error,
        "max_error": max_error,
        "l2_error": l2_error,
        "bound": proposed_bound,
        "S_energy_over_N3": float(s @ s) / n**3,
        "R_energy_over_N3": float(r @ r) / n**3,
        "diff_energy_over_N3": float(diff @ diff) / n**3,
    }


if __name__ == "__main__":
    print("N,U,V,vaughan_error,max_error,l2_error,bound,S2_over_N3,R2_over_N3,diff2_over_N3")
    for n in (2, 16, 64, 256, 1024, 4096):
        row = check(n)
        print(",".join(str(row[k]) for k in (
            "N", "U", "V", "vaughan_error", "max_error", "l2_error", "bound",
            "S_energy_over_N3", "R_energy_over_N3", "diff_energy_over_N3"
        )))
