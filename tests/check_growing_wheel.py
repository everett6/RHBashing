"""Numerical audit of the growing-wheel decomposition; not an RH proof."""
from fractions import Fraction
from math import factorial, gcd, log

import numpy as np

from check_fourier_reduction import mangoldt


PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31)


def wheel_for(n):
    q = 1
    primes = []
    for p in PRIMES:
        if q * p > n:
            break
        q *= p
        primes.append(p)
    phi = sum(gcd(k, q) == 1 for k in range(q))
    return q, len(primes), Fraction(q, phi)


def windows(a, n, h):
    prefix = np.concatenate(([0.0], np.cumsum(a)))
    x = np.arange(n, 2 * n)
    return prefix[x + 1] - prefix[x - h + 1]


def check(n):
    h = n // 2
    top = 2 * n - 1
    q, w, rho_exact = wheel_for(n)
    rho = float(rho_exact)
    lam = mangoldt(top + 1)
    coprime = np.zeros(top + 1, dtype=bool)
    coprime[1:] = [gcd(k, q) == 1 for k in range(1, top + 1)]

    c = np.zeros(top + 1)
    c[1:] = (lam[1:] - rho) * coprime[1:]
    s = windows(lam, n, h) - h
    centered = windows(c, n, h)
    coprime_count = windows(coprime.astype(float), n, h)
    delta = rho * coprime_count - h
    noncoprime_lam = lam * (~coprime)
    noncoprime_lam[0] = 0.0
    p_term = windows(noncoprime_lam, n, h)

    decomposition_error = s - centered - delta - p_term
    assert np.max(np.abs(decomposition_error)) < 1e-8

    delta_bound = rho * (2**w)
    assert np.max(np.abs(delta)) <= delta_bound + 1e-9
    p_bound = log(q)
    p_bound_pass = bool(np.max(p_term) <= p_bound + 1e-9)
    assert p_bound_pass
    error = s - centered
    error_bound = delta_bound + p_bound
    error_bound_pass = bool(np.max(np.abs(error)) <= error_bound + 1e-9)
    assert error_bound_pass

    psi = float(np.sum(lam[1:top + 1]))
    coefficient_energy = float(c[1:top + 1] @ c[1:top + 1])
    coefficient_bound = (
        2 * log(2 * n) * psi
        + 2 * rho * top
        + 2 * rho * rho * (2**w)
    )
    assert coefficient_energy <= coefficient_bound + 1e-9 * max(1.0, coefficient_bound)

    assert rho_exact <= w + 1
    assert factorial(w + 1) <= q
    rho_wheel_bound = rho_exact * rho_exact * (2**w)
    assert rho_wheel_bound <= 6 * q

    return {
        "N": n,
        "q": q,
        "w": w,
        "rho": rho_exact,
        "max_error": float(np.max(np.abs(error))),
        "error_bound": error_bound,
        "error_energy": float(error @ error),
        "coefficient_energy": coefficient_energy,
        "coefficient_bound": coefficient_bound,
        "max_delta": float(np.max(np.abs(delta))),
        "delta_bound": delta_bound,
        "max_P": float(np.max(p_term)),
        "P_bound": p_bound,
        "P_bound_pass": p_bound_pass,
        "error_bound_pass": error_bound_pass,
        "coefficient_bound_pass": coefficient_energy <= coefficient_bound + 1e-9 * max(1.0, coefficient_bound),
    }


if __name__ == "__main__":
    print("N,q,w,rho,maxerror,bound,errorenergy,coefficientenergy,coefficientbound,")
    print("maxdelta,deltabound,maxP,logq,P_bound_pass,error_bound_pass")
    rows = [check(n) for n in (16, 64, 256, 1024, 4096)]
    for row in rows:
        print(",".join(str(row[k]) for k in (
            "N", "q", "w", "rho", "max_error", "error_bound", "error_energy",
            "coefficient_energy", "coefficient_bound", "max_delta", "delta_bound",
            "max_P", "P_bound", "P_bound_pass", "error_bound_pass"
        )))
    print("Numerical checks only; failures of stated finite bounds are reported above.")
