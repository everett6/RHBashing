"""Finite numerical check of the Type I Fourier grouping and large-sieve bound."""
import cmath
from fractions import Fraction
import math

import numpy as np

from check_fourier_reduction import mangoldt
from check_global_vaughan import mobius


def floor_cuberoot(n):
    q = int(round(n ** (1.0 / 3.0)))
    while (q + 1) ** 3 <= n:
        q += 1
    while q**3 > n:
        q -= 1
    return q


def make_coefficients(u, v):
    mu = mobius(v + 1)
    lam = mangoldt(u + 1)
    qmax = u * v
    c = np.zeros(qmax + 1)
    for k in range(1, v + 1):
        for d in range(1, u + 1):
            c[k * d] += mu[k] * lam[d]
    return c


def direct_values(n, h, c):
    qmax = len(c) - 1
    out = []
    for x in range(n, 2 * n):
        value = Fraction(0, 1)
        for q in range(1, qmax + 1):
            if c[q]:
                value += Fraction.from_float(float(c[q])) * (
                    Fraction((x - h) % q, q) - Fraction(x % q, q)
                )
        out.append(float(value))
    return np.asarray(out)


def fourier_coefficients(c, h):
    coeffs = {}
    for q in range(2, len(c)):
        cq = float(c[q])
        if cq == 0:
            continue
        for j in range(1, q):
            freq = Fraction(j, q)
            z = cmath.exp(-2j * math.pi * j / q)
            shift = cmath.exp(-2j * math.pi * j * h / q)
            a = -cq / (q * (1.0 - z)) * (shift - 1.0)
            coeffs[freq] = coeffs.get(freq, 0j) + a
    return {freq: a for freq, a in coeffs.items() if abs(a) > 1e-14}


def check(n):
    h = n // 2
    u = v = floor_cuberoot(n)
    c = make_coefficients(u, v)
    assert all(abs(c[q]) <= math.log(q) + 1e-12 for q in range(1, len(c)))
    direct = direct_values(n, h, c)
    coeffs = fourier_coefficients(c, h)

    reconstructed = np.asarray([
        sum(a * cmath.exp(2j * math.pi * float(freq) * x)
            for freq, a in coeffs.items())
        for x in range(n, 2 * n)
    ])
    direct_error = float(np.max(np.abs(direct - reconstructed.real)))
    assert np.max(np.abs(reconstructed.imag)) < 1e-9
    assert direct_error < 1e-9

    energy = float(direct @ direct)
    coefficient_energy = float(sum(abs(a) ** 2 for a in coeffs.values()))
    qmax = u * v
    harmonic = np.concatenate(([0.0], np.cumsum(1.0 / np.arange(1, qmax + 1))))
    spectral_bound = math.log(2 * qmax)**2 * sum(
        harmonic[qmax // d]**2 for d in range(2, qmax + 1))
    assert coefficient_energy <= spectral_bound + 1e-9
    bound = (n + (u * v) ** 2) * coefficient_energy
    assert energy <= bound + 1e-9 * max(1.0, bound)
    return u, v, u * v, direct_error, energy, coefficient_energy, bound


if __name__ == "__main__":
    print("N,U,V,Q,max_reconstruction_error,E3_energy,coefficient_energy,large_sieve_bound")
    for n in (64, 216, 512, 1000):
        print(n, *check(n), sep=",")
