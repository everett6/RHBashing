"""Finite checks of the stated prime-energy bound and recurrence; not an RH proof.

These checks do not validate the source or hypotheses of an external theorem,
prove an asymptotic estimate, or establish any zero hypothesis.
"""
import math
from fractions import Fraction

import numpy as np

from check_fourier_reduction import mangoldt


GOLDEN = (1.0 + math.sqrt(5.0)) / 2.0


def exact_p(n):
    h = n // 2
    direct = sum((2 * x - h) ** 2 for x in range(n, 2 * n))
    formula = Fraction(79 * n**3, 12) - 5 * n**2 + Fraction(2 * n, 3)
    assert formula.denominator == 1
    assert direct == formula.numerator
    return direct


def amplitude_bound(n):
    log2n = math.log(2 * n)
    first = 9.39 * log2n**1.515 * math.exp(-0.8274 * math.sqrt(math.log(n / 2)))
    if n < 64:
        return first
    t = n / 2
    log_t = math.log(t)
    f = log_t**0.6 / math.log(log_t)**0.2
    second = 0.026 * log2n**1.801 * math.exp(-0.1853 * f)
    return min(first, second)


def prime_error_values(limit):
    lam = mangoldt(limit + 1)
    psi = np.cumsum(lam)
    n = np.arange(limit + 1, dtype=float)
    return psi - n


def variance(e, n):
    return float(e[n:2 * n] @ e[n:2 * n])


def difference_energy(e, n):
    h = n // 2
    d = e[n:2 * n] - e[n - h:2 * n - h]
    return float(d @ d)


def check_recurrence(e, n, label):
    assert n % 2 == 0
    left = math.sqrt(variance(e, n))
    right = math.sqrt(2.0 * variance(e, n // 2)) + GOLDEN * math.sqrt(difference_energy(e, n))
    tolerance = 1e-10 * max(1.0, left, right)
    assert left <= right + tolerance, (label, n, left, right)
    return left, right


def check_prime_bound(n):
    e = prime_error_values(2 * n - 1)
    f = difference_energy(e, n)
    p = exact_p(n)
    a = amplitude_bound(n)
    explicit = a * a * p
    assert f <= explicit + 1e-10 * max(1.0, explicit)
    return f, a, p, explicit


def check_power_mode():
    beta = 0.75
    nodes, weights = np.polynomial.legendre.leggauss(64)
    t_nodes = 1.5 + 0.5 * nodes
    g_nodes = (t_nodes**beta - (t_nodes - 0.5)**beta) ** 2
    integral = 0.5 * float(weights @ g_nodes)
    g1 = (1.0**beta - (1.0 - 0.5)**beta) ** 2
    g2 = (2.0**beta - (2.0 - 0.5)**beta) ** 2
    riemann_error_bound = g1 - g2

    # Check the recurrence and scaling identity for this growing toy,
    # not the prime bound. Its Mellin transform has a pole at beta>1/2.
    max_n = 512
    toy = np.zeros(2 * max_n + 1)
    toy[1:] = np.arange(1, 2 * max_n + 1, dtype=float) ** beta
    rows = []
    for n in (16, 32, 64, 128, 256, 512):
        check_recurrence(toy, n, "power toy")
        scaled_f = difference_energy(toy, n) / n ** (2 * beta + 1)
        assert scaled_f >= integral - 1e-11
        assert scaled_f <= integral + riemann_error_bound / n + 1e-11
        rows.append((n, scaled_f, integral, integral + riemann_error_bound / n))
    return beta, rows


if __name__ == "__main__":
    print("Exact P_N identity for dyadic N=8,...,4096:")
    for n in (2**k for k in range(3, 13)):
        p = exact_p(n)
        print(f"  N={n}: P_N={p}")

    print("Prime F(N) versus the stated explicit bound (floating-point check only):")
    for n in (8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096):
        f, a, p, bound = check_prime_bound(n)
        print(f"  N={n}: F={f:.12g}, A_N={a:.12g}, P_N={p}, A_N^2 P_N={bound:.12g}")

    print("Recurrence check for prime E(n)=psi(n)-n:")
    prime_e = prime_error_values(8192)
    for n in (16, 32, 64, 128, 256, 512, 1024, 2048, 4096):
        left, right = check_recurrence(prime_e, n, "prime E")
        print(f"  N={n}: sqrt(V)={left:.12g}, RHS={right:.12g}")

    print("Recurrence check for deterministic random real E(n):")
    rng = np.random.default_rng(20261025)
    random_e = np.zeros(8193)
    random_e[1:] = 0.1 * np.arange(1, 8193) + rng.standard_normal(8192)
    for n in (16, 30, 64, 100, 256, 512, 1024, 2048, 4096):
        left, right = check_recurrence(random_e, n, "random E")
        print(f"  N={n}: sqrt(V)={left:.12g}, RHS={right:.12g}")

    beta, toy_rows = check_power_mode()
    print(f"Power-mode scale check: beta={beta}, 64-node Gauss-Legendre integral for g(t) on [1,2].")
    for n, scaled_f, lower, upper in toy_rows:
        print(f"  N={n}: F/N^(2beta+1)={scaled_f:.12g}, integral={lower:.12g}, "
              f"upper Riemann bound={upper:.12g}")
    print("Quadrature and finite checks are numerical only; no asymptotic or zero claim is proved.")
