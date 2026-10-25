"""Exact periodic-model correlation and variance checks; not an RH proof."""
from fractions import Fraction
from math import gcd, prod


def prime_factors(n):
    factors = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            factors.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        factors.append(n)
    return factors


def periodic_values(q):
    phi = sum(gcd(n, q) == 1 for n in range(q))
    rho = Fraction(q, phi)
    nu = [rho if gcd(n, q) == 1 else Fraction(0) for n in range(q)]
    return phi, rho, nu


def singular_series(q, h):
    result = Fraction(1)
    for p in prime_factors(q):
        result *= Fraction(p * (p - 1 if h % p == 0 else p - 2), (p - 1) ** 2)
    return result


def direct_autocorrelation(nu, h):
    q = len(nu)
    return sum((nu[n] * nu[(n + h) % q] for n in range(q)), Fraction(0)) / q


def check(q, h):
    assert q >= 2 and h >= 1 and prod(prime_factors(q)) == q
    phi, rho, nu = periodic_values(q)
    for shift in range(min(q, 12)):
        assert direct_autocorrelation(nu, shift) == singular_series(q, shift)

    formula = h * (rho - 1)
    for shift in range(1, h):
        formula += 2 * (h - shift) * (singular_series(q, shift) - 1)

    # Scale each centered value by phi to use integer prefix sums.
    centered_scaled = [q if gcd(n, q) == 1 else 0 for n in range(q)]
    centered_scaled = [value - phi for value in centered_scaled]
    repeated = [centered_scaled[i % q] for i in range(q + h)]
    prefix = [0]
    for value in repeated:
        prefix.append(prefix[-1] + value)
    window_numerators = [prefix[x + h] - prefix[x] for x in range(q)]
    direct_variance = Fraction(sum(value * value for value in window_numerators), q * phi * phi)
    assert formula == direct_variance

    omega = len(prime_factors(q))
    bound = rho * rho * (4**omega)
    assert direct_variance <= bound
    return {
        "q": q,
        "H": h,
        "phi": phi,
        "rho": rho,
        "V": direct_variance,
        "bound": bound,
    }


if __name__ == "__main__":
    rows = [check(q, h) for q in (2, 6, 30, 210) for h in (4, 16, 64)]
    print(f"{len(rows)} exact Fraction checks passed: autocorrelation, variance identity, bound.")
    representative = next(row for row in rows if row["q"] == 210 and row["H"] == 64)
    print("q=210,H=64 representative:")
    print("  phi=", representative["phi"], "rho=", representative["rho"])
    print("  V=", representative["V"], "=", float(representative["V"]))
    print("  rho^2*4^omega=", representative["bound"], "=", float(representative["bound"]))
    print("Exact finite model only; no RH statement is proved.")
