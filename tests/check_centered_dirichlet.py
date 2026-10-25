"""High-precision finite-cutoff identities; numerical checks, not an RH proof."""
try:
    import mpmath as mp
except ImportError as exc:
    raise SystemExit("This check requires mpmath; install it in the active Python environment.") from exc


def mobius_table(limit):
    mu = [0] * (limit + 1)
    mu[1] = 1
    for n in range(1, limit + 1):
        for multiple in range(2 * n, limit + 1, n):
            mu[multiple] -= mu[n]
    return mu


def mangoldt_table(limit):
    lam = [mp.mpf(0)] * (limit + 1)
    for p in range(2, limit + 1):
        if all(p % q for q in range(2, int(p**0.5) + 1)):
            value = mp.log(p)
            power = p
            while power <= limit:
                lam[power] = value
                power *= p
    return lam


def make_cutoffs(u, v):
    limit = max(u, v)
    mu = mobius_table(v)
    lam = mangoldt_table(u)

    def m(s):
        return mp.fsum(mp.mpf(mu[k]) / mp.power(k, s) for k in range(1, v + 1))

    def ell(s):
        return mp.fsum(lam[d] / mp.power(d, s) for d in range(1, u + 1))

    a = mp.fsum(mp.mpf(mu[k]) / k for k in range(1, v + 1))
    c = mp.fsum(mp.mpf(mu[k]) * mp.log(k) / k for k in range(1, v + 1))
    l1 = ell(mp.mpf(1))
    d = 1 + c + a * l1
    return m, ell, a, c, l1, d


def values(s, cutoff):
    m, ell, a, c, l1, d = cutoff
    z = mp.zeta(s)
    zp = mp.diff(mp.zeta, s)
    zeta_log_derivative = -zp / z
    bcal = (1 - m(s) * z) * (zeta_log_derivative - ell(s))
    expanded = zeta_log_derivative + m(s) * zp + (m(s) * ell(s)) * z - ell(s)
    ccal = (
        zeta_log_derivative
        + (m(s) - a) * zp
        + (m(s) * ell(s) - d) * z
        - ell(s)
    )
    centered_relation = bcal - a * zp - d * z
    return bcal, expanded, ccal, centered_relation


def assert_close(x, y, tol=None):
    if tol is None:
        tol = mp.mpf("1e-60")
    scale = max(mp.mpf(1), abs(x), abs(y))
    assert abs(x - y) <= tol * scale, (x, y, abs(x - y))


def check_pair(u, v):
    cutoff = make_cutoffs(u, v)
    for s in (mp.mpf(2), mp.mpc(2, 1)):
        bcal, expanded, ccal, centered_relation = values(s, cutoff)
        assert_close(bcal, expanded)
        assert_close(ccal, centered_relation)

    one_eps = [mp.mpf("1e-3"), mp.mpf("1e-4"), mp.mpf("1e-5"), mp.mpf("1e-6")]
    near_one = [values(1 + eps, cutoff)[2] for eps in one_eps]
    one_steps = [abs(near_one[i + 1] - near_one[i]) for i in range(len(near_one) - 1)]
    assert all(one_steps[i + 1] < one_steps[i] for i in range(len(one_steps) - 1)), one_steps

    rho = mp.zetazero(1)
    zero_eps = [mp.mpf("1e-4"), mp.mpf("1e-6"), mp.mpf("1e-8")]
    residues = [eps * values(rho + eps, cutoff)[2] for eps in zero_eps]
    residuals = [abs(r + 1) for r in residues]
    assert all(residuals[i + 1] < residuals[i] for i in range(len(residuals) - 1)), residuals
    assert residuals[-1] < mp.mpf("1e-6"), residuals

    return {
        "U": u,
        "V": v,
        "near_one_values": near_one,
        "near_one_step_differences": one_steps,
        "rho": rho,
        "scaled_centered_at_zero": residues,
        "residue_residuals": residuals,
    }


if __name__ == "__main__":
    mp.mp.dps = 80
    print("Finite-cutoff analytic identities and local behavior (numerical only)")
    for u, v in ((1, 1), (3, 5), (8, 8)):
        row = check_pair(u, v)
        print(f"U={u}, V={v}")
        print("  centered values at 1+eps:", [mp.nstr(x, 16) for x in row["near_one_values"]])
        print("  successive differences:", [mp.nstr(x, 8) for x in row["near_one_step_differences"]])
        print("  first zero:", mp.nstr(row["rho"], 25))
        print("  eps*Ccal(rho+eps):", [mp.nstr(x, 18) for x in row["scaled_centered_at_zero"]])
        print("  residue residuals from -1:", [mp.nstr(x, 8) for x in row["residue_residuals"]])
    print("Numerical checks only; no assertion proves simplicity or RH.")
