# Growing moduli: exact local covariance and the unresolved prime residual

**Status:** the deterministic divisibility-model covariance is bounded
unconditionally, including moduli growing with N. The full correlation
upper bound for the actual prime coefficients is NOT proved.
The local-density model used here is standard; no novelty claim is made.

## 1. Uniform normalization, not a fixed-modulus limit

Let N≥8 be dyadic, H=N/2, and I={N,...,2N−1}. Let q≤N be squarefree,
w=ω(q), and ρ=q/φ(q). Define, for 1≤n<2N,

    ν_q(n)=ρ·1_{gcd(n,q)=1},
    c_q(n)=(Λ(n)−ρ)·1_{gcd(n,q)=1},
    W_q(x)=Σ_(x−H<n≤x)c_q(n),
    S(x)=ψ(x)−ψ(x−H)−H.

The dependence q=q(N) is allowed. Put

    δ_q(x)=Σ_window ν_q(n)−H,
    P_q(x)=Σ_(window,gcd(n,q)>1)Λ(n).

Then exactly S=W_q+δ_q+P_q. Inclusion–exclusion gives

    |δ_q(x)|≤ρ2^w.

The ratio of the window endpoints is at most 2. For every prime p dividing
q, the window therefore contains at most one power of p. Hence

    0≤P_q(x)≤Σ_(p|q)log p=log q,
    |S(x)−W_q(x)|≤ρ2^w+log q.                    (1)

These bounds have no hidden q-dependent constants.

Order the w distinct prime divisors as p₁<...<p_w. Since p_j≥j+1,

    ρ=Π_j p_j/(p_j−1)≤Π_j(j+1)/j=w+1,
    q≥(w+1)!.

The factorial bound implies w=O(log N/log log N), uniformly q≤N.
Consequently the right side of (1) is N^{o(1)}, uniformly over all such q:

    ||S−W_q||₂²≤N(ρ2^w+log q)²=N^{1+o(1)}.      (2)

In particular ||S−W_q||₂=o(N). The normalization is still negligible at
RH scale even when q is comparable with N.

## 2. The coefficient diagonal stays controlled uniformly

Let T=2N−1. Inclusion–exclusion also bounds the coprime count by

    #{n≤T:gcd(n,q)=1}≤T/ρ+2^w.

Using Λ²≤log(2N)Λ and Chebyshev's bound,

    Σ_(n≤T)c_q(n)²
       ≤2log(2N)ψ(T)+2ρT+2ρ²2^w.

The factorial inequalities imply ρ²2^w≤6q: the ratio
(w+1)²2^w/(w+1)! is at most 6 (its maximum is attained at w=2, and it
decreases thereafter). Also ρ≤w+1=O(log N). Therefore

    Σ_(n<2N)c_q(n)²≪Nlog(2N)                     (3)

with an absolute constant, uniformly q≤N. The window-kernel diagonal is
at most H times this, hence O(N²log N).

This is a coefficient/diagonal estimate. It does not bound the full
window energy, which includes off-diagonal correlations.

## 3. The complete covariance of the periodic model is explicit

For an integer shift h, define

    S_q(h)=Π_(p|q) [ p/(p−1),                    if p|h;
                     p(p−2)/(p−1)²,             if p does not divide h ].

The Chinese remainder theorem counts residue pairs avoiding zero modulo
every p|q and gives the EXACT period-average identity

    (1/q)Σ_(x mod q)ν_q(x)ν_q(x+h)=S_q(h).

The mean of ν_q is 1 and S_q(0)=ρ. Expanding the square of a window sum
therefore proves, for every integer H≥1,

    (1/q)Σ_(x mod q)δ_q(x)²
      = H(ρ−1)+2Σ_(h=1)^(H−1)(H−h)(S_q(h)−1)
      ≤ρ²4^w.                                    (4)

The upper bound follows from the pointwise inclusion–exclusion estimate.
Thus (4) rigorously accounts for cancellation between all shift classes
in the finite-prime periodic model. No conjecture about prime pairs is
used. For q≤N the bound is N^{o(1)}.

The displayed covariance identity averages over a full period. For the
actual interval I, no period alignment is needed: the pointwise bound gives

    Σ_(x∈I)δ_q(x)²≤Nρ²4^w=N^{1+o(1)}.            (5)

One must NOT replace correlations of Λ by S_q(h). Formula (4) is exact
for ν_q only. That substitution would introduce an unproved prime-pair
estimate, rather than prove it.

## 4. Uniform version of the positive shift-class obstruction

Choose a prime p not dividing q. Let O_(qp) denote the signed off-diagonal
aggregate of c_q at shifts divisible by qp, with the same overlap kernel
K as in the previous note. Let D_q be its nonnegative diagonal.

The exact residue-square identity gives O_(qp)+D_q as the sum over x∈I
of the squared window sums in all residue classes modulo qp.
Select the φ(q) classes divisible by p and coprime to q. Together their
window coefficient sum is

    −H/p+e_x,  |e_x|≤E:=ρ2^w+log p.

This uses inclusion–exclusion on multiples of p, and at most one p-power
in the window. Cauchy–Schwarz over the selected classes yields

    O_(qp) ≥ [N/φ(q)](H/p−E)_+²−D_q.             (6)

Unlike the previous fixed-wheel asymptotic, (6) is uniform in its parameters.

For any fixed 0<η<1, suppose q≤N^(1−η). A prime p not dividing q can be
chosen with p=N^{o(1)}: among the first w+1 primes one is missing, and
Bertrand's postulate gives p_{w+1}≤2^{w+1}=N^{o(1)}.
Then E=N^{o(1)} and H/p=N^{1−o(1)}. In view of (3), (6) implies

    O_(qp) ≥ (1−o(1))N³/[4φ(q)p²]
             ≥ N^(2+η−o(1)).                     (7)

Thus a modulus polynomially smaller than N still leaves a positive
superquadratic shift-class aggregate. Separate RH-scale upper bounds on
every progression remain impossible there. This does not refute the full
signed bound, which permits cancellation across progressions.

At q=N^{1−o(1)} the displayed main term need not dominate the diagonal.
That is where this obstruction loses force; it is not a proof that the
required upper bound then holds.

## 5. Attempted transfer of a published uniformity estimate

I checked Matomäki–Radziwiłł–Shao–Tao–Teräväinen, *Higher uniformity of
arithmetic functions in short intervals II. Almost all intervals*,
Inventiones mathematicae 244 (2026), 967–1091:
https://doi.org/10.1007/s00222-026-01408-6 .
The paper uses a related coprimality approximant and gives logarithmic
savings in short-interval sums, with exceptional sets. Its approximant
parameters differ from the q≤N setup here.

The displayed H/log^A X scale does not supply the needed square-root
scale for H comparable with N. Even an everywhere estimate N/log^A N
would give only N³/log^(2A)N upon squaring and summing N windows. For
every fixed A and B,

    [N³/log^(2A)N]/[N²log^B N]
       =N/log^(2A+B)N →∞.

Allowing any fixed A does not authorize A=A(N); constants and thresholds
may depend on A. Exceptional intervals would require additional estimates
as well. This comparison concerns the displayed quantitative estimates,
not a claim that the paper's methods can never be strengthened.

## 6. What has and has not been bounded

* All periodic-model covariance: proved bound (4), with (5) on our interval.
* Error in growing-wheel normalization: proved bound (2).
* Actual residual coefficient diagonal: proved bound (3).
* Superquadratic positive progression contributions for q≤N^(1−η): proved (7).
* Full window energy of c_q: no RH-scale upper bound established.

The remaining quantity is Σ_(x∈I)|W_q(x)|². Equivalently, it is the full
signed off-diagonal aggregate of the actual residual c_q, after including
the controlled diagonal. Neither the periodic calculation nor the
published log-saving estimate fills this gap.

## 7. Verification

`check_periodic_model.py` verified the autocorrelation product, variance
identity and bound for 12 pairs q∈{2,6,30,210}, H∈{4,16,64}, using exact
rational arithmetic. For q=210,H=64 the variance is exactly 617/48.

`check_growing_wheel.py` verified the prime-residual decomposition,
normalization inequalities, factorial estimates and coefficient-energy
bound for N=16,64,256,1024,4096. Those tests use floating-point logarithms;
the uniform claims are proved above, not inferred from the data.

Independent mathematical review checked (1)–(7). This is internal review,
not external peer review or formal verification. RH is not proved here.
