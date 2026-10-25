# Unconditional bounds for small and boundary separations

Status: partial estimates for the Vaughan-centered coefficients, not RH.
Use c_N(n)=B_N(n)−d_N+a_V log n with U=V=floor(N^(1/3)),
B_N(n)=Σ_(rk=n,r>U,k>V) μ(k)Σ_(d|r,d>U)Λ(d),
a_V=Σ_(k≤V)μ(k)/k, and
d_N=1+Σ_(k≤V)μ(k)log(k)/k+a_VΣ_(d≤U)Λ(d)/d.
Use the same window-overlap kernel K. Do not confuse this sequence with
the parity-normalized coefficients in `residue-class-cancellation.md`.

We have |c_N(n)|≤τ(n)log(2N)+O(log²(2N)). An elementary shifted-divisor
bound gives, uniformly for 1≤h<2N,

    Σ_(n<2N−h)|c_N(n)c_N(n+h)|
        ≪ Nlog⁴(2N) σ_{−1}(h),

where σ_{−1}(h)=Σ_(d|h)1/d. To justify the divisor estimate, bound τ(n)
by twice its divisors ≤sqrt(2N), and do the same for n+h. Divisibility
by d and e is compatible only if gcd(d,e)|h. Compatible solutions have
count at most 2N/lcm(d,e)+1. The main terms are bounded by
O(Nlog²(2N)Σ_(g|h)1/g), and the +1 errors by O(N). Cross terms involving
the O(log²N) part use Στ(n)≪NlogN. This proves the stated bound.

Since K≤H=N/2 and Σ_(h≤L)σ_{−1}(h)≤ζ(2)L, the absolute contribution
from shifts 1≤h≤L is

    O(N²Llog⁴(2N)).

Taking L=log N keeps this within O(N²log⁵N). This improves the logarithm
in the coefficient-energy-only estimate, but does not remove the factor L.

At the opposite boundary K(n,n+h)≤H−h. For H−M≤h<H the contribution
is bounded by

    O(Nlog⁴(2N)·M·[M+log(2N)]).

Indeed, for any interval of M shifts below 2N, summing 1/d over multiples
gives Σσ_{−1}(h)≤ζ(2)M+O(log(2N)). Taking M=floor(sqrt N) therefore
controls this boundary band by O(N²log⁴N).

An independent coefficient-energy argument gives the explicit bound
M(M+1)Σ_n c_N(n)² for the same absolute boundary contribution. The
script `check_correlation_bands.py` checks this bound and the full signed
kernel identity at six scales. The middle band of O(N) shifts remains
uncontrolled at RH scale; these estimates do not save a power there.
