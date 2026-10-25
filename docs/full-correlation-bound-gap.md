# Actual prime correlations: an upper bound and its unresolved gap

**Status: the requested RH-scale upper bound is not proved.** This note
records a valid, weaker upper bound for the actual prime correlations and
an independently audited consequence of a hypothetical off-line zero.
Neither statement is presented as a new result in the literature.

## 1. The quantity being bounded

Let N≥8 be dyadic, H=N/2, T=2N−1, and I={N,...,2N−1}. Write

    E(t)=ψ(t)−t,
    S(x)=E(x)−E(x−H),
    F(N)=Σ_(x∈I) S(x)².

For any squarefree integer 2≤q≤N, including q depending on N, put

    ρ_q=q/φ(q),  w=ω(q),
    c_q(n)=(Λ(n)−ρ_q)·1_{gcd(n,q)=1},
    W_q(x)=Σ_(x−H<n≤x)c_q(n).

Define the overlap kernel

    K_N(n,m)=[min(2N−1,n+H−1,m+H−1)−max(N,n,m)+1]_+,

and the diagonal and full signed off-diagonal sums

    D_q(N)=Σ_(n=1)^T K_N(n,n)c_q(n)²,
    C_q(N)=2Σ_(h=1)^(H−1) Σ_(n=1)^(T−h)
                         K_N(n,n+h)c_q(n)c_q(n+h).

Expanding the square gives exactly

    ||W_q||₂²=D_q(N)+C_q(N).                       (1)

The previous note proves, with constants uniform over q≤N,

    ||S−W_q||₂≤√N B_q,  B_q=ρ_q2^w+log q=N^{o(1)},
    0≤D_q(N)≪N²log(2N).                            (2)

The target is a uniform upper bound C_q(N)≪N²log^A(2N) for some fixed
finite A, or C_q(N)≪_εN^{2+ε} for every ε>0. Equations (1)–(2) show
that the latter target is equivalent to the corresponding bound for F.
The dyadic/Mellin argument in the earlier notes then implies RH.

## 2. A proved upper bound, at a weaker scale

Define

    A₁(N)=9.39 log(2N)^1.515 exp(−0.8274√log(N/2)),
    f(t)=(log t)^(3/5)/(log log t)^(1/5),
    A₂(N)=0.026 log(2N)^1.801 exp(−0.1853 f(N/2)),
    A_N=A₁(N)                         if N<64,
    A_N=min(A₁(N),A₂(N))              if N≥64,
    P_N=79N³/12−5N²+2N/3.

Then the following bounds are unconditional:

    F(N)≤A_N²P_N,                                  (3)
    C_q(N)≤(A_N√P_N+√N B_q)²−D_q(N).              (4)

The external input is Johnston–Yang, *Some explicit estimates for the
error term in the prime number theorem*, Theorems 1.1 and 1.4:
https://arxiv.org/html/2204.01980v2 . Their respective estimates are

    |E(t)|≤9.39 t log(t)^1.515 exp(−0.8274√log t)        (t≥2),
    |E(t)|≤0.026 t log(t)^1.801 exp(−0.1853 f(t))       (t≥23).

We use established estimates, without asserting that these constants are
currently optimal.

To prove (3), both window endpoints lie in [N/2,2N]. Bound each increasing
logarithmic factor at 2N and each decreasing exponential factor at N/2.
For the second estimate, f is increasing on t≥23, since its logarithmic
derivative has the sign of 3−1/log log t. Consequently |E(t)|≤A_N t
throughout the endpoint range. Thus

    |S(x)|≤A_N(2x−H).

The exact finite sum

    Σ_(x=N)^(2N−1)(2x−H)²=79N³/12−5N²+2N/3

proves (3). The triangle inequality and (1)–(2) prove (4).

In particular, for some absolute c>0, these bounds give

    F(N)≪N³ exp(−c (log N)^(3/5)/(log log N)^(1/5)). (5)

The same upper-bound scale holds for C_q(N), uniformly in q≤N. Indeed,
the extra norm in (2) is N^{1/2+o(1)}, whereas the bound for √F is of
size N^{3/2−o(1)}; the former is smaller and is absorbed after reducing c.
This comparison concerns the upper-bound expressions, not a lower bound
on F.

For every fixed A, the ratio of the right side of (5) to N²log^A N
tends to infinity. Its logarithm is

    log N−c (log N)^(3/5)/(log log N)^(1/5)−A log log N →∞.

Thus (3)–(5) do not establish the requested scale. This is a bound on the
actual primes, but it is a consequence of existing PNT estimates rather
than a new solution of the correlation problem.

## 3. An exact quantitative test for a proposed proof

Suppose ζ has a zero z₀=β+iγ of multiplicity m with β>1/2. Set

    g=(1+√5)/2,
    κ_β²=(1−2^(−2β−1))/(2β+1),
    L=limsup_(N→∞, N dyadic) F(N)/N^(2β+1).

Then

    L≥[m log 2 (1−2^(−β))/(g |z₀| κ_β)]²>0.       (6)

Moreover, for every choice of squarefree q=q(N)≤N,

    limsup_(N dyadic) C_q(N)/N^(2β+1)=L.            (7)

These are conditional statements about what an off-line zero would
force. They do not assert that such a zero exists or prove its absence.

### 3.1. Recover the ordinary error from the window error

Put V(N)=Σ_(n=N)^(2N−1)E(n)². The following recurrence uses only linear
algebra and holds for every real sequence E:

    √V(N)≤√(2V(N/2))+g√F(N).                      (8)

For completeness, write M=N/2. For M≤t<2M set
d₁=E(t+M)−E(t) and d₂=E(t+2M)−E(t+M). Then

    (E(t+M),E(t+2M))=(E(t),E(t))+(d₁,d₁+d₂).

The Euclidean operator norm of the matrix [[1,0],[1,1]] is g. Taking
the norm over all t gives (8); the squared norm of all the d₁,d₂ is F(N).

If L=∞, (6) is immediate. Suppose L<∞. Divide (8) by N^(β+1/2).
The multiplier of the preceding normalized V term is 2^(−β)<1.
Iterating this contraction and taking upper limits gives

    limsup √V(N)/N^(β+1/2)≤g√L/(1−2^(−β)).        (9)

For the actual error E(t), E(n+u)=E(n)−u for 0≤u<1. Therefore, if
I_E(N)=∫_N^(2N)E(t)²dt, the norm difference satisfies

    |√I_E(N)−√V(N)|≤√(N/3).

It follows that (9) also holds with V replaced by I_E.

### 3.2. A boundary pole gives a positive lower bound for the limsup

Let

    b_j=∫_(2^j)^(2^(j+1)) |E(t)| t^(−β−1)dt.

Cauchy–Schwarz and the preceding continuous bound give

    limsup b_j≤κ_β g√L/(1−2^(−β))=:B.             (10)

The definition of κ_β is exactly the integral
κ_β²=∫_1²u^(−2β−2)du. In particular b_j is bounded. Hence

    M(s)=∫_1^∞ E(t)t^(−s−1)dt

converges absolutely and locally uniformly on Re(s)>β. Initially on
Re(s)>1, partial summation of the von Mangoldt Dirichlet series yields

    M(s)=−ζ′(s)/(sζ(s))−1/(s−1).                  (11)

The two sides agree by meromorphic continuation throughout Re(s)>β;
the left side is holomorphic there. At the boundary point z₀, the
meromorphic right side has residue −m/z₀. Thus, approaching horizontally
from the right,

    lim_(ε↓0) εM(z₀+ε)=−m/z₀.                    (12)

On the other hand, for ε>0,

    |εM(z₀+ε)|≤εΣ_(j≥0)2^(−jε)b_j.

For any η>0, (10) bounds all sufficiently late b_j by B+η. The finite
initial segment vanishes after multiplication by ε, and

    lim_(ε↓0) εΣ_(j≥0)2^(−jε)=1/log 2.

Combining this with (12) gives m/|z₀|≤B/log 2. Substitution of B proves
(6). This Abelian limit is essential: merely saying that M is holomorphic
on the open half-plane does not by itself rule out a pole on its boundary.

### 3.3. No admissible growing wheel hides this contribution

By (2), the normalized norm difference satisfies

    |√F(N)−||W_q||₂|/N^(β+1/2)≤B_q/N^β→0.

The upper limits of F(N)/N^(2β+1) and ||W_q||₂²/N^(2β+1) are therefore
equal, including when they are infinite. Also

    D_q(N)/N^(2β+1)≪N^(1−2β)log(2N)→0.

Subtract this vanishing diagonal in (1) to obtain (7). In particular,
conditional on an off-line zero, the full signed correlation itself has
positive excursions of order N^(2β+1) along arbitrarily large dyadic N.
No conclusion that this happens at every dyadic N is being asserted.

## 4. Why the attempted proof does not close

The periodic covariance estimate controls the sums of ν_q−1, with
ν_q(n)=ρ_q·1_{gcd(n,q)=1}. For r=Λ−ν_q, write R for its window sum and
δ for the window sum of ν_q−1. The exact relation is

    F=||R||₂²+2〈R,δ〉+||δ||₂².

A bound on ||δ||₂ does not supply an upper bound on ||R||₂. Its small
coefficient diagonal does not supply that bound either: correlations
between different coefficients remain. The CRT computation proves no
such cancellation estimate for Λ. Replacing the actual prime correlation
by the periodic covariance at that step would assume the missing result.

Statement (6) makes this gap quantitative, but is not a contradiction:
the unconditional upper bound (5) can accommodate N^(2β+1) for any
fixed β<1. The RH-scale upper bound has not been derived.

## 5. Verification

`tests/check_full_correlation_bounds.py` checks the finite polynomial P_N
exactly and compares the actual prime energy with (3) at finite sizes.
It also checks (8) for prime, deterministic random, and power sequences.
All 10 exact polynomial checks, 10 finite prime-bound comparisons, and
24 recurrence checks passed. Six further numerical checks compare the
power-sequence energy with its predicted scaled integral; their quadrature
is not an interval-arithmetic certification.
These are algebra and implementation checks. They do not establish the
asymptotic estimate, the cited external theorem, or any assertion about
zeta zeros. The mathematical proofs are the arguments above.

An independent internal audit checked (3)–(12), including the explicit
constant in (6) and the moving-modulus transfer in (7). This is neither
formal verification nor external peer review.
