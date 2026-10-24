# Fourier reduction audit — 24 October 2026

**Status: a proved reduction and finite consistency checks, not a proof of RH.**
The universal low-frequency bound below has not been established. This note
does not claim novelty. Independent model review is not formal verification.

## Definitions

Write E(t)=ψ(t)−t, where ψ(t)=Σ_{n≤t}Λ(n). For even N set

    F(N) = Σ_{x=N}^{2N−1} |E(x)−E(x−N/2)|²,
    V(N) = Σ_{x=N}^{2N−1} |E(x)|².

All discrete sums over x use integers. Let L=4N, H=N/2 and define a(n)=Λ(n)−1
for 1≤n<2N, zero elsewhere on Z/LZ. Use the unnormalized forward DFT
â(k)=Σ_n a(n)exp(−2πikn/L) and inverse normalization 1/L.
Let W(k)=Σ_{r=0}^{H−1}exp(−2πikr/L) and g=a*1_[0,H−1].
On I={N,…,2N−1}, g(x)=E(x)−E(x−H).

Put K=ceil(sqrt(N log(2N))) and d(k)=min(k,L−k), 0≤k<L.
Let h be the inverse transform of â(k)W(k) restricted to d(k)≥K.
Let P be the inverse transform restricted to d(k)<K and even k, and put
Q(N)=Σ_{x∈I}|P(x)|². This is precisely the even low-frequency energy.

## 1. Unconditional Fourier tail and parity reduction

For k≠0, the geometric sum and sin(πd/L)≥2d/L give
|W(k)|≤L/(2d(k)). Also

    Σ_n |a(n)|² ≪ N log(2N),

since Λ(n)²≤log(2N)Λ(n) and Chebyshev's estimate gives ψ(2N)≪N.
Parseval therefore proves

    ||h||²_(Z/LZ) ≤ (4N²/K²) Σ_n |a(n)|² ≪ N².

For x∈I, g(x+2N)=0: all arguments of the moving window lie above the
support of a. The even low modes are unchanged by this shift and odd low
modes change sign. Adding g(x) and g(x+2N) gives the exact identity

    g(x) = 2P(x)+h(x)+h(x+2N).

The two index blocks are disjoint, so the norm of the last two terms on I
is at most sqrt(2)||h||. Thus

    |sqrt(F(N))−2sqrt(Q(N))| ≪ N.

In particular Q(N)=O(N²) if and only if F(N)=O(N²), with uniform constants.
This is not an estimate of Q itself.

## 2. Full proof that the dyadic variance bound implies RH

Write N=2M. For t=M,…,2M−1 set
d1(t)=E(t+M)−E(t), d2(t)=E(t+2M)−E(t+M).
Their combined squared norm is F(2M). The vector of E values on [2M,4M)
is the sum of (E(t),E(t)) and (d1(t),d1(t)+d2(t)). The matrix
[[1,0],[1,1]] has Euclidean operator norm φ=(1+sqrt(5))/2. Therefore

    sqrt(V(N)) ≤ sqrt(2V(N/2))+φ sqrt(F(N)).

Assuming F(N)≤CN² for every dyadic N, u_N=sqrt(V(N))/N satisfies
u_N≤u_(N/2)/sqrt(2)+φ sqrt(C). Iteration proves V(N)≪N².

For real t∈[n,n+1), E(t)=E(n)−(t−n), so
∫_n^(n+1)|E(t)|²dt≤2|E(n)|²+2/3.
Thus ∫_X^(2X)|E(t)|²dt≪X² on dyadic blocks. For σ>1/2,
Cauchy–Schwarz yields

    ∫_X^(2X)|E(t)| t^(−σ−1) dt ≪_σ X^(1/2−σ).

Summing dyadic blocks shows J(s)=∫_1^∞ E(t)t^(−s−1)dt is absolutely
convergent and holomorphic for Re(s)>1/2, locally uniformly there.
For Re(s)>1, partial summation gives

    −ζ'(s)/ζ(s) = s/(s−1)+sJ(s).

The right side is meromorphic in Re(s)>1/2 with its only possible pole at 1.
Uniqueness of meromorphic continuation makes this an identity throughout
that half-plane. A zero ρ of ζ there, of multiplicity m, would give the left
side a pole with residue −m, impossible on the right. Hence ζ has no zeros
with Re(s)>1/2. The functional equation reflects every nontrivial zero left
of 1/2 to one right of 1/2. All nontrivial zeros therefore lie on the line.

Conversely, under RH, Cramér's mean-square bound gives
∫_X^(2X)E(t)²dt≪X². The reverse unit-interval comparison gives V(N)≪N².
Using |a−b|²≤2|a|²+2|b|² and covering [N/2,3N/2) by dyadic blocks
then gives F(N)≪N². See Brent–Platt–Trudgian, *The mean square of the
error term in the prime number theorem*, Journal of Number Theory 238
(2022), https://arxiv.org/abs/2008.06140.

Consequently, with constants uniform over dyadic N,

    RH ⇔ F(N)=O(N²) ⇔ Q(N)=O(N²).

The conditional mean-square theorem is used only for RH ⇒ the bound;
it cannot establish the bound unconditionally.

## 3. Why the Fourier estimates cannot supply the missing bound

Replace a by +1 on 1,…,N and −1 on N+1,…,2N−1, with the same padding.
Its coefficient energy is 2N−1 and its coefficient sum is 1. For x=N+t,
the moving sum is H−2t for 0≤t<H and −H for H≤t<N. Hence exactly

    F_synthetic(N)=Σ_(t=0)^(H−1)(H−2t)²+H³=N³/6+N/3.

For this sequence the same tail bound improves to O(N²/log N).
The parity identity then gives
2sqrt(Q_synthetic)=sqrt(F_synthetic)+O(N/sqrt(log N)), and therefore

    Q_synthetic(N)=(1/24+o(1))N³.

It follows that support, small total sum, coefficient second moments,
and the displayed Fourier identities do not imply the target bound.
This counterexample is not a counterexample to RH: it is not Λ−1.

## 4. Corrections to earlier claims

* Littlewood's oscillation theorem over all x does not by itself put large
  errors at x=2N−1 with dyadic N. The earlier dyadic endpoint obstruction
  was unsupported. A separate subsequence argument would be required.
* Frequencies α=j/(2N) tending to zero are not automatically in a chosen
  zero major arc. An arc width must be specified; convergence to zero alone
  is insufficient.
* Rewriting a divisor sum as its DFT proves an identity, not cancellation.
  The elementary Type II Parseval estimate does not control the required
  localized, weighted combination with Type I terms.
* Negative covariance is not established merely by decomposing an energy.
  An exact restatement of Q(N)≤CN² is still the same unproved assertion.
* Numerical agreement on finite N verifies implementation consistency, not
  a uniform bound for all dyadic N or a proof of RH.

## 5. Reproducible verification

Run `python3 tests/check_fourier_reduction.py` with NumPy installed.
The script builds Λ by a prime-power sieve and independently compares direct
moving sums with the Fourier parity reconstruction. It checks the explicit
tail inequality and the exact synthetic energy. Tests ran successfully for
N=16,64,128,512,2048 for both prime and synthetic coefficients.
The largest absolute reconstruction discrepancy was 6.82e−13.
At N=2048 the synthetic Q/N³ was 0.0416664700168, consistent with 1/24.
These are floating-point checks, not interval certification.

## Outstanding requirement

An unconditional proof of Q(N)≪N² for the actual von Mangoldt coefficients
has not been obtained. No completed RH proof follows from this report.
