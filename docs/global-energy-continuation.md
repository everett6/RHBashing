# Continuation: an RH-equivalent global energy criterion

**Status: proved equivalences and algebraic identities; RH is not proved.**
This continues `fourier-audit.md`. No novelty claim is made. The results below
correct an unnecessarily restrictive conclusion in the earlier discussion:
localized cross-frequency cancellation is not logically required for an RH
criterion if a subpolynomial loss is allowed.

## 1. A weaker bound is sufficient

Keep F,V,Q and the cutoff from the audit. For any fixed δ≥0, suppose
Q(N)≪N^(2+δ) uniformly over dyadic N. The parity approximation implies
F(N)≪N^(2+δ). Divide the recurrence

    sqrt(V(N)) ≤ sqrt(2V(N/2))+φ sqrt(F(N))

by N^(1+δ/2). Its contraction factor is 2^(−(1+δ)/2)<1, so
V(N)≪N^(2+δ). The real-variable mean square has the same exponent.
Cauchy–Schwarz on [X,2X] now gives

    ∫_X^(2X)|E(t)| t^(−σ−1)dt ≪ X^((1+δ)/2−σ).

Thus the Mellin integral is holomorphic for Re(s)>(1+δ)/2. As in the
audit, set Z(s)=(s−1)ζ(s); its logarithmic derivative there equals
−1−s∫_1^∞E(t)t^(−s−1)dt. This rules out zeros in that half-plane.
For 0≤δ<1, functional-equation symmetry confines nontrivial zeros to

    (1−δ)/2 ≤ Re(ρ) ≤ (1+δ)/2.

Consequently Q(N)≪_ε N^(2+ε) for EVERY ε>0 implies RH: for any proposed
zero with real part β>1/2 choose ε<2β−1. Constants may depend on ε but
not N. A bound for just one fixed positive ε does not prove RH.

The converse follows from the conditional log-free bound in the audit.
In particular proving Q(N)≪N²(log N)^A for any fixed A≥0 would suffice.
A bound Q(N)≪N²L(N), with 1≤L(N)=N^o(1), is also sufficient and is
implied by RH. No such unconditional bound is established here.

## 2. Restore the global Parseval route

Define the full-period energy, without a frequency cutoff,

    D(N)=Σ_(x=0)^(4N−1)|g(x)|²
        =(1/(4N))Σ_(k=0)^(4N−1)|â(k)W(k)|².

The equality is Parseval. Since F(N)≤D(N), a global estimate
D(N)≪_ε N^(2+ε) for every ε>0 implies RH by section 1's recurrence
and Mellin argument applied directly to F.

Conversely, under RH the classical pointwise estimate is
|E(t)|≪sqrt(t) log²(2t) for t≥1. Set T=2N−1 and E(0)=0. For
0≤x<4N, define

    u_x=min(max(x−H,0),T),  v_x=min(x,T).

There is no nonzero circular wrap contribution, because the support of a
is in [1,T] and H=N/2. Therefore exactly

    g(x)=E(v_x)−E(u_x).

Both endpoints are in [0,T]. Under RH this is uniformly
O(sqrt(N)log²(2N)); summing its square over 4N sites gives

    D(N)≪N²log⁴(2N).

We have proved the equivalences

    RH ⇔ D(N)≪N²log⁴(2N)
       ⇔ D(N)≪_ε N^(2+ε) for every ε>0,

where all bounds are uniform over dyadic N. The forward pointwise input is
classical, not newly proved here. References: NIST DLMF §25.16.4,
https://dlmf.nist.gov/25.16#E4 (the equivalent epsilon-loss formulation);
Brent–Platt–Trudgian's introduction for the log-squared pointwise estimate,
https://arxiv.org/abs/2008.06140.

This is a positive diagonal spectral criterion: it has no cross-frequency
terms. The earlier claim that global Parseval must be abandoned was too
strong. The earlier objection concerned a LOG-FREE O(N²) estimate. It
does not exclude logarithmic or arbitrary-small-power losses.

Do not reverse the reasoning circularly: the proof that RH implies this
estimate cannot be used as an unconditional proof of the estimate.

## 3. Exact arithmetic correlation target

Set a_n=Λ(n)−1 for 1≤n≤T. Since two intervals of H consecutive integers
shifted by d overlap in H−|d| sites for |d|<H, expansion of the full energy
gives the exact identity

    D(N)=H Σ_(n=1)^T a_n²
       +2 Σ_(d=1)^(H−1)(H−d) Σ_(n=1)^(T−d) a_n a_(n+d).

This can also be obtained by Parseval and the triangular autocorrelation
of the interval window. The diagonal is O(N²log N), already within the
acceptable logarithmic-loss scale. The unresolved task can thus be stated
as an upper bound of O(N²log⁴ N) on the signed off-diagonal aggregate.
There is no claim that individual shifted correlations obey such a bound.
Positivity of D provides a lower bound on this aggregate, not the required
upper bound. Taking absolute values loses relevant arithmetic information.

## 4. Vaughan identity, with all terms present

Use Dirichlet convolution *, let 1(n)=1, and let ε_0 be the convolution
identity (ε_0(1)=1, ε_0(n)=0 for n>1). Then μ*1=ε_0 and
log=1*Λ. Truncate μ at V and Λ at U. Splitting μ*log gives

    Λ = μ_≤V*log + Λ_≤U − μ_≤V*1*Λ_≤U + μ_>V*1*Λ_>U.

Indeed μ_>V*1=ε_0−μ_≤V*1, which proves the identity immediately.
The last term B equals μ_>V*b_U, with
b_U(r)=Σ_(d|r,d>U)Λ(d). Since Λ≥0 and Σ_(d|r)Λ(d)=log r,
0≤b_U(r)≤log r. Therefore |B(n)|≤τ(n)log n.

Elementary divisor estimates give Σ_(n<2N)|B(n)|²≪N log⁵(2N).
Restrict B to 1≤n<2N and extend it by zero. For a moving interval of
length H~N, Young's inequality only gives

    ||B*1_[0,H−1]||²_2 ≤ H²||B||²_2 ≪ N³log⁵(2N).

This is one power of N above the required scale. Using the equivalent
N^ε divisor estimates yields N^(3+ε), not N^(2+ε). An exact convolution
identity does not produce the missing saving. The −1 in a=Λ−1 must also
be retained when assembling the complete discrepancy.

## 5. Verification and limitations

`python3 tests/check_global_vaughan.py` requires NumPy and the companion
`check_fourier_reduction.py`. It checks 12 Vaughan decompositions at
N=16,64,128, including endpoint cutoffs; maximum coefficient discrepancy
was 2.66e−15. Five cases at N=16,64,128,512,2048 compare direct linear
convolution, global Fourier energy, and the triangular correlation sum.
All assertions passed; largest absolute correlation discrepancy was
9.74e−8 against energy about 1.15e6. These are floating-point consistency
checks, not asymptotic evidence sufficient for RH or formal verification.

The genuine change in the research target is that subpolynomial losses
permit global positive Fourier energy. The required arithmetic estimate
remains unproved. This is a reduction of the proof target, not an RH proof.
