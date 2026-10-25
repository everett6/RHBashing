# Actual-prime obstruction to separate Type II bounds

Status: unconditional elementary lemma, not a proof of RH. No novelty claim.
This tests the proposed step of bounding Type II window energies separately
and uniformly in Vaughan cutoffs. It disproves that particular step.

## Exact specialization

In the Vaughan identity from `global-energy-continuation.md`, take U=V=1.
Since Λ(1)=0 and μ(1)=1, its four terms reduce to

    Λ(n)=log n+B(n),  B=μ_{>1}*log=Λ−log.

This follows directly from μ*log=Λ. These are the actual arithmetic
coefficients; no synthetic replacement of Λ is being made.

For even N put H=N/2, I={N,...,2N−1}, and define

    L(x)=Σ_(r=0)^(H−1) log(x−r),
    T(x)=Σ_(r=0)^(H−1) B(x−r),
    S(x)=Σ_(r=0)^(H−1)(Λ(x−r)−1).

Then S(x)=(L(x)−H)+T(x) exactly.

## Theorem and proof

Uniformly for x in I, every x−r is between N/2 and 2N. Thus

    L(x)=H log N+O(N).

Nonnegativity of Λ and Chebyshev's estimate ψ(2N)=O(N) give
0≤Σ_r Λ(x−r)≤ψ(2N)=O(N). Consequently

    T(x)=−H log N+O(N)

uniformly in x. Squaring and summing over N sites proves

    Σ_(x∈I)|T(x)|² = (1/4)N³log²N+O(N³log N).

In particular, for any fixed 0<ε<1 this is not O(N^(2+ε)).
Therefore there cannot be an RH-scale bound of that form on the separate
Type II window energy that holds for all choices U,V≥1.

Put P(x)=L(x)−H. The same uniform estimates yield

    Σ P(x)² = (1/4)N³log²N+O(N³log N),
    Σ T(x)² = (1/4)N³log²N+O(N³log N),
    2Σ P(x)T(x) = −(1/2)N³log²N+O(N³log N).

The leading terms cancel in Σ S². Indeed, directly S(x)=O(N), so the
complete energy is O(N³), which is stronger than merely summing the
displayed error terms. Neither calculation proves an RH-scale bound.

This exhibits necessary mixed-term cancellation for this decomposition,
not just a failure of a particular application of Cauchy–Schwarz.

## Scope and next mathematical requirement

U=V=1 is a degenerate but valid choice. The theorem does NOT rule out an
estimate for a carefully selected pair of cutoffs growing with N. It does
NOT contradict the global diagonal Fourier criterion: a sum of absolute
squares of complete Fourier coefficients differs from splitting each
coefficient into arithmetic pieces and bounding their energies separately.

Subtracting the deterministic log term from T cannot finish the proof:
T+L−H is exactly S again. A useful next inequality must control the residual
after these large terms cancel, without assuming the desired prime-error
bound in the process. No such inequality is proved in this note.

## Independent implementation checks

Run `python3 tests/check_type_ii_obstruction.py` with NumPy and the two
companion check modules present. It constructs B from μ_{>1}*log and
independently compares it with Λ−log; then it verifies the moving-sum and
covariance identities. All four cases N=64,256,1024,4096 passed.

At N=4096 the ratio Σ T²/(N³log²N) was 0.20371018; the proved asymptotic
limit is 0.25. The complete energy was approximately 1,344,985.47.
Floating-point cancellation limits precision in the expanded energy.
The proof of the limit is the uniform estimate above, not these data.
