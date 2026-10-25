# Removing the Type I terms at the required error scale

**Status:** unconditional approximation theorem; its remaining Type II
estimate is unproved. This is not a proof of RH and no novelty is claimed.
Unlike a statement that mixed covariance must be controlled, this note
specifies the deterministic cancellation and bounds all discarded terms.

## 1. Statement

For even N≥2 set H=N/2 and U=V=floor(N^(1/4)). For N≤x<2N put a=x−H.
All x are integers. Define

    A = Σ_(k≤V) μ(k)/k,
    C = Σ_(k≤V) μ(k)log(k)/k,
    L = Σ_(d≤U) Λ(d)/d,
    J(x) = ∫_(x−H)^x log t dt
         = x log x−x−(x−H)log(x−H)+(x−H),
    b_U(r) = Σ_(d|r,d>U) Λ(d),
    T(x) = Σ_(r>U,k>V,a<rk≤x) b_U(r)μ(k).

The smooth centering function for T is

    m(x) = H(1+C+AL)−A J(x).

Let R(x)=T(x)−m(x), and let S(x)=ψ(x)−ψ(x−H)−H. Then

    |S(x)−R(x)| ≤ 2V log(2N)+Vψ(U).                 (1)

The estimate is unconditional, uniform in x, and retains exact endpoints.
Since ψ(U)≪U and U=V=floor(N^(1/4)), its right side is O(sqrt N).
Consequently

    | (Σ_(x=N)^(2N−1) S(x)²)^(1/2)
      −(Σ_(x=N)^(2N−1) R(x)²)^(1/2) | ≪ N.        (2)

Thus the original RH-scale mean-square problem is equivalent to bounding
the energy of this explicitly centered Type II sum. Type I errors have
already been bounded at a sufficient scale; they are not an additional
unresolved covariance term in this formulation.

## 2. Proof from the exact Vaughan identity

For Dirichlet convolution *, write

    Λ = μ_≤V*log + Λ_≤U − μ_≤V*1*Λ_≤U
        + μ_>V*1*Λ_>U.

This follows from log=1*Λ and μ*1 equal to the convolution identity.
The last term is precisely the coefficient whose window sum is T(x).
Since U≤H≤a and the interval is (a,x], the Λ_≤U term vanishes on it.
Let I1,I3 be the window sums of the first and third terms. Then exactly

    S(x)=I1(x)−I3(x)+T(x)−H.                       (3)

### The logarithmic Type I term

For 1≤k≤V, put α=a/k and β=x/k. We have 1≤α≤β. Partial summation gives

    Σ_(α<m≤β) log m −∫_α^β log t dt
      = −{β}log β+{α}log α+∫_α^β {t}/t dt
      = η_k,

where {t}=t−floor(t). The endpoints in the discrete sum are important.
Taking absolute values bounds |η_k| by
log β+log α+log(β/α)=2log β≤2log(2N).
Substitution t=y/k in the integral therefore yields

    I1(x)=Σ_(k≤V)μ(k)Σ_(a/k<m≤x/k)log m
         = A J(x)−HC+e1(x),
    |e1(x)|≤2Vlog(2N).                             (4)

No cancellation in μ has been assumed here.

### The divisor Type I term

Counting multiples gives exactly

    I3(x)=Σ_(k≤V,d≤U)μ(k)Λ(d)
            [floor(x/(kd))−floor(a/(kd))].

The bracket differs from H/(kd) by
{a/(kd)}−{x/(kd)}, of absolute value less than 1. Hence

    I3(x)=HAL+e3(x),
    |e3(x)|≤Σ_(k≤V)|μ(k)|Σ_(d≤U)Λ(d)≤Vψ(U).       (5)

Combining (3)–(5) gives S=R+e1−e3 and proves (1).
Since Vlog(2N)=O(sqrt N) and UV≤sqrt N, Chebyshev's estimate proves
the uniform square-root error. The Euclidean norm of the error on N
windows is O(N), proving (2) by the reverse triangle inequality.

## 3. Consequence and precise remaining assertion

The earlier audit proved that F(N)=Σ S(x)²≪N² uniformly over dyadic N
is equivalent to RH. By (2), we therefore have the criterion

    RH ⇔ Σ_(x=N)^(2N−1)|T(x)−m(x)|²≪N²
          uniformly over dyadic N.                 (6)

Likewise the bound ≪_εN^(2+ε) for every ε>0 is equivalent to RH, as shown
by the dyadic recurrence and Mellin argument in the continuation note.
A bound with any fixed power of log N in place of N^ε would suffice.

Equation (6) is NOT proved unconditionally here. The proved result is
the deterministic centering and its error bound (1). Computing m(x)
requires only μ up to V and Λ up to U; it contains no unknown prime
counting error at scale N.

This resolves a particular design issue: one need not seek a bound on the
uncentered T, which can contain a large deterministic component, nor leave
all Type I–Type II covariances unspecified. The unresolved estimate is
now a variance bound around the explicitly displayed function m(x).
Its proof still requires arithmetic cancellation not supplied by Vaughan's
identity alone.

## 4. Why the cutoff choice matters

The proof works more generally for integer U,V≤N/2, with the same bound
2Vlog(2N)+Vψ(U). Therefore sufficient conditions for the RH-scale error are

    Vlog(2N)=O(sqrt N),   UV=O(sqrt N).

Balanced quarter-power cutoffs satisfy both. The elementary argument with
U=V=N^(1/3) would instead leave a pointwise error O(N^(2/3)), giving an
energy bound O(N^(7/3)); that error cannot be discarded at the target scale.
This is a limitation of this error estimate, not an impossibility theorem
for other treatments of the Type I sums.

For every fixed 0<α<1/2, the choices U≈N^α and V≈N^(1/2−α) also meet
the asymptotic conditions (with constants depending on α).

## 5. Attempt to finish by absolute values

The elementary estimate |B(n)|≤τ(n)log n, for B=μ_>V*b_U, bounds its
coefficient energy by O(Nlog⁵N). Young's inequality for an H-term window
then bounds Σ|T(x)|² only by O(N³log⁵N). Subtracting m after this inequality
cannot recover the lost signed cancellation. This is why neither (1) nor
an absolute-value divisor estimate proves (6).

The next substantive requirement is an estimate for the centered bilinear
sum itself, uniform over dyadic N. A finite table of small values, or an
independence assumption about the Möbius signs, cannot replace that proof.

## 6. Verification

An independent mathematical audit checked the quadrature, signs, endpoints,
cutoff conditions, and norm implication. This is not formal verification.
Run `python3 tests/check_centered_type_ii.py` with NumPy and the existing
companion modules. Six cases N=2,16,64,256,1024,4096 passed:

* Exact Vaughan reconstruction from independently built convolutions.
* The identity S−R=e1−e3 for the two Type I errors.
* The explicit error bound in (1), including the smallest even N.

At N=4096, U=V=8, max|S−R| was 16.7301 and the explicit bound was
198.042. The error norm was 375.063. Floating-point checks confirm the
implementation on these inputs; the proof of (1) is the argument above.

Related files: `fourier-audit.md`, `global-energy-continuation.md`,
`type-ii-obstruction.md`.
