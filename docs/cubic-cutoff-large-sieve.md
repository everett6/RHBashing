# Averaging the Type I error permits cubic-root cutoffs

**Status:** unconditional error estimate using the classical additive large
sieve. The final centered Type II estimate is still unproved. This note
improves the error treatment in this investigation; it makes no claim of
priority or a new result in the literature.

## 1. The improvement

Keep the definitions of S,T,m,R from `centered-type-ii-reduction.md`, but
now take U=V=floor(N^(1/3)), with even N≥8. The same exact identity gives

    S(x)−R(x)=e1(x)−e3(x),  N≤x<2N,

where |e1(x)|≤2V log(2N) and

    e3(x)=Σ_(k≤V,d≤U) μ(k)Λ(d)
                  ({(x−H)/(kd)}−{x/(kd)}),  H=N/2.

The previous pointwise estimate |e3|≤Vψ(U) was too large at these cutoffs.
Instead we prove the mean-square estimate

    Σ_(x=N)^(2N−1)|e3(x)|²
       ≪ (N+Q²)Q log²(2Q),  Q=UV.                 (1)

Since Q≤N^(2/3), the right side is O(N²log²(2N)). Also
Σ|e1|²≪NV²log²(2N)≪N^(5/3)log²(2N). Hence

    ||S−R||_2 ≪ Nlog(2N).                        (2)

The loss is logarithmic, which is acceptable for the every-epsilon RH
criterion proved in `global-energy-continuation.md`.

## 2. Fourier expansion of the floor error

Write e(t)=exp(2πit) and group products kd=q:

    c(q)=Σ_(kd=q,k≤V,d≤U) μ(k)Λ(d).

An especially useful bound is

    |c(q)|≤Σ_(d|q)Λ(d)=log q.                     (3)

There is no extra divisor-count factor: every divisor d determines k=q/d.
For integer x, the q-periodic sequence {x/q} has nonconstant Fourier
coefficients

    −1/[q(1−e(−j/q))],  1≤j<q.

This follows by differentiating a finite geometric sum. The constant
coefficient disappears from {(x−H)/q}−{x/q}. Reduce j/q to b/d with
(b,d)=1, 1≤b<d. All occurrences of that same frequency must be combined.
The resulting coefficient is exactly

    β_(b/d)=−[e(−bH/d)−1]/[1−e(−b/d)]
                    ·Σ_(l≤Q/d)c(dl)/(dl).         (4)

Thus e3(x)=Σ_(d≤Q)Σ_((b,d)=1,1≤b<d)β_(b/d)e(bx/d).
The q=1 summand contributes zero. Formula (4) avoids incorrectly treating
repeated rational frequencies as independent.

Let h_j=Σ_(l≤j)1/l. Using (3) without discarding the dependence on d,

    |Σ_(l≤Q/d)c(dl)/(dl)| ≤ log(2Q) h_floor(Q/d)/d.

Moreover |1−e(−b/d)|≥4 min(b,d−b)/d and |e(−bH/d)−1|≤2. Therefore

    |β_(b/d)|≪log(2Q) h_floor(Q/d)/min(b,d−b).

The sum of the squares of 1/min(b,d−b) over 1≤b<d is bounded by an
absolute constant. The remaining harmonic factors satisfy

    Σ_(d≤Q) h_floor(Q/d)²
      = Σ_(l,m≤Q) floor(Q/max(l,m))/(lm)
      ≤ Q Σ_(l,m≥1) 1/[lm max(l,m)] ≪ Q.

The last series converges, since it is at most
2Σ_(m≥1)h_m/m². Consequently

    Σ_(b/d reduced)|β_(b/d)|²≪Qlog²(2Q).           (5)

Keeping the harmonic factor before summing over d saves two powers of
log Q compared with bounding every h_floor(Q/d) by 1+log Q.

## 3. Apply the large sieve

Distinct reduced fractions with denominators at most Q have circular
spacing at least Q^(−2). The dual additive large-sieve inequality, on an
interval of N consecutive integers, bounds the energy of their exponential
sum by O(N+Q²) times the squared coefficient norm. Applying it to (5)
proves (1).

This uses a classical theorem, not an unproved assertion about prime or
Möbius correlations. A source giving the additive inequality and duality
is K. S. Kedlaya, *Notes on analytic number theory*, Chapter 15:
https://kskedlaya.org/ant/chap-largesieve.html .

## 4. What this changes in the RH target

With these cubic-root cutoffs the unknown bilinear sum T has both factors
r,k>N^(1/3), and because rk<2N, both are O(N^(2/3)). The earlier
quarter-power cutoff allowed a wider factor range. This narrows the
bilinear problem without leaving a power-sized Type I error.

The explicit centering remains

    m(x)=H(1+C+AL)−A∫_(x−H)^x log t dt,
    A=Σ_(k≤V)μ(k)/k,
    C=Σ_(k≤V)μ(k)log(k)/k,
    L=Σ_(d≤U)Λ(d)/d.

For R=T−m, equation (2) and the earlier dyadic/Mellin proof give

    RH ⇔ Σ_(x=N)^(2N−1)|R(x)|²≪N²log²(2N)
             uniformly over dyadic N.             (6)

Forward: RH gives ||S||_2≪N, so (2) gives the displayed bound for R.
Reverse: the displayed bound and (2) give F(N)=||S||²≪N²log²(2N),
which is O_ε(N^(2+ε)) for every ε>0, and hence implies RH.

The finite small N excluded by N≥8 can be absorbed into the constants.
We have proved the error theorem (1), not the arithmetic estimate (6).

## 5. A separate analytic route that fails

For fixed finite U,V define the Dirichlet polynomials
M_V(s)=Σ_(k≤V)μ(k)k^(−s), L_U(s)=Σ_(d≤U)Λ(d)d^(−s).
The Type II coefficient series initially on Re(s)>1 is

    Bcal(s)=(1−M_V(s)ζ(s))(−ζ'(s)/ζ(s)−L_U(s)).

Set A=M_V(1), C=−M_V'(1), D=1+C+A L_U(1). Subtracting the coefficient
main term D−A log n yields

    Ccal(s)=−ζ'/ζ+(M_V−A)ζ'+(M_V L_U−D)ζ−L_U.

At s=1, write t=s−1. The three potentially singular terms have residues
1, C, and −1−C respectively, so the singularity is removable. At any
zeta zero ρ of multiplicity m, all terms except −ζ'/ζ are holomorphic.
Thus Ccal has residue −m at ρ, independent of finite U,V.

This rules out removing a hypothetical off-line zero merely by tuning a
finite cutoff. It does not rule out bounds with cutoffs depending on N:
those are a family of transforms, and no uniform limiting argument follows
from this pole calculation alone. This route supplies no RH proof.

## 6. Verification and remaining work

An independent mathematical audit checked the Fourier coefficients,
grouping, separation, and large-sieve application. The implementation
`tests/check_type_i_large_sieve.py` compared the floor-sum error directly
against its grouped Fourier reconstruction for N=64,216,512,1000.
All checks passed; maximum absolute reconstruction error was 3.12e−12.
The numerical large-sieve inequality passed in all four cases.

A separate high-precision check script for the Dirichlet-series identities
was written but could not run because mpmath is absent from the available
Python runtimes. Only its syntax was checked; no numerical validation of
that script is claimed. The pole argument above is analytic.

The next required estimate is the centered Type II energy in (6). The
large sieve applied here controls the Type I remainder only. Applying the
same conclusion to T without an additional argument would be an error.
