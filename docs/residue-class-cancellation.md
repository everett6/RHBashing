# A necessary cancellation between residue classes

**Status:** unconditional arithmetic obstruction and exact normalizations.
This does not prove RH or improve the upper bound for the full prime-error
energy. It shows that a proposed separate bound for each shift class is
false, even after normalization by any fixed finite set of primes.
No claim of novelty in the literature is made.

## 1. Setup

Let N≥8 be dyadic, H=N/2, I={N,...,2N−1}, and

    S(x)=ψ(x)−ψ(x−H)−H.

All windows are (x−H,x], containing H integers. Their union is the integer
interval H<n<2N. Define the overlap kernel

    K(n,m)=[min(2N−1,n+H−1,m+H−1)−max(N,n,m)+1]_+.

It counts windows containing both n and m. In particular K≤H, and K=0
when |n−m|≥H. All coefficient sequences below are restricted to 1≤n<2N.

## 2. Remove the deterministic parity component exactly

Let a(n)=Λ(n)−1 and r(n)=Λ(n)−2·1_{n odd}. Their difference is alternating,
so every even-length H-window has the same sum for a and r.

Set c₂(n)=(Λ(n)−2)·1_{n odd}, with window sum W₂(x). The only even prime
power in H<n<2N is N itself. Therefore exactly

    S(x)=W₂(x)+(log 2)·1_{N≤x≤N+H−1},
    ||S−W₂||₂²=H(log 2)².                         (1)

All odd-shift correlations of c₂ vanish. For r, every nonzero odd-shift
pair with K>0 includes the exceptional site n=N. The absolute total of
these terms is at most

    2H log 2 · Σ_(n<2N,n odd)|Λ(n)−2| ≪ N²,

using Chebyshev's bound ψ(2N)≪N. This estimate applies to r, not to the
original a; the alternating adjustment redistributes correlation terms.

The c₂ diagonal is

    D₂=Σ_n c₂(n)²K(n,n) ≪ N²log(2N),              (2)

because Λ(n)²≤log(2N)Λ(n). Thus parity alone leaves only even shifts,
with an acceptable diagonal and negligible window error.

## 3. The shifts divisible by six have a large positive aggregate

Define the SIGNED aggregate

    O₆=2Σ_(1≤h<H,6|h)Σ_(n<2N−h)c₂(n)c₂(n+h)K(n,n+h).

For each a modulo 6 define

    C_a(x)=Σ_(x−H<n≤x,n≡a mod6)c₂(n).

Expanding the squares gives the exact identity

    O₆=Σ_(x∈I)Σ_(a mod6)C_a(x)²−D₂.             (3)

This includes all same-residue pairs and subtracts every diagonal once.
Focus on residue class 3. These are odd multiples of 3, so Λ can be
nonzero only at powers of 3. Every window has endpoint ratio at most 2,
and therefore contains at most one power of 3. The number of integers in
that residue class is H/6+O(1). Hence

    C_3(x)=−H/3+e_x,  |e_x|≤2+log 3.             (4)

All other squares in (3) are nonnegative. Squaring (4), summing N windows,
and using (2) proves

    O₆ ≥ NH²/9−O(NH)−D₂
       ≥ N³/36−O(N²log(2N)).                     (5)

In particular O₆ is positive and Ω(N³) for all sufficiently large dyadic N.
This is a lower bound for a SIGNED aggregate, not an absolute-value bound.

Consequences:

* One cannot bound each congruence class of shifts separately by an upper
  bound O(N²log^A N), for any fixed A, in this representation.
* One cannot prove simultaneous nonpositivity of every individual shift
  contribution at every sufficiently large N.
* The remaining shift classes must cancel this large positive aggregate
  if the full RH-scale estimate is to hold.

The positive shift may depend on N. This argument does not prove that one
particular fixed shift has a positive contribution for all large N.

## 4. Normalize modulo six, with a controlled error

Define

    c₆(n)=(Λ(n)−3)·1_{gcd(n,6)=1},
    W₆(x)=Σ_(x−H<n≤x)c₆(n).

Let P₂₃(x) be the Mangoldt weight of powers of 2 or 3 in the window, and
let δ(x)=3·#{n in the window:gcd(n,6)=1}−H. Then exactly

    S(x)=W₆(x)+P₂₃(x)+δ(x).                       (6)

Here 0≤P₂₃≤log6, because each window contains at most one power of each
prime. Since H mod6 is 2 or 4, δ lies in {−2,1} or {−1,2}, respectively.
Consequently

    |S(x)−W₆(x)|≤2+log6,
    ||S−W₆||₂²≤N(2+log6)².                       (7)

This removes the deterministic class-3 component from the previous
representation. It does not establish the residual energy bound.

## 5. Every fixed finite normalization leaves another obstruction

Let q be a fixed squarefree integer containing 2, let φ(q) be Euler's
totient, and set ρ=q/φ(q). Define

    c_q(n)=(Λ(n)−ρ)·1_{gcd(n,q)=1},
    W_q(x)=Σ_(x−H<n≤x)c_q(n).

Periodicity of the coprime indicator gives a discrepancy at most q after
normalization. Prime powers with primes dividing q contribute at most
Σ_(p|q)log p=log q per window. Thus

    |S(x)−W_q(x)|≤q+log q.                        (8)

For fixed q this is an acceptable O_q(sqrt N) norm error. The diagonal
D_q=Σ c_q(n)²K(n,n) is O_q(N²log(2N)).

Choose any fixed prime p not dividing q, and group shifts divisible by qp:

    O_(qp)=2Σ_(1≤h<H,qp|h)Σ_n c_q(n)c_q(n+h)K(n,n+h).

As in (3), this is the sum of squared window sums in residue classes
modulo qp, minus D_q. There are φ(q) residue classes a satisfying
p|a and gcd(a,q)=1. In each such class the only possible nonzero Λ values
are at powers of p. Thus its window sum equals

    −ρH/(qp)+O_{q,p}(1)
      =−H/(φ(q)p)+O_{q,p}(1).

Retaining only these φ(q) nonnegative squares proves

    O_(qp) ≥ NH²/(φ(q)p²)−O_{q,p}(N²log(2N))
           = N³/(4φ(q)p²)−O_{q,p}(N²log(2N)).     (9)

For every fixed q,p, this is a positive Ω_{q,p}(N³) aggregate at large N.
For example q=2,p=3 gives (5); q=6,p=5 gives a leading lower bound N³/200.

Both q and p are FIXED in this theorem. Constants and thresholds may
depend on them. The proof does not rule out a normalization q=q(N),
and does not supply the uniform estimates such a proposal would need.

This is a substantive restriction on a prospective proof: after any fixed
normalization, separate RH-scale upper bounds for every shift progression
are impossible. Cross-class signed cancellation must be retained, or a
different argument must account for these deterministic contributions.

## 6. Exact finite sign certificates

For N=16, the c₂ contribution at shift h=2 lies strictly between −77 and
−76, while the contribution at h=6 lies strictly between 37 and 38.
These signs are certified by `tests/certify_shift_sign.py`, using rational
arithmetic only. For 0≤t≤1/3, the program encloses the logarithm using

    log((1+t)/(1−t))=2Σ_(j≥0)t^(2j+1)/(2j+1).

After m terms its positive tail is at most
2t^(2m+1)/((2m+1)(1−t²)). Range reduction puts every required logarithm
in this interval. Exact interval multiplication and summation certify the
two rational brackets. This is a finite sign proof, not an RH proof.

## 7. Numerical checks and current boundary

The parity identities passed at N=16,64,256,1024,4096. The residue-square
formula and modulus-six correction also passed at those scales; direct
shift expansion was compared independently at N=16,64.
At N=4096, O₆ was approximately 2,808,766,518; its explicit class-3 lower
formula was approximately 1,837,343,557. These floating-point values are
implementation checks. The all-scale lower bound is the proof in (3)–(5).

Tests: `check_parity_reduction.py`, `check_residue_obstruction.py`, and
`certify_shift_sign.py`. `check_correlation_bands.py` separately audits
the earlier Vaughan-centered kernel decomposition; its coefficients differ
from the c₂ and c₆ coefficients here.

No bound on the full signed even-shift aggregate has been proved. The
remaining problem is the cancellation between large arithmetic shift
classes, with sufficient precision to leave O(N² times fixed logarithms).
