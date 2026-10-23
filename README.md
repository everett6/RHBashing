# RHBashing

CUDA exploration of the Riemann zeta function on the critical line, targeted at
the NVIDIA RTX 5070 desktop (Blackwell, `sm_120`). This is a finite-range
**candidate scanner**, not an RH proof or certified zero counter. It cannot
search all numbers or detect zeros off the critical line. Sampling can miss
closely spaced zeros and zeros with even multiplicity.

It approximates Hardy's real-valued function
`Z(t) = exp(i theta(t)) zeta(1/2 + i t)` using the Riemann–Siegel main sum and
leading correction. The omitted correction is asymptotically `O(t^-3/4)`;
this is not a certified error bound. Signs close to zero can be wrong.
The supported exploratory range is `200 <= t <= 1e9`, using double precision.

## RTX 5070 build and run

Install an NVIDIA driver and CUDA Toolkit **12.8 or later** supporting Blackwell,
and GNU Make plus a CUDA-supported C++ compiler. No API keys or services are used.

```sh
make cuda
make test-cuda  # numerical checks and CPU/GPU bracket parity on your GPU
./build/rhscan-cuda 200 10000 0.05 > candidates.csv
```

Arguments are `START END STEP`, finite numbers with `END > START` and `STEP > 0`.
Rows contain sign changes of the approximation, with endpoints and approximate
Z values. The final endpoint is always sampled. Output streams in bounded
batches; timing goes to stderr. Smaller steps can reveal additional candidates,
but cannot certify completeness. Cost grows with sample count and `sqrt(t)`.
The example is a starting point, not a measured performance claim.

```sh
# Independent high-precision checks (can be much slower than the scan).
python3 -m venv /tmp/rhbashing-venv
/tmp/rhbashing-venv/bin/pip install mpmath==1.3.0
/tmp/rhbashing-venv/bin/python tools/verify.py candidates.csv > checked.csv
```

The checker reevaluates endpoints using `mpmath.siegelz` at 50 decimal digits,
then bisects sign changes. It fails if a candidate does not bracket a sign
change. These are numerical checks, not interval proofs. They cannot discover
zeros the scanner missed. A complete zero count needs additional certified
methods such as Turing's method, absent here.

## CPU reference and tests

```sh
make cpu
./build/rhscan-cpu 200 240 0.05 > candidates.csv
make test  # needs mpmath==1.3.0 in the active Python environment
```

CPU and CUDA use the same formula and CLI. CPU evaluation uses OpenMP. Tests
compare values against independent high-precision values, exercise removable
correction singularities, and check known zero brackets.

## CUDA design

One warp evaluates each sample: lanes sum every 32nd term and reduce with warp
shuffles. A contiguous `(log(n), 1/sqrt(n))` table is computed once on the host,
copied once, and reused across batches. Blocks have eight warps. Device buffers
and pinned host output are reused. Only scalar Z values return to the CPU, which
emits brackets across batch boundaries. Memory is bounded by 65,536 samples plus
the term table. Batches shrink at high t to limit work per kernel launch on
desktop display GPUs. There are no atomics or per-sample device allocations.

Consumer Blackwell has limited FP64 throughput. Double precision and accurate
trigonometric functions are intentional; `--use_fast_math` is not enabled.
Benchmark on the actual card before expanding a search. This is not an
implementation of high-altitude algorithms such as Odlyzko–Schönhage, nor a
replacement for a certified zeta computation package.

Formula reference: NIST DLMF §25.10, <https://dlmf.nist.gov/25.10>.
The theta expansion and leading correction are in `src/riemann.hpp`.
