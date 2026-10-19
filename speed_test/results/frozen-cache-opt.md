# Frozen-bound cache optimization (2026-10-19)

The matched baseline contains the same concurrent gamma/generalized-exact-case
updates as the optimized version, with only this optimization reversed in an
isolated source snapshot. Both use wave speeds frozen per implicit solve. The optimized version
caches B per freeze generation and actual stage dt, stores B flux as scalar
identity multipliers, and caches HLL weights in Freeze. A/C, local LU and X
continue to refresh at every Newton state. No tolerances or method were changed.

Both versions were built on node2315 with the supplied `intel-build.sh`, using
`BUILD_DIR=build-opt-freeze-matched-before` and `BUILD_DIR=build-opt-freeze-final`.
The baseline script ran from an isolated snapshot. Temporary builds, snapshots,
drivers and raw test outputs were deleted at the user's request after validation;
the commands, binary hashes and measured summaries below are retained as the
final report. Referenced build names describe the completed runs, not existing
build directories.
Compiler: IntelLLVM 2025.2.0; Release/native; CPU: Xeon Platinum 8268.

Executable SHA256:

- Before: `b9dcd51c3489de37e0c9418d21e8efd80f725f0c45ac2d0b50337f26abd16f89`
- After: `6dce15509c5ae3f1875199f61943ab11e1ec75c04aea49fca60be5123efb8105`

## Correctness

`python3 code_test/test_solver.py --exe build-opt-freeze-final/euler2d
--output code_test/results/frozen-opt-suite` passed both existing tests.
The converged BE final DG/HDG maximum nodal difference was 1.67e-16.

One-off comparisons used 8x8 cells, CFL .1, T=.04, Newton/linear tolerance
1e-12, both formulations, saved coefficients and independent update checks.
They exercised BE/ILU, BE/full, DIRK/partial, DIRK/element and BE/Q8.
All succeeded; maximum matched before/after coefficient difference was 2.17e-15.
These tests cover repeated stages/steps and clipped final steps, and do not
extend the maintained two-test code suite. The temporary per-case outputs have
been cleaned up.

## Timing and process memory

Workload: 16x16 Q4, BE, CFL=1, one step, exactly three Newton updates,
partial/block-Jacobi, linear/Newton tolerances 1e-12, timing-only.
OMP/MKL/OpenBLAS/BLIS/vecLib thread counts were all one. Each version/form
was warmed once, then measured seven times, alternating version and form order.
The following are medians; time is elapsed seconds per accepted update,
not process CPU time. RSS covers the complete solver process.

| Form | Before mean Newton (s) | After mean Newton (s) | Before/after | Before peak RSS (KiB) | After peak RSS (KiB) |
|---|---:|---:|---:|---:|---:|
| DG | .3820115225 | .3837176546 | .99555 | 52428 | 52568 |
| HDG | .2289353748 | .2304172828 | .99357 | 48108 | 47668 |

On this small matched case DG update time increased about .45% and HDG about
.65%; no elapsed-time speedup is established. Earlier measurements exhibited
substantial drift and also straddled concurrent source updates; those preliminary
measurements were excluded from benefit claims and have been cleaned up. These small timing
changes and RSS differences are local observations, not a large-case performance
or memory-capacity guarantee. Additional HLL weights cost memory on DG;
HDG additionally benefits from compact B storage. B flux payload itself is
1/16 of the former payload; total RSS does not have that reduction factor.

Raw samples, resource JSONs and one-off drivers have been removed as requested;
only this final summary remains. Repeating the measurements requires rebuilding
the versions and running new cases.
