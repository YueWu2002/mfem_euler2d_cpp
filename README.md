# MFEM Euler: equivalent DG and condensed HDG

Serial C++17 comparison of SSP-DIRK3 DG and its algebraically equivalent
frozen-HLL HDG reformulation. The principal timing metric is the **mean elapsed
seconds per accepted Newton update inside the time loop**. Numerical tests and
benchmarks are explicit commands, never build-time side effects.

## Project objectives

1. **Verify correctness and DG/HDG equivalence.** With the same DIRK scheme,
   HLL bounds, quadrature, and tolerances, compare recovered body-DOF updates,
   final solutions, conservation, and convergence. Equivalence is up to solver
   tolerances and roundoff: DG solves the body system; HDG solves the condensed
   face system and recovers the body update by back-substitution.
2. **Measure and improve HDG efficiency.** Establish speedups through isolated
   single-thread paired tests at equal accuracy, reporting mean accepted
   Newton-update time and converged total-solve time, including all
   setup, condensation, recovery, and common DG-residual work. Report the
   configurations where HDG is faster or slower; DOF reduction alone is not
   a timing result.

Convergence tests use the production DIRK for both methods; other integrators
are not comparison baselines. Independent convergence cases may run in
parallel, but performance comparisons must run in isolation. See
[RESULTS.md](RESULTS.md) for evidence and current test status.

## Discretization and case

The equation of state is created in `src/main.cpp` as `EulerGas gas(3.)`.
`EulerGas` owns protected, immutable `gamma` and `gamma_minus_one` members;
its pressure, flux, Jacobian, and HLL methods read them directly. `Gamma()`
provides read-only access. `Problem` owns a copy shared by its DG/HDG paths,
with no process-global gas state. The ideal-gas kernels support gamma > 1;
the rotated exact-solution formulas also use variable gamma. The selected test
problem supplies gamma when constructing `EulerGas` (currently 3 in main);
gamma is not a CLI option and stays constant for that problem's run.
CSV output reports the configured value.

- Euler simple wave, problem-specified gamma=3 and default angle 45 degrees. Change `--alpha-deg` in
  `(0,90)`; the centered domain is
  `[-1/cos(alpha),1/cos(alpha)] x [-1/sin(alpha),1/sin(alpha)]`.
- `nx x ny` affine quadrilateral Qk elements, native MFEM Gauss-Legendre nodal
  basis; degrees `0..8`, including Q8. MFEM supports arbitrary-order L2 square
  elements; the upper bound 8 is this application's supported range.
- Exactly `(k+1)^2` tensor-product GL volume points and `k+1` GL face points.
  Face Pk traces use the nodal basis at these same points.
- Periodicity is imposed by Cartesian face pairing, including opposite domain
  boundaries. Exported MFEM geometry is the centered rectangular mesh; its mesh
  topology itself is not converted to a periodic MFEM mesh.
- Einfeldt/Roe HLL bounds (`--safety 1.1`) are fixed throughout each implicit
  solve, including all Newton iterations and line-search trials. Backward
  Euler uses the beginning-of-step state; each DIRK stage uses its initial
  guess (the previous stage solution after stage zero). Stage zero reuses
  the beginning-of-step CFL estimate. DG and HDG use the same policy.
  Newton differentiates the residual with these bounds held constant; flux
  Jacobians still follow the current iterate. GMRES and line search make this
  an inexact, damped Newton solve of the fixed-bound stage equation.
- The default `ssp-dirk3` uses a four-stage third-order DIRK tableau.
  `--integrator be` retains backward Euler as a diagnostic.
- Initial data use an accurate local L2 projection. Exact smooth errors are
  evaluated before shock formation, `T < 4/((gamma+1)*pi)`; the default endpoint is `0.2`.

See [MATHEMATICS.md](MATHEMATICS.md) for equations, inverse-mass weighting, and
DG/HDG update equivalence; [DESIGN.md](DESIGN.md) for storage and performance.

## Build

On the OSCAR Intel nodes, use the supplied optimized build script:

```sh
./intel-build.sh
```

Otherwise use a configured double-precision MFEM installation/build:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DMFEM_DIR=/path/to/mfem/lib/cmake/mfem -DEULER_NATIVE=ON
cmake --build build --config Release -j
```

`EULER_NATIVE=ON` uniformly enables `-march=native` for
Intel/IntelLLVM/GCC/Clang/AppleClang; IntelLLVM/GCC also compile SIMD
directives without enabling threaded OpenMP. Build MFEM itself in Release.
The optional `--solver direct` path requires MFEM built with MKL PARDISO.
No tests execute during configure or build.

`make clean` or `cmake --build build --target clean` also removes `lint.log`
from the project root and the build directory. Reconfigure an existing build
directory once to pick up this clean rule.

## Run

```sh
# Converge every implicit step through T=0.2; weighted right block Jacobi GMRES.
./build/euler2d --nx 16 --ny 16 --k 3 --cfl 1 --T 0.2 --form both

# A controlled timing workload: 3 time steps, exactly 3 Newton updates each.
./build/euler2d --nx 16 --ny 16 --k 3 --cfl 1 --steps 3 --newton-steps 3 --form both --timing-only

# Matrix-free Jacobian/Schur multiplication, direct BSR block ILU(0) PC.
./build/euler2d --nx 16 --ny 16 --k 3 --cfl 1 --steps 3 --precond block-ilu0 --assembly partial

# Fully assembled multiplication with the same residual/quadrature.
./build/euler2d --nx 16 --ny 16 --k 3 --dt 0.001 --steps 1 --assembly full

# PARDISO reference for a small one-step update check only.
./build-intel/euler2d --nx 8 --ny 8 --k 2 --steps 1 --dt 1e-4 \
  --newton-steps 1 --solver direct --form hdg --check-update

# Independent DG equation and one-update comparison, outside benchmark mode.
./build/euler2d --nx 8 --ny 8 --k 2 --steps 1 --newton-steps 1 --form hdg --check-update
```

`--cfl C` overrides `--dt`; define `h=min(hx,hy)` and
`amax=max_faces,q(max(-SL,SR))`, **including safety**, so `dt=C*h/amax`.
There is no hidden `2k+1` factor. `--steps N` overrides `--T` and does not clip
the last step. In T mode the last step is clipped, so its actual CFL can be
smaller. Every chosen endpoint must remain below `4/((gamma+1)*pi)`.

`--newton-steps N` performs exactly N damped Newton updates per DIRK stage
(four stages per time step for the default method);
this is a controlled workload, not a claim of nonlinear convergence. Without
it, use `--max-newton` (30) and `--newton-tol` (1e-12). The convergence norm is
the physical cell dual L2 residual divided by dt and by sqrt(domain area).
Both paths use the same positivity/Armijo line search on the DG residual.

GMRES uses `--restart 20`, `--max-krylov 1500`, `--linear-tol 1e-12`.
Right PCs are `block-jacobi` (default) and `block-ilu0`. For ILU, use
`--ilu-ordering none` (default) or `mdf` (one-time MDF on the unweighted fixed
topology graph). The graph, ordering and ILU update schedule are prepared before
the time loop; numerical assembly/factorization stay inside Newton timing.
`--assembly partial`
(default) caches quadrature Jacobians for DG and factored local Schur actions
for HDG; `element` stores local matrix actions; `full` also assembles global
CSR and uses it for multiplication. `--matrix-free` selects `partial`.
Production ILU receives cell/face BSR blocks directly, with no intermediate
scalar CSR on the matrix-free path. HDG B/C use shared Qk tensor face restrictions and compact flux
blocks; cached X=A^{-1}B avoids repeated local solves during GMRES. Jacobi and
production HDG local LU overwrite their matrices and retain partial pivoting.
HDG uses one streamed cell A/LU workspace, D stores only diagonal entries, and
Schur Jacobi evaluates only matching face panels. No local solve is repeated
to achieve these memory reductions.
DG matrix-free volume multiplication uses explicit one-dimensional tensor-line
derivative contractions (O((k+1)^3)); face restriction/lifting is O((k+1)^2).
The exact HDG local inverse is not assumed to be tensor separable.
High-order ILU factors can still be large; these changes do not guarantee that
256x256 Q8 fits in 64 GiB. See DESIGN.md for storage and timing details.

`--solver dense` is an explicit diagnostic limited to 5000 global unknowns.
`--solver direct` means MKL PARDISO. Direct modes require full assembly, select
it automatically if unspecified, and reject explicit matrix-free options. Use
them for small correctness references, not production timings.

## Timings and files

`--output PREFIX` defaults to `code_test/results/run`. DG/HDG append `_dg`/`_hdg`.

| Output | Contents |
|---|---|
| `*_summary.csv` | Integrator/stages, final time, DOFs, counts, `newton_s`, **`mean_newton_s`**, errors, conservation drift |
| `*_iterations.csv` | One accepted update per row; DIRK stage, residual, damping, Krylov count, component times, `iteration_s` |
| `*_comparison.csv` | Final nodal and density differences; ratio of mean Newton times |
| `*.mesh`, `*.gf`, `*_coefficients.csv` | Centered geometry and native GL nodal solution |
| `*_verification.csv`, `*_weighted_gmres.csv` | Explicitly requested algebra/basis-invariance diagnostics |

`iteration_s` includes residual/forward evaluation, Jacobian refresh, local
factorization/condensation, any global assembly, PC construction, Krylov work,
recovery, required residual validation and line search. It also includes
ordinary bookkeeping and allocations that still occur within that update.
Only storage preparation occurs before the time loop. Initial projection,
freeze/CFL selection once per time step, convergence-only terminal checks,
final errors/export, CSV formatting/writes and optional independent
`--check-update` reference solves are excluded from the mean Newton timer.
`time_loop_s` is a separate wall-clock diagnostic and includes step overhead.

`--timing-only` suppresses iteration CSV and solution files; summary and paired
comparison remain. It still computes the required true residual checks.
Do not use `--check-update` for benchmarking: subtracting its elapsed time does
not remove cache effects. Fixed workloads can eventually reach rounding-level
residuals; avoid requesting many superfluous Newton updates.

## Tests and paired benchmarks (run manually)

`code_test/` contains only a smoke test and a backward-Euler, one-step,
converged-Newton DG/HDG update-equivalence test;
`accuracy_test/` contains final-time error and convergence-order studies; and
`speed_test/` contains timing and resource experiments. Each driver defaults to
the output area owned by its category.

All test categories use `--newton-tol 1e-12` and `--linear-tol 1e-12`.
The solver defaults and the accuracy driver's defaults match these settings;
the correctness and speed drivers pass them explicitly. Fixed-update speed
tests still use their prescribed Newton update counts, so the Newton tolerance
does not terminate those updates; the linear tolerance remains active.
Existing binaries must be rebuilt to acquire changed C++ defaults, and stored
results retain the tolerances actually used when they were generated.

```sh
python code_test/test_solver.py --exe ./build/euler2d
python speed_test/benchmark.py --exe ./build/euler2d --nx 16 --ny 16 --k 3 --cfl 1
python accuracy_test/convergence.py --exe ./build-intel/euler2d --grids 8 16 32 --T .2 --cfl .1
# Independent cases may run concurrently; resulting times are not benchmarks:
python accuracy_test/convergence.py --exe ./build-intel/euler2d --grids 32 64 128 256 --T .2 --cfl .5 --jobs 8
```

The test driver uses only the Python standard library. Its smoke case uses one
CFL-selected DIRK step. Its only equivalence case uses one CFL-selected
backward-Euler step, lets Newton converge normally, and checks every recovered
HDG update against the independent DG update.

The benchmark driver warms each path once, alternates execution order, fixes
thread counts for common BLAS libraries, checks update counts and reports the
median of each run's mean Newton time. Failures stop the comparison rather
than being included as zero-time samples. High-order local LU memory/cost can
still dominate; measure before drawing DG/HDG speed conclusions.

Speed tests use one backward-Euler step (`--integrator be --steps 1`) and
exactly three Newton updates (`--newton-steps 3`), for both DG and HDG.
This controlled workload compares mean accepted-Newton-update wall time;
three updates do not imply nonlinear convergence. Accuracy tests retain DIRK.
Slurm cases write to `speed_test/results/be_newton3/job_<ID>/k<k>_<precond>/`.
Submit with `bash speed_test/submit.sh`; each submission creates its own
manifest under `speed_test/results/be_newton3/`. After completion, run
`python speed_test/collect.py` for the CSV/Markdown summaries in that directory,
or add `--update-results` to update the backward-Euler section of `RESULTS.md`.
Existing DIRK measurements retain their original paths and configuration.

## Validation status

The current Intel/MKL build has been compiled and numerically exercised. See
[RESULTS.md](RESULTS.md) for exact commands, tolerances, errors, convergence
rates, and timings. Generated CSV/log artifacts under `code_test/results/`,
`accuracy_test/`, and `speed_test/results/` are ignored.
