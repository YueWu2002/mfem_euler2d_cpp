# MFEM Euler: equivalent DG and condensed HDG

Serial C++17 comparison of backward-Euler DG and its algebraically equivalent
frozen-HLL HDG reformulation. The principal timing metric is the **mean elapsed
seconds per accepted Newton update inside the time loop**. No numerical tests or
performance benchmarks were run for this refactor.

## Discretization and case

- Euler gamma=3 simple wave, default angle 45 degrees. Change `--alpha-deg` in
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
- Einfeldt/Roe HLL bounds are frozen from the previous time level, including
  `--safety 1.1`. Physical fluxes and their derivatives use the current state.
- Initial data use an accurate local L2 projection. Exact smooth errors are
  evaluated before shock formation, `T < 1/pi`; the default endpoint is `0.2`.

See [MATHEMATICS.md](MATHEMATICS.md) for equations, inverse-mass weighting, and
DG/HDG update equivalence; [DESIGN.md](DESIGN.md) for storage and performance.

## Build

Use a configured **double-precision, serial MFEM** installation/build:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DMFEM_DIR=/path/to/mfem/lib/cmake/mfem -DEULER_NATIVE=ON
cmake --build build --config Release -j
```

`EULER_NATIVE=ON` enables `-march=native` for GCC/Clang. Build MFEM itself in
Release too. Native fixed-size tensors are exposed by `mfem-performance.hpp`.
The optional sparse direct solver requires MFEM built with SuiteSparse.
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

# Matrix-free Jacobian/Schur multiplication, assembled native block ILU(0) PC.
./build/euler2d --nx 16 --ny 16 --k 3 --cfl 1 --steps 3 --precond block-ilu0 --assembly partial

# Fully assembled multiplication with the same residual/quadrature.
./build/euler2d --nx 16 --ny 16 --k 3 --dt 0.001 --steps 1 --assembly full

# Explicitly selected direct reference; HDG globally solves only its Schur system.
./build/euler2d --nx 4 --ny 4 --k 2 --steps 1 --solver direct --form hdg

# Independent DG equation and one-update comparison, outside benchmark mode.
./build/euler2d --nx 4 --ny 3 --k 2 --steps 1 --newton-steps 1 --form hdg --check-update
```

`--cfl C` overrides `--dt`; define `h=min(hx,hy)` and
`amax=max_faces,q(max(-SL,SR))`, **including safety**, so `dt=C*h/amax`.
There is no hidden `2k+1` factor. `--steps N` overrides `--T` and does not clip
the last step. In T mode the last step is clipped, so its actual CFL can be
smaller. Every chosen endpoint must remain below `1/pi`.

`--newton-steps N` performs exactly N damped Newton updates per time step;
this is a controlled workload, not a claim of nonlinear convergence. Without
it, use `--max-newton` (25) and `--newton-tol` (1e-10). The convergence norm is
the physical cell dual L2 residual divided by dt and by sqrt(domain area).
Both paths use the same positivity/Armijo line search on the DG residual.

GMRES uses `--restart 20`, `--max-krylov 1500`, `--linear-tol 1e-11`.
Right PCs are `block-jacobi` (default) and `block-ilu0`. For ILU, use
`--ilu-ordering none` (default, fixed block ordering) or `mdf`. MDF's numerical
ordering need not be invariant under a change of basis. `--assembly partial`
(default) caches quadrature Jacobians for DG and factored local Schur actions
for HDG; `element` stores local matrix actions; `full` also assembles global
CSR and uses it for multiplication. `--matrix-free` selects `partial`.
Native ILU always needs a global CSR PC matrix, even with matrix-free action.

`--solver dense` is an explicit diagnostic limited to 5000 global unknowns.
Direct modes require full assembly, select it automatically if unspecified,
and reject explicit matrix-free options.

## Timings and files

`--output PREFIX` defaults to `results/run`. DG/HDG append `_dg`/`_hdg`.

| Output | Contents |
|---|---|
| `*_summary.csv` | Actual final time, DOFs, flags, counts, `newton_s`, **`mean_newton_s`**, errors, conservation drift |
| `*_iterations.csv` | One accepted update per row; residual, damping, Krylov count, linear residuals, stage times, `iteration_s` |
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

```sh
python scripts/test_solver.py --exe ./build/euler2d --quick
python scripts/test_solver.py --exe ./build/euler2d --full-time
python scripts/benchmark.py --exe ./build/euler2d --nx 16 --ny 16 --k 3 --cfl 1 --steps 3 --newton-steps 3
```

The test driver uses only the Python standard library. It checks k=0..8 algebra,
dense nonorthogonal basis changes, DG/HDG one-update equivalence, partial/element/
full actions, both PCs, CFL/angle/geometry, fixed update counts, timing accounting,
conservation and CLI rejection; `--full-time` additionally advances to 0.2.
The optional CTest entry invokes `--verify` only when CTest is explicitly run.
SuiteSparse direct mode is not required by the Python tests.

The benchmark driver warms each path once, alternates execution order, fixes
thread counts for common BLAS libraries, checks update counts and reports the
median of each run's mean Newton time. Failures stop the comparison rather
than being included as zero-time samples. High-order local LU memory/cost can
still dominate; measure before drawing DG/HDG speed conclusions.

## Refactor validation status

All application C++ translation units passed `g++ -fsyntax-only` against the
provided `C:/Users/20333/Desktop/mfem-master` headers, using a temporary minimal
serial/double configuration in `tmp/headercheck`. Both Python drivers passed
AST syntax parsing without execution. This is an interface/syntax check, not
a linked MFEM build, a numerical test, or a performance measurement. The optional
SuiteSparse branch still requires checking with an actual SuiteSparse-enabled
MFEM installation. Historical `RESULTS.md` and `lint.log` are not current evidence.
