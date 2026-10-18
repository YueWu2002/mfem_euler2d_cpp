# Implementation and performance contracts

## Module responsibilities

| Module | Responsibility |
|---|---|
| `euler.hpp` | Types, dimensions, reusable workspace interfaces and ownership |
| `physics.cpp` | Gamma=3 flux/Jacobian, Einfeldt speeds, rotated exact wave |
| `problem.cpp` | MFEM geometry/GL basis, residuals, local blocks, metrics |
| `partial.cpp` | Refreshable DG quadrature Jacobian cache and analytic action |
| `operators.cpp` | Element/factored actions, diagonal panels, optional CSR |
| `condensation.cpp` | Preallocation, reusable cell LU, Schur action and recovery |
| `linear.cpp` | Persistent weighted right GMRES and MFEM LU/ILU/direct PCs |
| `stepper.cpp` | Time/Newton loops, common residual criterion, stage timers |
| `verification.cpp` | Independent algebra and dense basis-change diagnostics |
| `scripts/test_solver.py` | Manually invoked regression suite |
| `scripts/benchmark.py` | Warmup, alternating pairs and median mean-Newton time |

## Layout and reusable storage

Cell layout is `(cell, native MFEM GL node, conservative component)`; face
layout is `(global periodic face, GL node, conservative component)`.
`State=mfem::TVector<4>` and `Matrix4=mfem::TMatrix<4,4>` have inline contiguous
storage. TMatrix uses MFEM's column-major convention; access it with `(i,j)`.
Their entries are explicitly initialized before accumulation. Degree-dependent
blocks use native MFEM DenseMatrix/Vector, with allocation separate from refresh.

Before the time loop, allocate:

- residual, RHS, correction, trial and line-search residual vectors;
- partial DG's physical Jacobians and diagonal block storage;
- HDG A/B/C/D blocks, restrictions, X/y, LU storage, Schur and recovery scratch;
- GMRES Arnoldi basis, Hessenberg/Givens arrays and all application vectors;
- block-Jacobi diagonal storage and pivoted LU buffers.

Numerical refresh, all factorizations, Jacobian evaluation and necessary matrix
products still occur at the current Newton state **inside** the Newton timer.
Preallocation does not move numerical assembly work out of the comparison.
The default partial/block-Jacobi path reuses the above buffers through the run.
Non-owning MFEM Vector views access contiguous Jacobi blocks without gather,
scatter or per-block heap storage. DenseMatrixInverse pointer constructors
allocate factors without factoring an uninitialized matrix; `Factor()` copies
and refactors the current values in those same buffers.

Partial DG caches both directional volume flux Jacobians and HLL face
Jacobians. Block-Jacobi volume assembly reuses these caches instead of
re-evaluating the state and Jacobian. Volume/face evaluation uses a precomputed
mapping to native MFEM node order and sparse GL restriction/gradient tables:
volume values are direct nodal reads, volume gradient work is O((k+1)^3), and
face evaluation uses only its k+1 nonzero basis functions per face point.
Mass inverse application in the collocated path is a direct diagonal scaling.
The partial DG diagonal builder applies these inverse-mass row weights while
accumulating sparse face terms directly into the cell block, avoiding two
full dense face-scratch clears/scales for every face.

HDG B and D depend on frozen wave speeds and dt, not on the Newton state.
A freeze-generation counter refreshes them once per time step; A and C refresh
at every Newton state. Nodal trace collocation lets A/B/C face terms accumulate
directly in their destination entries, without a dense face matrix or a second
mass-scaling pass. This counter is valid because production changes dt
only immediately after Freeze, before the first HDG refresh of that step.
Public callers changing dt themselves must Freeze again before reuse.

Schur diagonal extraction multiplies only matching face panels of C and X.
Its inner loop streams contiguous columns rather than strided dot products.
Partial Schur action stores no C*X product; only explicit element/full modes do.
Recovery reuses cached X/y and scratch, with no global recovery matrix.

## Krylov and ownership

Weighted restarted GMRES uses twice-modified Gram-Schmidt and Givens rotations.
Its dot product represents the inverse-mass weak residual metric after undoing
stored-equation row scaling (see MATHEMATICS.md). It recomputes the true residual
at restart/candidate convergence. Right PC factors stay fixed during one solve. Schur tolerance refinement
continues from the previous correction with the same factors; it does not
repeat assembly or PC setup.
Since the PC is fixed and linear, the final Arnoldi linear combination is
preconditioned once; a second full `Z` basis is unnecessary. This reduces memory
without changing right-preconditioned GMRES algebra. There is no left-Jacobi
preconditioning that would change the minimized norm.

Each Run owns Hybrid before Condensed. Condensed's operators borrow stable
pointers to Hybrid C and Condensed X; they are destroyed before their owners.
The block vectors are not resized to different lengths after PrepareHybrid.
DenseMatrixInverse holds pointers to stable DenseMatrix objects. Replacing
these vectors or retaining an action after its workspace dies is invalid.

Reusable scratch and mutable physical-Jacobian temporaries make a Problem/
workspace single-threaded and non-reentrant. A future OpenMP implementation
needs per-thread scratch and race-free face accumulation, not simply a pragma
around the current face loop. No MPI support is included.

## What still allocates and why

The optional explicit DG element/full paths rebuild their local operator
objects. Optional CSR assembly builds and reduces contributions into a fresh
matrix. Native MFEM BlockILU::SetOperator reconstructs its symbolic block pattern
(including sets in the supplied MFEM source); retaining the BlockILU object does
not avoid this internal work. The supplied BlockILU::Mult also creates a temporary block-sized Vector on each
application. Sparse direct solver internals also manage their
own allocations. These costs remain inside Newton assembly/setup timings.
`--check-update` and verification deliberately allocate independent reference
objects; use them for correctness, not production benchmarking.

For maximum serial performance, start with Release MFEM/application and
`--assembly partial --precond block-jacobi --timing-only`. Whether ILU's reduced
Krylov count compensates for its assembly/setup cost is a measured question.
The native ILU route is retained as requested; a custom block-CSR ILU with a
cached symbolic pattern would be a separate implementation beyond native MFEM's
current API. It could reduce setup allocation but should be validated separately.

## Timing boundary

The mean-Newton numerator starts before that update's residual evaluation and
ends after accepted line search/state update. It includes Jacobian construction,
condensation, required assembly/PC factorization, weighted GMRES, recovery and
common true-residual checks. Optional independent update-reference time is
subtracted and reported only as a diagnostic value; do not benchmark with it.
The accepted line-search forward residual is reused at the next Newton state;
it is not evaluated twice, and its work remains charged to the update that
computed it. CSV conversion/write occurs after the update timer. `--timing-only` removes
iteration output entirely. Freeze/CFL and convergence-only final residual
checks are per-step overhead, represented in `time_loop_s`, not assigned to an
accepted Newton update. Mean is `newton_s/newton_updates`, never an assumed
constant iteration count or a total-runtime ratio.

## Further performance work requiring measurement or changed assumptions

1. Benchmark MFEM local LU with a **single-thread** LAPACK/BLAS backend against
   its built-in kernels, especially k=5..8. Library choice can dominate local LU.
2. Profile time by stage before specializing fixed k or fusing residual and
   Jacobian evaluation. A fused pass could remove state/pressure recomputation,
   but must retain separate residual-only evaluation for line search.
3. Cache explicit CSR symbolic structure and numerical assembly slots. This
   helps full/ILU runs, but needs an assembly API that preserves block structure
   and explicit zero entries across Newton states.
4. Tensor/operator-specific local factorizations or approximate local solves
   could reduce high-k cost. An approximate A^{-1} would change exact condensation
   and update equivalence unless used only as a PC; it is **not enabled**.
5. Overintegration, LGL nodes, reduced face traces, modified stabilization,
   mixed precision, fast-math, and looser tolerances can change the numerical
   method or equivalence guarantees. They are **not enabled** for speed claims.

The refactor was not numerically tested or benchmarked at the user's request.
Historical RESULTS.md and lint.log predate these changes and are not validation
or performance evidence for this version.
