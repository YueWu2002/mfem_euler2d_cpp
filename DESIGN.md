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
| `linear.cpp` | Persistent weighted right GMRES, numeric PC refresh and direct solves |
| `block_ilu.cpp` | Fixed BSR graph/order/update schedule, in-place block ILU(0) |
| `stepper.cpp` | SSP-DIRK3/BE stages, Newton loops, common residual criterion and timers |
| `verification.cpp` | Independent algebra and dense basis-change diagnostics |
| `code_test/test_solver.py` | Smoke test and converged backward-Euler update-equivalence test |
| `accuracy_test/convergence.py` | Final-time convergence-order study |
| `speed_test/benchmark.py` | Warmup, alternating pairs and median mean-Newton time |

Convergence diagnostics may run independent parameter cases concurrently, but
the application has no distributed MPI mesh/vector/halo implementation. An MPI
launcher alone duplicates serial solves; each independent process must receive
distinct parameters and output paths. DG/HDG timing comparisons remain isolated
single-thread runs. Correctness and equal accuracy are prerequisites for any
HDG speedup claim; all setup, condensation, and recovery costs must be counted.

## Layout and reusable storage

`EulerGas` owns protected const gamma and gamma-minus-one members, initialized
by its constructor in main. `Gamma()` is a read-only accessor. Each `Problem`
owns an immutable gas copy; DG/HDG use that same equation of state. Kernels
are const methods, with no mutable global gas state or per-call gamma argument.
`--verify` also tests independent gas objects at gamma=1.4 and 1.67 without
altering the production gas configuration.

Cell layout is `(cell, native MFEM GL node, conservative component)`; face
layout is `(global periodic face, GL node, conservative component)`.
`State=mfem::TVector<4>` and `Matrix4=mfem::TMatrix<4,4>` have inline contiguous
storage. TMatrix uses MFEM's column-major convention; access it with `(i,j)`.
Their entries are explicitly initialized before accumulation. Degree-dependent
blocks use native MFEM DenseMatrix/Vector, with allocation separate from refresh.

Before the time loop, allocate:

- residual, RHS, correction, trial and line-search residual vectors;
- partial DG's physical Jacobians and diagonal block storage;
- HDG diagonal D, tensor B/C couplings, restrictions, X/y, one cell A/LU buffer;
- GMRES Arnoldi basis, Hessenberg/Givens arrays and all application vectors;
- block-Jacobi diagonal storage (overwritten by LU) and pivot arrays;
- four stage slopes and the DIRK affine-base/state work vectors.

Numerical refresh, all factorizations, Jacobian evaluation and necessary matrix
products still occur at the current Newton state **inside** the Newton timer.
Preallocation does not move numerical assembly work out of the comparison.
The default partial/block-Jacobi path reuses the above buffers through the run.
Non-owning MFEM Vector views access contiguous Jacobi blocks without gather,
scatter or per-block heap storage. Jacobi and production HDG cell factors use
MFEM LUFactors directly on their matrix storage, with partial pivoting. Each
Newton refresh overwrites the previous LU with the new matrix before factoring.
No second dense LU buffer or additional intra-block reordering is retained.
Production HDG reuses a single cell matrix and pivot array: assemble A, factor,
solve for X/y, then overwrite the workspace for the next cell. Each cell is
assembled and factored exactly once per Newton update. All-cell A storage is
retained only by the independent diagnostic/reference API.
Diagnostic const-input condensation copies its input to preserve reference data;
this is not the production path.

Partial DG caches both directional volume flux Jacobians and HLL face
Jacobians. Block-Jacobi volume assembly reuses these caches instead of
re-evaluating the state and Jacobian. Volume/face evaluation uses a precomputed
mapping to native MFEM node order and sparse GL restriction/gradient tables:
volume values are direct nodal reads, volume gradient work is O((k+1)^3), and
face evaluation uses only its k+1 nonzero basis functions per face point.
The DG matrix-free volume action explicitly evaluates pointwise flux-Jacobian
products into two reusable cell buffers, then applies the one-dimensional
weak derivative contractions D_x^T and D_y^T on tensor lines. A coordinate-based
map translates lexicographic tensor indices to MFEM's native DOF/quadrature
ordering; no assumption about MFEM's ordering is needed. Two shared q-by-q
derivative tables suffice for this action, where q=k+1. The full shape tables
remain available for independent reference assembly and verification.
Mass inverse application in the collocated path is a direct diagonal scaling.
Weighted dot products, cell norms, and four-component diagonal mass scaling
carry explicit OpenMP-SIMD directives. IntelLLVM uses `-qopenmp-simd` with
AVX-512-capable `-march=native`, without enabling threaded OpenMP.
The partial DG diagonal builder applies these inverse-mass row weights while
accumulating sparse face terms directly into the cell block, avoiding two
full dense face-scratch clears/scales for every face.

HDG B and D depend on frozen wave speeds and dt. Bounds now follow each
nonlinear residual state, so a freeze-generation counter invalidates B and D
after every bound refresh; A and C refresh at every quasi-Newton state.
Nodal trace collocation lets A face terms accumulate directly into the cell
matrix. B/C store only four-component flux blocks per face GL node and borrow
the shared normal-direction interpolation tables. C first restricts along a
tensor line, then applies the flux block; B applies its flux block and lifts
along the same line with diagonal mass scaling. Thus B/C storage per cell is
O(k+1), and their vector actions are O((k+1)^2), rather than dense
O((k+1)^3) work. B is expanded directly into the cached X buffer only for
the multi-RHS local solve; no persistent dense B or C is stored.
This counter is valid because production changes dt
only alongside Freeze, before the HDG refresh at that nonlinear state.
Public callers changing dt themselves must Freeze again before reuse.

Schur Jacobi extraction applies only the matching face restriction C_s to the
matching X_s panel, forming four t-by-t products instead of a 4t-by-4t product.
Off-diagonal face products are not evaluated for Jacobi. D is stored as a vector
and applied pointwise, with no dense face D matrices. Partial Schur action
stores no global collection of C*X products; explicit element/full modes do.
Recovery reuses cached X/y and scratch, with no global recovery matrix.
The exact local X=A^{-1}B is generally dense: variable Euler flux Jacobians and
face stabilization do not give a separable Kronecker inverse. We retain X and
exact pivoted local LU. Tensor acceleration applies to DG derivatives and
HDG B/C restrictions/lifts, not an approximate replacement of A^{-1}.

## Krylov and ownership

Weighted restarted GMRES uses twice classical Gram-Schmidt and Givens rotations.
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
Tensor couplings also borrow Problem's interpolation and inverse-mass tables.
Replacing these vectors or retaining an action after its workspace dies is invalid.

Reusable scratch and mutable physical-Jacobian temporaries make a Problem/
workspace single-threaded and non-reentrant. A future OpenMP implementation
needs per-thread scratch and race-free face accumulation, not simply a pragma
around the current face loop. No MPI support is included.

## What still allocates and why

The optional explicit DG element/full paths rebuild their local operators;
explicit CSR and direct solvers still allocate their own storage. Partial DG
CSR remains a two-pass streamed reference/full-action path. Production
matrix-free ILU does not use this CSR: it fills cell/face BSR blocks directly.
Before the time loop, the fixed mesh graph, block permutation, elimination
update positions, factor buffers, pivots and solve scratch are prepared once.
Newton assembly fills the same BSR storage, and numerical ILU overwrites it.
Diagonal LU is also in place, avoiding a second diagonal copy. The only dropped
products are those outside the fixed block graph (ILU(0)). All numerical refresh
and factorization remain in Newton timing. Full-action GMRES still needs its
separate matrix for multiplication because PC storage is overwritten by ILU.
`--check-update` and verification deliberately allocate independent reference
objects; use them for correctness, not production benchmarking.

For maximum serial performance, start with Release MFEM/application and
`--assembly partial --precond block-jacobi --timing-only`. Whether ILU's reduced
Krylov count compensates for its assembly/setup cost is a measured question.
The default ordering `none` preserves the original cell/face ordering. `mdf`
now runs MFEM's MDF routine once on the unweighted topology graph, before the
time loop; it is no longer numerical MDF recomputed from each Newton matrix.
This uses the MDF symbol supplied by the pinned MFEM build (declared locally
because MFEM does not expose it in its public header). A future MFEM upgrade
must preserve that symbol or replace the ordering adapter. Neither ordering
depends on the Newton state. Native MFEM ILU is retained as a verification
reference, not as the production numerical-refresh path.

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
3. Cache the optional full/direct CSR pattern as well. Production ILU already
   uses a cached BSR graph and elimination schedule.
4. Tensor/operator-specific local factorizations or approximate local solves
   could reduce high-k cost. An approximate A^{-1} would change exact condensation
   and update equivalence unless used only as a PC; it is **not enabled**.
5. Overintegration, LGL nodes, reduced face traces, modified stabilization,
   mixed precision, fast-math, and looser tolerances can change the numerical
   method or equivalence guarantees. They are **not enabled** for speed claims.

The current implementation has been compiled and exercised on the target Intel
environment. `RESULTS.md` records the exact commands and separates verified
facts from performance conclusions.
