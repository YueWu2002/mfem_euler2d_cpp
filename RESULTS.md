# DG/HDG correctness, convergence, and performance

## Purpose and acceptance criteria

The project has two linked objectives:

1. Establish that condensed HDG is an algebraically equivalent implementation
   of DG: compare recovered body-DOF updates, final body solutions, conservation,
   and convergence under identical discretization and solver tolerances.
   Equality means agreement within linear/nonlinear tolerances and roundoff,
   not bitwise identity. HDG solves the condensed face system and back-substitutes;
   DG solves the body system.
2. Demonstrate and quantify HDG acceleration where the measurements support it.
   Compare isolated single-thread mean accepted quasi-Newton-update time and
   total solution time at equal accuracy. Include Jacobian construction,
   condensation, preconditioner setup, global solve, recovery, and residual/
   line-search work. Report both speedups and slowdowns, with configurations.

Both methods use the same DIRK, HLL bounds, quadrature, time-step policy,
DG residual norm, and stopping criteria. PARDISO is only a small update
reference; production runs use weighted right-GMRES with block Jacobi or
block ILU(0). Convergence studies validate accuracy, not performance.

## Current algorithm and regression checks

The last-step flag is set after selecting dt by testing `time+dt >= T`;
when true, dt becomes `T-time` and the loop stops after that step.
`T=0` performs no steps; explicit `--steps` overrides T.

HLL bounds follow every nonlinear residual state, including line-search
trials, and are frozen only in each quasi-Newton linear solve. The code-test
scope is intentionally limited to a smoke run and a backward-Euler one-step
Newton-converged DG/HDG update-equivalence check:

```sh
python code_test/test_solver.py --exe build-intel/euler2d
```

Both cases select dt from CFL; no code test uses a fixed short dt.

## Current accuracy sweep: k=2, SSP-DIRK(4,3), tolerance 1e-13

Status: configured for a fresh run; no results for this configuration have been
collected yet. All five grids must be rerun with the updated code.

| category | setting |
|---|---|
| space | k=2 (Q2 GL body nodes, P2 GL traces) |
| time integration | SSP-DIRK(4,3): four stages, third order (`--integrator ssp-dirk3`) |
| final time | T=0.2; final step clipped to T |
| grids | N=32,64,128,256,512; NxN cells for both DG and HDG |
| CFL | 2.0 |
| Newton relative tolerance | 1e-13; maximum 25 updates per stage |
| GMRES relative tolerance | 1e-13; restart 20; maximum 1500 iterations |
| operator / preconditioner | partial assembly, matrix-free weighted right-GMRES, block Jacobi |
| execution | one Intel node, one CPU, single-threaded; cases run sequentially |
| Slurm resources | batch partition, constraint intel, 60G memory, 1 day including build; no node pinned |
| build and results | fresh build in `accuracy_test/results/job_JOBID/build/`; results in its parent directory |

`accuracy_test/job.sh` invokes:

```sh
python -u accuracy_test/convergence.py --exe "${BUILD_ROOT}/euler2d" \
  --grids 32 64 128 256 512 --T 0.2 --cfl 2.0 \
  --newton-tol 1e-13 --linear-tol 1e-13 --jobs 1 \
  --output "${RUN_ROOT}"
```

Data collection retains per-update iteration CSVs, per-form summaries,
four-component exact-solution L2 errors and convergence rates, all four final
DG-minus-HDG component L2 differences, joint solution L2 and maximum nodal
difference, total accepted Newton updates across all stages, conservation drift,
and timing/log files. The existing physical residual norms and roundoff safeguards
remain in effect; 1e-13 specifies the requested solver tolerances.

Final-time exact-solution errors and rates (`L2 error (rate)`):

| N | form | rho | rhou | rhov | E |
|---:|:---:|---:|---:|---:|---:|
| 32 | DG | pending | pending | pending | pending |
| 32 | HDG | pending | pending | pending | pending |
| 64 | DG | pending | pending | pending | pending |
| 64 | HDG | pending | pending | pending | pending |
| 128 | DG | pending | pending | pending | pending |
| 128 | HDG | pending | pending | pending | pending |
| 256 | DG | pending | pending | pending | pending |
| 256 | HDG | pending | pending | pending | pending |
| 512 | DG | pending | pending | pending | pending |
| 512 | HDG | pending | pending | pending | pending |

Final DG-minus-HDG solution differences:

| N | rho L2 | rhou L2 | rhov L2 | E L2 | joint solution L2 | max nodal |
|---:|---:|---:|---:|---:|---:|---:|
| 32 | pending | pending | pending | pending | pending | pending |
| 64 | pending | pending | pending | pending | pending | pending |
| 128 | pending | pending | pending | pending | pending | pending |
| 256 | pending | pending | pending | pending | pending | pending |
| 512 | pending | pending | pending | pending | pending | pending |

Total accepted Newton updates over all DIRK stages:

| N | DG | HDG | absolute difference |
|---:|---:|---:|---:|
| 32 | pending | pending | pending |
| 64 | pending | pending | pending |
| 128 | pending | pending | pending |
| 256 | pending | pending | pending |
| 512 | pending | pending | pending |

## Previous Q2 DIRK results: tolerance 1e-12

The following records use the previous code and tolerance 1e-12. They are not
results of the current 1e-13 sweep above.

```sh
python accuracy_test/convergence.py --exe build-intel/euler2d \
  --grids 32 64 128 256 --T .2 --cfl 2.0 \
  --newton-tol 1e-12 --linear-tol 1e-12 --jobs 4 \
  --output accuracy_test/convergence_order/dirk3_cfl2_tol1e12_T02
```

All parameters for this production convergence run are fixed as follows:

| category | setting |
|---|---|
| PDE/case | 2D Euler, gamma=3, alpha=pi/4, periodic rectangle, T=0.2 |
| space | NxN affine cells, Q2 GL body nodes, P2 GL trace nodes |
| quadrature | 3x3 tensor GL volume rule; 3-point GL face rule |
| numerical flux | Einfeldt/Roe HLL, safety coefficient 1.1 |
| HLL update | recomputed at every nonlinear residual/trial state; frozen only in each quasi-Newton linear solve |
| time integration | four-stage third-order DIRK; final step clipped exactly to T |
| time-step selection | CFL=2.0 using safety-scaled maximum HLL bound and min(hx,hy) |
| nonlinear stopping | DG physical weak-residual dual norm; relative tolerance 1e-12; maximum 25 updates per stage |
| linear solver | matrix-free weighted right-preconditioned GMRES |
| GMRES metric | physical dual norm described below |
| right preconditioner | block Jacobi; body-cell blocks for DG, trace-face blocks for HDG |
| GMRES controls | relative tolerance 1e-12, restart 20, maximum 1500 iterations |
| assembly mode | partial/operator action; no global Jacobian or Schur matrix for Mult |
| forms | DG body system and condensed HDG trace system with body recovery |
| execution | one thread per case; OMP/MKL/OpenBLAS/BLIS/vecLib thread limits all 1 |
| concurrency | up to 4 independent cases; concurrent wall times are not performance benchmarks |
| output | no mesh/solution fields; per-update CSV, per-case summary, aggregate convergence CSV |

The stored DG residual is `r_s=M^{-1}r_w`, where `r_w` is the physical weak
residual and M is the physical block mass matrix. `CellNorm` computes

```text
sqrt(r_s^T M r_s) = sqrt(r_w^T M^{-1} r_w),
```

which is exactly the requested DG dual norm. Weighted GMRES uses the same
metric for the DG system. The condensed HDG solve uses its corresponding trace
metric, then validates the recovered correction with the common DG dual norm;
Newton stopping and line search for both forms always use that DG dual norm.
This is a weighted right-preconditioned formulation, algebraically equivalent
to applying the paired Gram transform before a conventional preconditioned
iteration; it is not an unweighted Euclidean residual test.

The tables below are filled only from complete CFL=2.0, tolerance=1e-12 runs:

Final-time exact-solution errors and rates for DG and HDG
(each entry will be `L2 error (rate)`):

| N | form | rho | rhou | rhov | E |
|---:|:---:|---:|---:|---:|---:|
| 32 | DG | 1.00539387e-3 (-) | 7.41333053e-4 (-) | 7.41333053e-4 (-) | 6.88933246e-4 (-) |
| 32 | HDG | 1.00539387e-3 (-) | 7.41333053e-4 (-) | 7.41333053e-4 (-) | 6.88933246e-4 (-) |
| 64 | DG | 1.51106794e-4 (2.73412) | 1.12145056e-4 (2.72476) | 1.12145056e-4 (2.72476) | 9.26221620e-5 (2.89493) |
| 64 | HDG | 1.51106794e-4 (2.73412) | 1.12145056e-4 (2.72476) | 1.12145056e-4 (2.72476) | 9.26221620e-5 (2.89493) |
| 128 | DG | pending | pending | pending | pending |
| 128 | HDG | pending | pending | pending | pending |
| 256 | DG | pending | pending | pending | pending |
| 256 | HDG | pending | pending | pending | pending |

Final DG-minus-HDG solution differences:

| N | rho L2 | rhou L2 | rhov L2 | E L2 | joint solution L2 | max nodal |
|---:|---:|---:|---:|---:|---:|---:|
| 32 | 3.79099009e-15 | 3.62107223e-15 | 3.60769294e-15 | 3.69794418e-15 | 7.36029955e-15 | 8.77076189e-15 |
| 64 | 5.42860492e-15 | 5.29281312e-15 | 5.36151655e-15 | 5.54160835e-15 | 1.08138293e-14 | 1.45439216e-14 |
| 128 | pending | pending | pending | pending | pending | pending |
| 256 | pending | pending | pending | pending | pending | pending |

Total accepted Newton updates over all DIRK stages:

| N | DG | HDG | absolute difference |
|---:|---:|---:|---:|
| 32 | 49 | 49 | 0 |
| 64 | 72 | 72 | 0 |
| 128 | 121 | pending | pending |
| 256 | pending | pending | pending |

The completed paired runs require identical total Newton update counts. Pending
entries are not inferred from interrupted iteration files.

Here `DG-HDG solution L2` is the joint physical-domain norm
`sqrt(||delta rho||_L2^2 + ||delta rhou||_L2^2 + ||delta rhov||_L2^2 +
||delta E||_L2^2)`, evaluated with the same high-order rule used for exact
solution errors. The comparison CSV also records all four component L2 norms.

`max nodal` is the largest absolute DG-minus-HDG coefficient over every cell,
every local tensor GL body node, and all four conservative components. Since
the Q2 basis is GL nodal, these coefficients are the solution values at those
cell-local nodes; discontinuous copies on neighboring cells are each included.
It is an infinity norm, not an L2 norm.

Every component L2 entry reconstructs the cell Q2 polynomial and integrates
its squared error over the physical domain using MFEM quadrature order 24
(`max(24,2*k+8)`), then takes the square root. DG--HDG component differences
use the same reconstruction and rule. These are absolute norms without domain
or exact-solution normalization. See MATHEMATICS.md for the formulas and for
the distinction from the DG dual residual norm used by the solvers.

Matching exact-solution error norms alone do not prove identical solutions;
the separate body-update and final-solution comparisons provide that evidence.
Per-case iteration and summary CSVs are written as the sweep progresses;
the aggregate convergence CSV is written at completion.

## Retained historical paired performance measurements

<!-- speed-test-results:start -->
### CPU time per Newton step

Configuration: 256x256, CFL=5.0, SSP-DIRK3, four stages, three fixed
Newton updates per stage (12 measured updates), matrix-free partial
assembly, restarted weighted right GMRES, and one process/thread.

Peak RSS is measured per solver process over its complete lifetime,
including setup; missing measurements are shown as —, not estimates.

| k | preconditioner | CFL | DG CPU time per Newton step (s) | HDG CPU time per Newton step (s) | DG/HDG | DG peak RSS (GiB) | HDG peak RSS (GiB) | node | job_id |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| — | — | — | — | — | — | — | — | — | — |

### Measured process peak RSS

Failed processes are retained here; their RSS is not a completed-run estimate.

| k | preconditioner | form | status | peak RSS (GiB) | wall (s) | node | job ID |
|---:|---|---|---|---:|---:|---|---:|
| — | — | — | no measurements | — | — | — | — |
<!-- speed-test-results:end -->

OSCAR's published hardware table lists batch nodes from 190 GB through 1540 GB
per node. Thus 1540 GB is the largest advertised batch-node capacity, not a
uniform per-job entitlement: an actual request is limited by the selected
node, account/QoS policy and availability. Slurm `--mem` is requested per node.
The current speed script pins node2321 and requests 64G; Slurm control queries
were unavailable while this note was updated, so this document does not claim
an unverified node2321 capacity or account-specific maximum.

These measurements predate the stage-state HLL update and must be repeated
before claiming current-method speedups. They are retained because both
favorable and unfavorable results inform optimization; unrelated historical
convergence tables have been removed.

Environment: Intel Xeon Platinum 8268, IntelLLVM 2025.2, MFEM commit
`54d7f64513`, Release `-march=native`, MKL sequential, all BLAS/OpenMP
thread counts fixed to one. Workload: 16x16 Q2, CFL=.1, three time steps,
three updates per DIRK stage, warmup and alternating DG/HDG execution.
Ratios below are DG time / HDG time; a ratio greater than one favors HDG.

| PC | repetitions | DG mean update (s) | HDG mean update (s) | ratio | DG loop (s) | HDG loop (s) | ratio |
|---|---:|---:|---:|---:|---:|---:|---:|
| block Jacobi | 5 | .0148320 | .0276546 | .5363 | .536279 | .997925 | .5374 |
| block ILU(0) | 3 | .0784877 | .0488920 | 1.6053 | 2.82789 | 1.76253 | 1.6044 |

Raw data: `speed_test/results/historical/benchmark_k2/`,
`speed_test/results/historical/benchmark_k2_ilu0/`.
A fixed-update benchmark measures controlled work, not necessarily a
fully converged solution; converged total-solve timings must be reported
separately. Fewer trace DOFs alone are not evidence of acceleration.

## Data hygiene

Only CFL-coupled space-time convergence studies are retained. RK4, fixed-short-
dt, fixed-mesh time scans, and incomplete single-grid pilot data were deleted.
