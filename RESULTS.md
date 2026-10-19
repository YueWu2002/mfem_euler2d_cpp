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
   Compare isolated single-thread mean accepted Newton-update time and
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

HLL bounds are now fixed throughout each implicit solve: one backward-Euler
step or one DIRK stage, including Newton iterations and line-search trials.
Existing measurements below predate this policy and need rerunning to validate
the modified method. The code-test
scope is intentionally limited to a smoke run and a backward-Euler one-step
Newton-converged DG/HDG update-equivalence check:

```sh
python code_test/test_solver.py --exe build-intel/euler2d
```

Both cases select dt from CFL; no code test uses a fixed short dt.

All future tests use Newton and GMRES tolerances 1e-12. Existing results below
retain their actual run settings; they have not been rerun with the new defaults.
Fixed-update speed tests retain their prescribed update counts.

The next accuracy sweep uses the same k=2, SSP-DIRK(4,3), CFL=2.0,
T=0.2 and N=32,64,128,256,512 settings, with both tolerances changed to
1e-12 and at most 30 Newton updates per stage. No results at this new
configuration have been collected yet.

## Collected accuracy sweep: k=2, SSP-DIRK(4,3), tolerance 1e-13

Status: job 7584466 exited with failure after completing N=32,64,128.
N=256 and N=512 both reported `Newton failed to converge` during DG;
neither has a final DG/HDG summary or comparison. No full-sweep
`convergence.csv` was produced. Raw output is retained in
`accuracy_test/results/job_7584466/`; missing final results are not inferred.

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

Job 7584466 used the following command (the current submission script now uses
1e-12 for both tolerances):

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
| 32 | DG | 1.00539387e-3 (-) | 7.41333053e-4 (-) | 7.41333053e-4 (-) | 6.88933246e-4 (-) |
| 32 | HDG | 1.00539387e-3 (-) | 7.41333053e-4 (-) | 7.41333053e-4 (-) | 6.88933246e-4 (-) |
| 64 | DG | 1.51106794e-4 (2.73412) | 1.12145056e-4 (2.72476) | 1.12145056e-4 (2.72476) | 9.26221620e-5 (2.89493) |
| 64 | HDG | 1.51106794e-4 (2.73412) | 1.12145056e-4 (2.72476) | 1.12145056e-4 (2.72476) | 9.26221620e-5 (2.89493) |
| 128 | DG | 2.08487965e-5 (2.85753) | 1.51791945e-5 (2.88520) | 1.51791945e-5 (2.88520) | 1.40643951e-5 (2.71931) |
| 128 | HDG | 2.08487965e-5 (2.85753) | 1.51791945e-5 (2.88520) | 1.51791945e-5 (2.88520) | 1.40643951e-5 (2.71931) |
| 256 | DG | failed | failed | failed | failed |
| 256 | HDG | not reached | not reached | not reached | not reached |
| 512 | DG | failed | failed | failed | failed |
| 512 | HDG | not reached | not reached | not reached | not reached |

Observed convergence orders for the collected 1e-13 sweep:

For each component, p = log(e_coarse/e_fine) / log(N_fine/N_coarse),
using final-time L2 errors at fixed CFL=2.0 and T=0.2. These are
CFL-coupled space-time orders. DG and HDG agree to the five digits shown.

| refinement | form | rho order | rhou order | rhov order | E order |
|:---|:---:|---:|---:|---:|---:|
| 32 -> 64 | DG | 2.73412 | 2.72476 | 2.72476 | 2.89493 |
| 32 -> 64 | HDG | 2.73412 | 2.72476 | 2.72476 | 2.89493 |
| 64 -> 128 | DG | 2.85753 | 2.88520 | 2.88520 | 2.71931 |
| 64 -> 128 | HDG | 2.85753 | 2.88520 | 2.88520 | 2.71931 |

Orders involving N=256 or N=512 are unavailable because those runs failed
before producing final solutions.

Final DG-minus-HDG solution differences:

| N | rho L2 | rhou L2 | rhov L2 | E L2 | joint solution L2 | max nodal |
|---:|---:|---:|---:|---:|---:|---:|
| 32 | 3.99344942e-15 | 3.83026092e-15 | 3.83886988e-15 | 3.89280857e-15 | 7.77877995e-15 | 1.04916076e-14 |
| 64 | 5.44246970e-15 | 5.32318467e-15 | 5.30668096e-15 | 5.51713163e-15 | 1.07961278e-14 | 1.30451205e-14 |
| 128 | 6.69837806e-15 | 6.59156515e-15 | 6.54029160e-15 | 6.84283216e-15 | 1.33385444e-14 | 1.90958360e-14 |
| 256 | unavailable | unavailable | unavailable | unavailable | unavailable | unavailable |
| 512 | unavailable | unavailable | unavailable | unavailable | unavailable | unavailable |

Total accepted Newton updates over all DIRK stages:

| N | DG | HDG | absolute difference |
|---:|---:|---:|---:|
| 32 | 52 | 52 | 0 |
| 64 | 78 | 78 | 0 |
| 128 | 129 | 129 | 0 |
| 256 | unavailable | unavailable | unavailable |
| 512 | unavailable | unavailable | unavailable |

For the completed grids, requested and actual maximum CFL are both 2.0.
Time-step counts are 3, 5, 10 and mean dt values are 0.0666666667,
0.04, 0.02, respectively; all completed summaries report T=0.2.
The maximum recorded conservation drift is 1.20792265e-13.
DG/HDG total accepted Newton updates match exactly for these grids.
The incomplete sweep remains in its original job directory rather than being
archived as a complete experiment in `convergence_order/`.

## Backward-Euler Newton timing workload

New speed tests use one backward-Euler step and exactly three Newton updates
for each of DG and HDG, with matching CFL and solver settings. No measurements
for this workload have been collected yet. The DIRK results below remain historical.
Raw data: `speed_test/results/be_newton3/job_<ID>/k<k>_<precond>/`.
Collect with `python speed_test/collect.py --update-results`.

<!-- speed-test-be-results:start -->
Pending: backward Euler, 1 time step, 3 Newton updates per method.
<!-- speed-test-be-results:end -->

## Retained historical paired performance measurements

<!-- speed-test-results:start -->
### CPU time per Newton step

Collected 2026-10-19: 16 of 18 cases completed. Block-ILU(0) at k=7 and k=8
failed during DG with OUT_OF_MEMORY at a 120 GiB allocation; HDG was not run.
Failed-process RSS below is not a successful-run memory requirement.

CPU: Intel(R) Xeon(R) Gold 6242 CPU @ 2.80GHz (node1601), one thread.
Despite the historical table title, mean_newton_s measures elapsed wall time
per accepted Newton update, not process CPU time. Process CPU and wall times
are recorded separately in MEMORY_ESTIMATE.md and each case's resources.csv.

Configuration: 256x256, CFL=5.0, SSP-DIRK3, four stages, three fixed
Newton updates per stage (12 measured updates), matrix-free partial
assembly, restarted weighted right GMRES (restart dimension 20, maximum 1500 Krylov iterations), and one process/thread.

Peak RSS is measured per solver process over its complete lifetime,
including setup; missing measurements are shown as —, not estimates.

| k | preconditioner | CFL | DG CPU time per Newton step (s) | HDG CPU time per Newton step (s) | DG/HDG | DG peak RSS (GiB) | HDG peak RSS (GiB) | node | job_id |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | block-ilu0 | 5.0 | 0.63753157566396523 | 2.0773129234904943 | 0.30690204083105854 | 0.2716484069824219 | 0.569091796875 | node1601.oscar.ccv.brown.edu | 7583445 |
| 0 | block-jacobi | 5.0 | 1.2944837657507728 | 3.4374422644117448 | 0.3765834205137697 | 0.2038726806640625 | 0.3782844543457031 | node1601.oscar.ccv.brown.edu | 7583444 |
| 1 | block-ilu0 | 5.0 | 3.5944653986586501 | 5.5809946404750947 | 0.6440546229144315 | 1.2134208679199219 | 1.4156227111816406 | node1601.oscar.ccv.brown.edu | 7583443 |
| 1 | block-jacobi | 5.0 | 6.0816388836003581 | 7.5270784365879848 | 0.8079680496005509 | 0.6711692810058594 | 0.9375953674316406 | node1601.oscar.ccv.brown.edu | 7583442 |
| 2 | block-ilu0 | 5.0 | 15.661668416753772 | 14.593375791411079 | 1.0732039413369616 | 4.321250915527344 | 2.917133331298828 | node1601.oscar.ccv.brown.edu | 7583441 |
| 2 | block-jacobi | 5.0 | 27.856269659506665 | 24.902720978916022 | 1.118603452333232 | 1.7382240295410156 | 1.9668388366699219 | node1601.oscar.ccv.brown.edu | 7583440 |
| 3 | block-ilu0 | 5.0 | 40.146268260489116 | 24.858708220243546 | 1.614978055367987 | 11.949943542480469 | 5.264698028564453 | node1601.oscar.ccv.brown.edu | 7583439 |
| 3 | block-jacobi | 5.0 | 59.157091871153646 | 38.210232769404904 | 1.5482002485606678 | 3.8802223205566406 | 3.654285430908203 | node1601.oscar.ccv.brown.edu | 7583438 |
| 4 | block-ilu0 | 5.0 | 109.98157531642967 | 45.723570625073414 | 2.4053584138969044 | 27.354736328125 | 8.634872436523438 | node1601.oscar.ccv.brown.edu | 7583437 |
| 4 | block-jacobi | 5.0 | 144.84621904102275 | 82.978139214334121 | 1.7455949291280466 | 7.740577697753906 | 6.175926208496094 | node1601.oscar.ccv.brown.edu | 7583436 |
| 5 | block-ilu0 | 5.0 | 218.36511683250339 | 72.284765028669895 | 3.0209009705695866 | 54.7783203125 | 13.223461151123047 | node1601.oscar.ccv.brown.edu | 7583435 |
| 5 | block-jacobi | 5.0 | 251.45108258575462 | 125.85146847617034 | 1.9979987967591037 | 14.17373275756836 | 9.731624603271484 | node1601.oscar.ccv.brown.edu | 7583434 |
| 6 | block-ilu0 | 5.0 | 445.10608194175228 | 116.70809940740641 | 3.8138405492147482 | 99.36463165283203 | 19.210681915283203 | node1601.oscar.ccv.brown.edu | 7583433 |
| 6 | block-jacobi | 5.0 | 452.2425671913295 | 212.96592885597963 | 2.123544219588116 | 24.203170776367188 | 14.49472427368164 | node1601.oscar.ccv.brown.edu | 7583432 |
| 7 | block-jacobi | 5.0 | 794.52435649201561 | 353.14934146215506 | 2.2498253945552382 | 39.05031967163086 | 20.670978546142578 | node1601.oscar.ccv.brown.edu | 7583430 |
| 8 | block-jacobi | 5.0 | 1375.5149045302605 | 474.81393736684703 | 2.8969556204655404 | 60.13582992553711 | 28.433704376220703 | node1601.oscar.ccv.brown.edu | 7583428 |

### Measured process peak RSS

Failed processes are retained here; their RSS is not a completed-run estimate.

| k | preconditioner | form | status | peak RSS (GiB) | wall (s) | node | job ID |
|---:|---|---|---|---:|---:|---|---:|
| 0 | block-ilu0 | dg | completed | 0.2716484069824219 | 10.94 | node1601.oscar.ccv.brown.edu | 7583445 |
| 0 | block-ilu0 | hdg | completed | 0.569091796875 | 28.53 | node1601.oscar.ccv.brown.edu | 7583445 |
| 0 | block-jacobi | dg | completed | 0.2038726806640625 | 18.73 | node1601.oscar.ccv.brown.edu | 7583444 |
| 0 | block-jacobi | hdg | completed | 0.3782844543457031 | 44.61 | node1601.oscar.ccv.brown.edu | 7583444 |
| 1 | block-ilu0 | dg | completed | 1.2134208679199219 | 47.55 | node1601.oscar.ccv.brown.edu | 7583443 |
| 1 | block-ilu0 | hdg | completed | 1.4156227111816406 | 71.71 | node1601.oscar.ccv.brown.edu | 7583443 |
| 1 | block-jacobi | dg | completed | 0.6711692810058594 | 77.05 | node1601.oscar.ccv.brown.edu | 7583442 |
| 1 | block-jacobi | hdg | completed | 0.9375953674316406 | 94.65 | node1601.oscar.ccv.brown.edu | 7583442 |
| 2 | block-ilu0 | dg | completed | 4.321250915527344 | 194.77 | node1601.oscar.ccv.brown.edu | 7583441 |
| 2 | block-ilu0 | hdg | completed | 2.917133331298828 | 181.33 | node1601.oscar.ccv.brown.edu | 7583441 |
| 2 | block-jacobi | dg | completed | 1.7382240295410156 | 339.41 | node1601.oscar.ccv.brown.edu | 7583440 |
| 2 | block-jacobi | hdg | completed | 1.9668388366699219 | 304.31 | node1601.oscar.ccv.brown.edu | 7583440 |
| 3 | block-ilu0 | dg | completed | 11.949943542480469 | 494.4 | node1601.oscar.ccv.brown.edu | 7583439 |
| 3 | block-ilu0 | hdg | completed | 5.264698028564453 | 306.78 | node1601.oscar.ccv.brown.edu | 7583439 |
| 3 | block-jacobi | dg | completed | 3.8802223205566406 | 716.92 | node1601.oscar.ccv.brown.edu | 7583438 |
| 3 | block-jacobi | hdg | completed | 3.654285430908203 | 465.88 | node1601.oscar.ccv.brown.edu | 7583438 |
| 4 | block-ilu0 | dg | completed | 27.354736328125 | 1342.49 | node1601.oscar.ccv.brown.edu | 7583437 |
| 4 | block-ilu0 | hdg | completed | 8.634872436523438 | 560.01 | node1601.oscar.ccv.brown.edu | 7583437 |
| 4 | block-jacobi | dg | completed | 7.740577697753906 | 1748.07 | node1601.oscar.ccv.brown.edu | 7583436 |
| 4 | block-jacobi | hdg | completed | 6.175926208496094 | 1005.32 | node1601.oscar.ccv.brown.edu | 7583436 |
| 5 | block-ilu0 | dg | completed | 54.7783203125 | 2660.91 | node1601.oscar.ccv.brown.edu | 7583435 |
| 5 | block-ilu0 | hdg | completed | 13.223461151123047 | 882.3 | node1601.oscar.ccv.brown.edu | 7583435 |
| 5 | block-jacobi | dg | completed | 14.17373275756836 | 3031.77 | node1601.oscar.ccv.brown.edu | 7583434 |
| 5 | block-jacobi | hdg | completed | 9.731624603271484 | 1522.74 | node1601.oscar.ccv.brown.edu | 7583434 |
| 6 | block-ilu0 | dg | completed | 99.36463165283203 | 5409.98 | node1601.oscar.ccv.brown.edu | 7583433 |
| 6 | block-ilu0 | hdg | completed | 19.210681915283203 | 1419.94 | node1601.oscar.ccv.brown.edu | 7583433 |
| 6 | block-jacobi | dg | completed | 24.203170776367188 | 5448.5 | node1601.oscar.ccv.brown.edu | 7583432 |
| 6 | block-jacobi | hdg | completed | 14.49472427368164 | 2572.16 | node1601.oscar.ccv.brown.edu | 7583432 |
| 7 | block-ilu0 | dg | failed | 119.74918746948242 | 87.51 | node1601.oscar.ccv.brown.edu | 7583431 |
| 7 | block-jacobi | dg | completed | 39.05031967163086 | 9566.69 | node1601.oscar.ccv.brown.edu | 7583430 |
| 7 | block-jacobi | hdg | completed | 20.670978546142578 | 4259.16 | node1601.oscar.ccv.brown.edu | 7583430 |
| 8 | block-ilu0 | dg | failed | 119.7475700378418 | 89.5 | node1601.oscar.ccv.brown.edu | 7583429 |
| 8 | block-jacobi | dg | completed | 60.13582992553711 | 16552.07 | node1601.oscar.ccv.brown.edu | 7583428 |
| 8 | block-jacobi | hdg | completed | 28.433704376220703 | 5723.66 | node1601.oscar.ccv.brown.edu | 7583428 |
<!-- speed-test-results:end -->

OSCAR's published hardware table lists batch nodes from 190 GB through 1540 GB
per node. Thus 1540 GB is the largest advertised batch-node capacity, not a
uniform per-job entitlement: an actual request is limited by the selected
node, account/QoS policy and availability. Slurm `--mem` is requested per node.
The current speed script pins node2321 and requests 64G; Slurm control queries
were unavailable while this note was updated, so this document does not claim
an unverified node2321 capacity or account-specific maximum.

These measurements have not been revalidated with the current per-implicit-solve
fixed-HLL policy and must be repeated
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
