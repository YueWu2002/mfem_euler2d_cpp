### CPU time per Newton step

Configuration: 256x256, CFL=5.0, SSP-DIRK3, four stages, three fixed
Newton updates per stage (12 measured updates), matrix-free partial
assembly, restarted weighted right GMRES (restart dimension 20,
maximum 1500 Krylov iterations), and one process/thread.

Peak RSS is measured per solver process over its complete lifetime,
including setup; missing measurements are shown as —, not estimates.

### Memory capacity observations

- 2026-10-19: 58.99 GB was sufficient to run the `k=8`, `256x256`,
  SSP-DIRK(4,3), block-Jacobi case. This is a demonstrated sufficient
  allocation, not a measured peak-RSS value or a minimum-memory bound.

| k | preconditioner | CFL | DG CPU time per Newton step (s) | HDG CPU time per Newton step (s) | DG/HDG | DG peak RSS (GiB) | HDG peak RSS (GiB) | node | job_id |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| — | — | — | — | — | — | — | — | — | — |

### Measured process peak RSS

Failed processes are retained here; their RSS is not a completed-run estimate.

| k | preconditioner | form | status | peak RSS (GiB) | wall (s) | node | job ID |
|---:|---|---|---|---:|---:|---|---:|
| — | — | — | no measurements | — | — | — | — |
