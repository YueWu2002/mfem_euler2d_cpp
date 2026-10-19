#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "${PROJECT_ROOT}"
mkdir -p speed_test/slurm speed_test/results/be_newton3

MANIFEST="speed_test/results/be_newton3/submitted_jobs_$(date +%Y%m%dT%H%M%S)_$$.csv"
printf 'sequence,job_id,k,preconditioner\n' > "${MANIFEST}"

# Submission order only; scheduler start order is not guaranteed.
degrees=(8 7 6 5 4 3 2 1 0)
preconditioners=(block-jacobi block-ilu0)
sequence=0
for k in "${degrees[@]}"; do
  for preconditioner in "${preconditioners[@]}"; do
    sequence=$((sequence + 1))
    job_id=$(sbatch --parsable \
      --job-name="euler-k${k}-${preconditioner}" \
      --export="ALL,TEST_K=${k},TEST_PRECOND=${preconditioner}" \
      speed_test/job.sh)
    job_id=${job_id%%;*}
    printf '%d,%s,%d,%s\n' \
      "${sequence}" "${job_id}" "${k}" "${preconditioner}" \
      >> "${MANIFEST}"
  done
done

echo "Submitted ${sequence} jobs; manifest: ${MANIFEST}"
