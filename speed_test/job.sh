#!/usr/bin/env bash
#SBATCH --job-name=euler-speed
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodelist=node1601
#SBATCH --mem=200G
#SBATCH --time=10:00:00
#SBATCH --output=speed_test/slurm/%x-%j.out
#SBATCH --error=speed_test/slurm/%x-%j.err

set -euo pipefail
PROJECT_ROOT=${SLURM_SUBMIT_DIR}
cd "${PROJECT_ROOT}"

module purge
module load intel-oneapi-compilers intel-oneapi-mkl hpcx-mpi
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

test -n "${TEST_K:-}"
test -n "${TEST_PRECOND:-}"
CASE_ROOT="speed_test/results/be_newton3/job_${SLURM_JOB_ID}/k${TEST_K}_${TEST_PRECOND}"
srun --mpi=pmix --export=ALL --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
  python speed_test/run_case.py \
    --exe speed_test/build/euler2d \
    --k "${TEST_K}" --precond "${TEST_PRECOND}" --cfl 5.0 \
    --output "${CASE_ROOT}"
