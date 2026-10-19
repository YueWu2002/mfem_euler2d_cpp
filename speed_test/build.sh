#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SPEED_ROOT="${PROJECT_ROOT}/speed_test"
MFEM_BUILD=${MFEM_BUILD:-/oscar/home/ywu355/MFEM/build-intel-mpi}

module purge
module load intel-oneapi-compilers intel-oneapi-mkl hpcx-mpi cmake
export OMPI_CC=icx
export OMPI_CXX=icpx

cmake --fresh -S "${PROJECT_ROOT}" -B "${SPEED_ROOT}/build" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DMFEM_DIR="${MFEM_BUILD}" \
  -DEULER_NATIVE=ON \
  -DBUILD_TESTING=OFF
cmake --build "${SPEED_ROOT}/build" -j "${SLURM_CPUS_PER_TASK:-8}"
