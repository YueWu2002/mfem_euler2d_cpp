#!/usr/bin/env bash
set -euo pipefail

# Serial application build against the Intel/MKL MFEM build documented in
# build-intel-mfem.sh.  Override MFEM_BUILD or BUILD_DIR when needed.
module purge
module load intel-oneapi-compilers intel-oneapi-mkl hpcx-mpi cmake

export OMPI_CC=icx
export OMPI_CXX=icpx
: "${MFEM_BUILD:=/oscar/home/ywu355/MFEM/build-intel-mpi}"
: "${BUILD_DIR:=build-intel}"

cmake --fresh -S . -B "${BUILD_DIR}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DMFEM_DIR="${MFEM_BUILD}" \
  -DEULER_NATIVE=ON
  # enable hardward-aligned optimization for the compilation host CPU
cmake --build "${BUILD_DIR}" -j
