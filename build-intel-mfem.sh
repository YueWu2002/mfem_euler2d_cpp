####### install MFEM on Intel nodes on Brown University CCV OSCAR #######

# must use Intel CPUs, e.g., use "-C intel" when submitting sbatch jobs or "-f intel" when starting interactive/salloc jobs

# file structure:
#/oscar/home/ywu355/MFEM/   # current directory!
#├── mfem/                  # MFEM source
#├── build-intel-mpi/       # Intel + HPC-X + MKL build
#├── install-intel-mpi/     # install
#├── build-gcc-mpi/         # optional: GCC version build
#└── install-gcc-mpi/       # optional: GCC version install

# make sure the current directory is MFEM

# download mfem package as the mfem directory
rm -rf mfem
git clone https://github.com/mfem/mfem.git

# load modules
module purge # clear loaded modules
module load intel-oneapi-compilers intel-oneapi-mkl hpcx-mpi # MPI for Brown OSCAR
module load cmake
# for WSL2: source /opt/intel/oneapi/setvars.sh

# file structure
export MFEM_ROOT=/oscar/home/ywu355/MFEM
export MFEM_SRC=$MFEM_ROOT/mfem
export MFEM_BUILD=$MFEM_ROOT/build-intel-mpi
export MFEM_INSTALL=$MFEM_ROOT/install-intel-mpi

rm -rf "$MFEM_BUILD" "$MFEM_INSTALL"
mkdir -p "$MFEM_BUILD" "$MFEM_INSTALL"

# make sure that MFEM_SRC exists!

# select Intel compilers for HPC-X MPI compilers and check
export OMPI_CC=icx
export OMPI_CXX=icpx

# cmake configuration
cmake --fresh -S "$MFEM_SRC" -B "$MFEM_BUILD" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_INSTALL_PREFIX="$MFEM_INSTALL" \
  -DMFEM_USE_MPI=ON \
  -DMFEM_USE_METIS=ON \
  -DMFEM_FETCH_TPLS=ON \
  -DMFEM_USE_LAPACK=ON \
  -DBLA_VENDOR=Intel10_64lp_seq \
  -DMFEM_USE_MKL_PARDISO=ON \
  -DMFEM_USE_MKL_CPARDISO=ON \
  -DMKL_PARDISO_DIR="$MKLROOT" \
  -DMKL_CPARDISO_DIR="$MKLROOT" \
  -DMKL_LIBRARY_DIR=lib \
  -DMKL_MPI_WRAPPER_LIB=mkl_blacs_openmpi_lp64 \
  -DMFEM_USE_SIMD=ON \
  -DMFEM_USE_OPENMP=OFF

# build MFEM
cmake --build "$MFEM_BUILD" -j

# install MFEM
cmake --install "$MFEM_BUILD"

# register the default testing command
cmake -S "$MFEM_SRC" -B "$MFEM_BUILD" \
  -DMFEM_MPIEXEC="srun --mpi=pmix" \
  -DMFEM_MPIEXEC_NP="-n" \
  -DMFEM_MPI_NP=2

# build tests
cmake --build "$MFEM_BUILD" --target tests -j       # Build unit-tests
cmake --build "$MFEM_BUILD" --target examples -j    # Build examples
cmake --build "$MFEM_BUILD" --target miniapps -j    # Build miniapps

# execute tests
ctest --test-dir "$MFEM_BUILD" --output-on-failure -T test

# test ex1p
cmake --build "$MFEM_BUILD" --target ex1p -j
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
# export all Shell exported env vars to srun; choose working directory
srun --mpi=pmix -n 2 --chdir="$MFEM_BUILD/examples" ./ex1p -no-vis # general case
srun --mpi=pmix --export=ALL -n 2 --chdir="$MFEM_BUILD/examples" ./ex1p -no-vis # when using Open OnDemand Code Server (since it sets SLURM_EXPORT_ENV to None)
# for WSL2: mpirun -np 2 --chdir="$MFEM_BUILD/examples" ./ex1p -no-vis # general case