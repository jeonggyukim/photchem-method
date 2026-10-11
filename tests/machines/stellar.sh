# Machine preset: stellar (Princeton), Intel compilers, Slurm.
# Used through tests/machine.sh; the names it defines are listed in
# machines/macos-homebrew.sh. Run builds and runs inside a Slurm allocation
# (salloc or sbatch): mpirun and srun fail on the login node.
#
# INCOMPLETE: the modules are not filled in. Load the Intel compiler, Intel MPI and a
# parallel HDF5 that sets HDF5DIR (find their names with `module avail`), then remove
# this note. The two settings below are the ones known to work for Tigris on stellar.
module purge
: "${HDF5DIR:?tests/machines/stellar.sh: add the module loads, see its first lines}"

# Tigris is configured with --cxx=icpx; the harnesses use the same compiler.
export CXX=icpx
TIGRIS_LIB_FLAGS=(--cxx=icpx --lib_path="${HDF5DIR}/lib64" --include="${HDF5DIR}/include")
CMAKE=cmake
ATHENAK_CMAKE_FLAGS=(-D CMAKE_CXX_COMPILER="$CXX")
SUNDIALS_DIR=
PYTHON=python3
MAKE_JOBS=${SLURM_CPUS_ON_NODE:-8}

mpi_launch() {
  local n=$1
  shift
  srun -n "$n" "$@"
}
