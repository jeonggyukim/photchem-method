# Machine preset: syntax (KIAS), GPU cluster for AthenaK, same modules and $HOME as
# grammar. Used through tests/machine.sh; the names it defines are listed in
# machines/macos-homebrew.sh. GPU partitions reject CPU-only jobs: build Tigris and
# run CPU tests on grammar or on the a10 partition.
_opts=$-
set +u
if ! type module > /dev/null 2>&1; then
  export MODULEPATH=${MODULEPATH:-/opt/ohpc/pub/modulefiles}
  . /opt/ohpc/admin/lmod/lmod/init/bash > /dev/null
fi
module purge > /dev/null 2>&1
module load gnu12/12.2.0 openmpi4/4.1.5 fftw/3.3.10 cuda/12.6.3
case $_opts in *u*) set -u ;; esac
unset _opts

export CXX=g++
export OMPI_CXX=$CXX
TIGRIS_LIB_FLAGS=(--fftw_path="${FFTW_DIR}")
SUNDIALS_DIR=
PYTHON=$HOME/miniconda3/envs/pyathena/bin/python
MAKE_JOBS=${SLURM_CPUS_ON_NODE:-4}

# CMake 3.31 must come first on PATH for AthenaK's test suite.
export PATH=$HOME/libs/cmake/bin:$PATH

mpi_launch() {
  local n=$1
  shift
  mpirun -np "$n" "$@"
}
