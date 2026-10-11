# Machine preset: grammar (KIAS), CPU cluster, OpenHPC modules, Slurm.
# Used through tests/machine.sh; the names it defines are listed in
# machines/macos-homebrew.sh. Run builds and runs inside a Slurm job
# (sbatch or salloc), not on the login node.
#
# `module` is missing in ssh commands and Slurm jobs, so the Lmod init file is sourced
# directly. Its init reads FPATH unset, which `set -u` rejects, so nounset is off while
# Lmod runs. openmpi4 and fftw are visible only after gnu12 is loaded. HDF5 is not
# loaded: -hdf5 with -mpi needs a parallel HDF5, which has not been checked here.
_opts=$-
set +u
if ! type module > /dev/null 2>&1; then
  export MODULEPATH=${MODULEPATH:-/opt/ohpc/pub/modulefiles}
  . /opt/ohpc/admin/lmod/lmod/init/bash > /dev/null
fi
module purge > /dev/null 2>&1
module load gnu12/12.2.0 openmpi4/4.1.5 fftw/3.3.10
case $_opts in *u*) set -u ;; esac
unset _opts

export CXX=g++
export OMPI_CXX=$CXX
TIGRIS_LIB_FLAGS=(--fftw_path="${FFTW_DIR}")
SUNDIALS_DIR=
PYTHON=$HOME/miniconda3/envs/pyathena/bin/python
# The login node allows at most 4 make jobs; inside a job use the allocated cores.
MAKE_JOBS=${SLURM_CPUS_ON_NODE:-4}

# openmpi4/4.1.5 has no Slurm PMI support, so srun aborts every rank before MPI_Init;
# mpirun inside the allocation works.
mpi_launch() {
  local n=$1
  shift
  mpirun -np "$n" "$@"
}
