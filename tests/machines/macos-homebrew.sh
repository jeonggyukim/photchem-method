# Machine preset: macOS with Homebrew (gcc, open-mpi, hdf5-mpi, fftw, boost) and
# miniforge. Used through tests/machine.sh; see "Machines" in tests/README.md.
#
# Every preset defines the same names:
#   CXX, OMPI_CXX     C++ compiler for the harnesses and the one mpicxx calls
#   TIGRIS_LIB_FLAGS  library flags appended to every Tigris configure.py
#   CMAKE             the cmake that configures AthenaK
#   ATHENAK_CMAKE_FLAGS  device and compiler flags appended to every AthenaK cmake
#   SUNDIALS_DIR      SUNDIALS prefix (ode_solver/bdf_reference_check.cpp only)
#   PYTHON            the Python that has this repository installed
#   MAKE_JOBS         parallel make jobs
#   mpi_launch N CMD  how to start CMD on N MPI ranks

# Homebrew's open-mpi defaults to clang++; GCC builds need OMPI_CXX set.
export CXX=/opt/homebrew/bin/g++-16
export OMPI_CXX=$CXX

# hdf5-mpi, fftw and boost are keg-only, so their prefixes are given explicitly.
TIGRIS_LIB_FLAGS=(--hdf5_path=/opt/homebrew/opt/hdf5-mpi
                  --include=/opt/homebrew/opt/boost/include
                  --fftw_path=/opt/homebrew/opt/fftw)
CMAKE=cmake
# AthenaK on the CPU, serial Kokkos backend, with the same GCC (its objects are
# libstdc++; clang would link libc++).
ATHENAK_CMAKE_FLAGS=(-D CMAKE_CXX_COMPILER="$CXX")
SUNDIALS_DIR=/opt/homebrew/opt/sundials
PYTHON=/opt/homebrew/Caskroom/miniforge/base/envs/pyathena/bin/python
MAKE_JOBS=8

mpi_launch() {
  local n=$1
  shift
  mpirun -np "$n" "$@"
}
