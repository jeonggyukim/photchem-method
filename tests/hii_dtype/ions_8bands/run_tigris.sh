#!/bin/zsh
# The expanding H II region (GOW17 + O3,S3,N3 and He++, 8 bands, recombination photons
# on the spot, 8 ranks) to t = 1.0 (0.978 Myr), one HDF5 dump per 0.02:
#   athinput.hii_dtype_ions          no dust, 128^3 (the default)
#   athinput.hii_dtype_ions_dusty64  dust, molecular initial gas, 64^3; build with
#                                    -shld_co -shld_ci so that CO and C I are shielded
# usage: run_tigris.sh WORKDIR [INPUT [CONFIGURE_FLAG...]]
#   The build goes to WORKDIR (athena, build.log), the run to WORKDIR/run (run.out,
#   run.err, time.out, run.log). slice_movie.py reads WORKDIR/run.
#   WORKDIR = $PHOTCHEM_RUNS/M5_hii_dtype_ions/ions_8bands_static_conditions for the
#   default input, .../ions_8bands_dusty64 for the dusty one.
source "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)/tests/env.sh"
W=${1:A}
INPUT=${2:-athinput.hii_dtype_ions}
(( $# >= 2 )) && shift 2 || shift $#
extra=("$@")
T=${TIGRIS_DIR:?set TIGRIS_DIR to a Tigris checkout}
HERE=${0:A:h}
tigris_build $W -mpi -gow17 --photchem_ions=O3,S3,N3 --gow17_bands=8 --prob=hii -hdf5 \
  "${extra[@]}" || { echo "build failed"; exit 1; }
sed "s|\${TIGRIS_DIR}|$T|g" $HERE/$INPUT > $W/$INPUT
run_logged $W/run $W mpi_launch 8 ../athena -i ../$INPUT
tail -3 $W/run/run.out
