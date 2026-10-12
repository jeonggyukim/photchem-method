#!/bin/zsh
# The expanding H II region of athinput.hii_dtype_ions (GOW17 + O3,S3,N3 and He++,
# 8 bands, no dust, recombination photons on the spot, 128^3, 8 ranks) to t = 1.0
# (0.978 Myr), one HDF5 dump per 0.02.
# usage: run_tigris.sh WORKDIR
#   The build goes to WORKDIR (athena, build.log), the run to WORKDIR/run (run.out,
#   run.err, time.out, run.log). slice_movie.py reads WORKDIR/run.
#   WORKDIR = $PHOTCHEM_RUNS/M5_hii_dtype_ions/ions_8bands_static_conditions.
source "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)/tests/env.sh"
W=${1:A}
T=${TIGRIS_DIR:?set TIGRIS_DIR to a Tigris checkout}
HERE=${0:A:h}
tigris_build $W -mpi -gow17 --photchem_ions=O3,S3,N3 --gow17_bands=8 --prob=hii -hdf5 \
  || { echo "build failed"; exit 1; }
sed "s|\${TIGRIS_DIR}|$T|g" $HERE/athinput.hii_dtype_ions > $W/athinput.hii_dtype_ions
run_logged $W/run $W mpi_launch 8 ../athena -i ../athinput.hii_dtype_ions
tail -3 $W/run/run.out
