#!/bin/zsh
# run_tigris.sh with the NCR chemical network in place of GOW17 (8 ranks), for the NCR
# twin of a GOW17 input, such as athinput.hii_cloud_m5n20_64_ncr.
# usage: run_tigris_ncr.sh WORKDIR INPUT [CONFIGURE_FLAG...]
#   The build goes to WORKDIR (athena, build.log), the run to WORKDIR/run (run.out,
#   run.err, time.out, run.log).
source "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)/tests/env.sh"
W=${1:A}
INPUT=$2
shift 2
extra=("$@")
T=${TIGRIS_DIR:?set TIGRIS_DIR to a Tigris checkout}
HERE=${0:A:h}
tigris_build $W -mpi -ncr --prob=hii -hdf5 "${extra[@]}" || { echo "build failed"; exit 1; }
sed "s|\${TIGRIS_DIR}|$T|g" $HERE/$INPUT > $W/$INPUT
run_logged $W/run $W mpi_launch 8 ../athena -i ../$INPUT
tail -3 $W/run/run.out
