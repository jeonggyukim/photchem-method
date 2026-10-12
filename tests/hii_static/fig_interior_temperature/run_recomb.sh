#!/bin/zsh
# Two runs of the post-processed Stromgren sphere (GOW17 + O3,S3,N3, 7 bands, 64^3, 8 ranks)
# with <photchem_gow17>/recombination_ots = true and = false (the diffuse ray tracer on,
# 11 diffuse groups: 7 point bands and 4 recombination groups).
# usage: run_recomb.sh WORKDIR   (athinput.stromgren_recomb is copied there)
#   The build goes to WORKDIR (athena, build.log) and each run to WORKDIR/ots_true and
#   WORKDIR/ots_false (run.out, run.err, time.out, run.log).
source "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)/tests/env.sh"
W=${1:A}
T=${TIGRIS_DIR:?set TIGRIS_DIR to a Tigris checkout}
tigris_build $W -mpi -gow17 --photchem_ions=O3,S3,N3 --gow17_bands=7 \
   --ngroup_rayt_diffuse=11 --prob=photchem_postproc -hdf5 || { echo "build failed"; exit 1; }
sed "s|\${TIGRIS_DIR}|$T|g" ${0:A:h}/athinput.stromgren_recomb > $W/athinput.stromgren_recomb
run_logged $W/ots_true $W mpi_launch 8 ../athena -i ../athinput.stromgren_recomb
run_logged $W/ots_false $W mpi_launch 8 ../athena -i ../athinput.stromgren_recomb \
   rayt_diffuse/rayt_diffuse=true photchem_gow17/recombination_ots=false
grep -h "converged\|did not" $W/ots_true/cloud.postproc.txt $W/ots_false/cloud.postproc.txt
