#!/bin/zsh
# The post-processed static H II region (GOW17 + O3,S3,N3 and He++, 8 bands, 64^3,
# 8 ranks, 20 point rays per cell) for each spectrum, with
# <photchem_gow17>/recombination_ots = true and = false (the diffuse ray tracer on,
# 12 diffuse groups: 8 point bands and 4 recombination groups):
#   sb99   Starburst99 at 2 Myr (the Tigris table)
#   bb5e4  a 5e4 K blackbody, about the hottest O star (../blackbody/T5e4)
#   bb1e5  a 1e5 K blackbody (../blackbody/T1e5); not run by default
# usage: run_tigris.sh WORKDIR [SPECTRUM...]   (default: sb99 bb5e4)
#   The build goes to WORKDIR (athena, build.log), each run to WORKDIR/<spectrum>_ots
#   and WORKDIR/<spectrum>_no_ots (run.out, run.err, time.out, run.log). plot.py reads
#   WORKDIR = $PHOTCHEM_RUNS/T6_multi_ion/hii_static_stages. The blackbody spectra are
#   written by ../blackbody/run_cloudy_bb.sh, which also runs the Cloudy models.
source "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)/tests/env.sh"
W=${1:A}
shift
T=${TIGRIS_DIR:?set TIGRIS_DIR to a Tigris checkout}
HERE=${0:A:h}
tigris_build $W -mpi -gow17 --photchem_ions=O3,S3,N3 --gow17_bands=8 \
   --ngroup_rayt_diffuse=12 --prob=photchem_postproc -hdf5 || { echo "build failed"; exit 1; }
sed "s|\${TIGRIS_DIR}|$T|g" $HERE/athinput.stromgren_stages > $W/athinput.stromgren_stages
typeset -A SED
SED=(sb99 $T/inputs/tables/sed/sb99_Z014_GenevaV00_2Myr.txt
     bb1e5 $HERE/../blackbody/T1e5/bb_T1e5.txt
     bb5e4 $HERE/../blackbody/T5e4/bb_T5e4.txt)
spectra=("$@")
(( ${#spectra} )) || spectra=(sb99 bb5e4)
for s in "${spectra[@]}"; do
  run_logged $W/${s}_ots $W mpi_launch 8 ../athena -i ../athinput.stromgren_stages \
     photchem_gow17/sed_file=${SED[$s]}
  run_logged $W/${s}_no_ots $W mpi_launch 8 ../athena -i ../athinput.stromgren_stages \
     photchem_gow17/sed_file=${SED[$s]} rayt_diffuse/rayt_diffuse=true \
     photchem_gow17/recombination_ots=false
done
grep -H "converged\|did not" $W/*_ots/cloud.postproc.txt $W/*_no_ots/cloud.postproc.txt
