#!/bin/zsh
# The static H II region of ../fig_ion_stages (GOW17 + O3,S3,N3 and He++, 8 bands, 64^3,
# 8 ranks, 20 point rays per cell, recombination photons on the spot and carried by the
# diffuse ray tracer) with the heating and cooling by process written out:
#   <problem>/rates_by_process = true   every cell, as user output variables heat_*, cool_*
#   two probes                    cloud.probe0.txt and cloud.probe1.txt, the cells
#                                 holding r = 1.0 and 2.5 pc on the (1,1,1) diagonal
#                                 (cell centres at 1.14 and 2.44 pc): the terms over
#                                 T = 1e3-3e4 K with the cell's abundances and radiation
#                                 held fixed
# R_s = 2.93 pc for bb5e4, so the probe cells sit at 0.39 and 0.83 R_s.
# usage: run_tigris.sh WORKDIR [SPECTRUM...]   (default: bb5e4)
#   plot.py reads WORKDIR = $PHOTCHEM_RUNS/T6_multi_ion/hii_static_heat_cool.
source "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)/tests/env.sh"
W=${1:A}
shift
T=${TIGRIS_DIR:?set TIGRIS_DIR to a Tigris checkout}
HERE=${0:A:h}
tigris_build $W -mpi -gow17 --photchem_ions=O3,S3,N3 --gow17_bands=8 \
   --ngroup_rayt_diffuse=12 --prob=photchem_postproc -hdf5 || { echo "build failed"; exit 1; }
# Probe positions [pc]. The source is at the origin; a point at radius r on the (1,1,1)
# diagonal has x1 = x2 = x3 = r/sqrt(3).
d1=0.57735027  # probe 0: r = 1.0 pc
d2=1.44337567  # probe 1: r = 2.5 pc
sed -e "s|\${TIGRIS_DIR}|$T|g" \
    -e "/^<problem>/a\\
rates_by_process = true\\
n_probe = 2\\
probe0_x1 = $d1\\
probe0_x2 = $d1\\
probe0_x3 = $d1\\
probe1_x1 = $d2\\
probe1_x2 = $d2\\
probe1_x3 = $d2\\
probe_tmin = 1.0e3\\
probe_tmax = 3.0e4\\
probe_nt = 141
" $HERE/../fig_ion_stages/athinput.stromgren_stages > $W/athinput.heat_cool
typeset -A SED
SED=(sb99 $T/inputs/tables/sed/sb99_Z014_GenevaV00_2Myr.txt
     bb5e4 $HERE/../blackbody/T5e4/bb_T5e4.txt)
spectra=("$@")
(( ${#spectra} )) || spectra=(bb5e4)
for s in "${spectra[@]}"; do
  run_logged $W/${s}_ots $W mpi_launch 8 ../athena -i ../athinput.heat_cool \
     photchem_gow17/sed_file=${SED[$s]}
  run_logged $W/${s}_no_ots $W mpi_launch 8 ../athena -i ../athinput.heat_cool \
     photchem_gow17/sed_file=${SED[$s]} rayt_diffuse/rayt_diffuse=true \
     photchem_gow17/recombination_ots=false
done
grep -H "converged\|did not" $W/*_ots/cloud.postproc.txt $W/*_no_ots/cloud.postproc.txt
