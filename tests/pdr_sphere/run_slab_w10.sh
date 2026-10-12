#!/bin/zsh
# The Wolfire et al. (2010) cloud as their plane-parallel slab (athinput.slab_w10), at
# the mean density (f = 1) and at their clump density (f = 4), with the GOW17 binary
# run_w10.sh built for the sphere.
#
#   ./run_slab_w10.sh [F [NX]]   (after ./run_w10.sh <tigris checkout> gow17)
#     F: 1, 4 or both (default both); NX: cells along the beam (default 4096; 512
#     for f = 1 and 1024 for f = 4 already resolve both transitions, 0.009 mag a cell
#     at the surface)
# Runs go to $PHOTCHEM_RUNS/T8_pdr_sphere_w10/slab_f<F>_n<NX>, on 8 ranks of NX/8 cells.
source "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)/tests/env.sh"
set -e
HERE=${0:A:h}
BASE=${PHOTCHEM_RUNS:?the directory of the raw simulation output}/T8_pdr_sphere_w10
BUILD=$BASE/gow17
NX=${2:-4096}
for f in ${1:-1 4}; do
  run=$BASE/slab_f${f}_n$NX; rm -rf $run; mkdir -p $run
  len=$(( 40.96/f ))
  dx=$(( len/NX ))
  /bin/cp -f $BUILD/tigress_coolftn_ncr.txt $run/
  (cd $run && $PYTHON $HERE/make_slab_w10.py slab_f1.f32 $NX $len --factor $f)
  # cubic cells: rayt_diffuse spaces its rays by the x1 cell width, so a cross-section
  # narrower than two such widths holds fewer rays than cell rows and leaves rows unlit
  run_logged $run $BUILD mpi_launch 8 $BUILD/athena -i $HERE/athinput.slab_w10 \
    photchem/mode=gow17 mesh/x1max=$len mesh/nx1=$NX meshblock/nx1=$(( NX/8 )) \
    mesh/x2min=$(( -dx )) mesh/x2max=$dx mesh/x3min=$(( -dx )) mesh/x3max=$dx
  grep -h "converged after\|did not converge" $run/run.out
done
