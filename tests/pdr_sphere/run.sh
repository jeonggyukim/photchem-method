#!/bin/zsh
# A uniform sphere of 1000 H cm^-3, 2 pc in radius, in gas of 0.05 H cm^-3, lit from
# every side by one Draine field through the diffuse ray tracer, post-processed to
# convergence with the NCR and the GOW17 photochemistry.  The input is
# athinput.cloud_postproc here, a copy of inputs/photchem/athinput.cloud_postproc
# (tigris branch gow17).  GOW17 is built with -shld_co, so CO is shielded along the
# rays.  Each run takes about a minute on eight cores.
#
#   ./run.sh <tigris checkout> [ncr|gow17]
# Runs go to $PHOTCHEM_DATA_DIR/pdr_sphere/<mode>; reduce.py then writes the slices
# through the cloud centre that the figure script reads.
source "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)/tests/env.sh"
set -e
HERE=${0:A:h}
SRC=${1:?tigris checkout}
export TIGRIS_DIR=$SRC
DATA=${PHOTCHEM_DATA_DIR:?the directory holding the large reduced arrays}/pdr_sphere
INC=$(mpicxx --showme:incdirs | tr ' ' '\n' | head -1)
common=(-mpi -shld_ci -shld_cr -hdf5 -h5double --prob=photchem_postproc
        --include=$INC --cflag=-std=c++17)
for mode in ${2:-ncr gow17}; do
  if [[ $mode == gow17 ]]; then flags=(-gow17 -shld_co); args=(photchem/mode=gow17)
  else flags=(-ncr); args=(); fi
  run=$DATA/$mode; rm -rf $run; mkdir -p $run
  tigris_build $run $common $flags
  /bin/cp -f $SRC/inputs/tables/tigress_coolftn_ncr.txt $run/
  run_logged $run $run mpi_launch 8 ./athena -i $HERE/athinput.cloud_postproc $args
  grep -h "converged after\|did not converge" $run/run.out
done
