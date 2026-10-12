#!/bin/zsh
# The Wolfire, Hollenbach & McKee (2010) n ~ 1/r cloud (athinput.cloud_w10), lit from
# every side by G0' = 10 through the diffuse ray tracer, post-processed to convergence
# with the GOW17 and/or the NCR photochemistry.
#
#   ./run_w10.sh <tigris checkout> [gow17|ncr]   (default gow17)
# Runs go to $PHOTCHEM_RUNS/T8_pdr_sphere_w10/<mode>.
source "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)/tests/env.sh"
set -e
HERE=${0:A:h}
SRC=${1:?tigris checkout}
export TIGRIS_DIR=$SRC
DATA=${PHOTCHEM_RUNS:?the directory of the raw simulation output}/T8_pdr_sphere_w10
INC=$(mpicxx --showme:incdirs | tr ' ' '\n' | head -1)
common=(-mpi -shld_ci -shld_cr -hdf5 -h5double --prob=photchem_postproc
        --include=$INC --cflag=-std=c++17)
for mode in ${2:-gow17}; do
  if [[ $mode == gow17 ]]; then flags=(-gow17 -shld_co); args=(photchem/mode=gow17)
  else flags=(-ncr); args=(); fi
  run=$DATA/$mode; rm -rf $run; mkdir -p $run
  tigris_build $run $common $flags
  /bin/cp -f $SRC/inputs/tables/tigress_coolftn_ncr.txt $run/
  (cd $run && $PYTHON $HERE/make_density_w10.py w10_128.f32 128)
  run_logged $run $run mpi_launch 8 ./athena -i $HERE/athinput.cloud_w10 $args
  grep -h "converged after\|did not converge" $run/run.out
done
