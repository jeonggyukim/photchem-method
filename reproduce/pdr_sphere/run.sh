#!/bin/zsh
# A uniform sphere of 1000 H cm^-3, 2 pc in radius, in gas of 0.05 H cm^-3, lit from
# every side by one Draine field through the diffuse solver, post-processed to
# convergence with the NCR and the GOW17 photochemistry.  The input is
# athinput.cloud_postproc here, a copy of inputs/photchem/athinput.cloud_postproc
# (tigris branch gow17).  GOW17 is built with -shld_co, so CO is shielded along the
# rays.  Each run takes about a minute on eight cores.
#
#   ./run.sh <tigris checkout> [ncr|gow17]
# Runs go to ../../../data/pdr_sphere/<mode>; reduce.py then writes the mid-plane
# slices the figure script reads.
set -e
HERE=${0:A:h}
SRC=${1:?tigris checkout}
DATA=${HERE:h:h:h}/data/pdr_sphere
export OMPI_CXX=clang++
INC=$(mpicxx --showme:incdirs | tr ' ' '\n' | head -1)
common=(-mpi -shld_ci -shld_cr -hdf5 -h5double --prob=photchem_postproc
        --hdf5_path=/opt/homebrew/opt/hdf5-mpi --include=$INC
        --include=/opt/homebrew/opt/boost/include --cflag=-std=c++17)
for mode in ${2:-ncr gow17}; do
  if [[ $mode == gow17 ]]; then flags=(-gow17 -shld_co); args=(photchem/mode=gow17)
  else flags=(-ncr); args=(); fi
  (cd $SRC && python3 configure.py $common $flags > /dev/null && make clean > /dev/null \
     && make -j8 > /dev/null)
  run=$DATA/$mode; rm -rf $run; mkdir -p $run
  /bin/cp -f $SRC/bin/athena $SRC/inputs/tables/tigress_coolftn_ncr.txt $run/
  (cd $run && mpirun -np 8 ./athena -i $HERE/athinput.cloud_postproc $args > run.log 2>&1)
  grep -h "converged after\|did not converge" $run/run.log
done
