#!/bin/bash
#SBATCH -J hiidiag
#SBATCH -p normal
#SBATCH -N 2
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=64
#SBATCH --mem-per-cpu=2G
#SBATCH -t 3:00:00
#SBATCH -o %x_%j.out
# Which cells set dt_user at 256^3, n_H = 1e3, dust, molecular start, 128 ranks:
# 100 cycles each at cfl_photchem = 10 and 100, with the per-cell timescale outputs
# (t_cool, t_Hp, t_H2, t_CO) in every dump and a history row every cycle.
set -e
TREE=/gpfs/jeonggyukim/tigris-photchem/tree_7838f2cbd_wip
RUN=$PWD
set +u
source ~/ai-notes/scripts/grammar_modules.sh
set -u
if [ ! -x $RUN/athena ]; then
  cd $TREE
  python3 configure.py --prob=hii -gow17 --photchem_ions=O3,S3,N3 --gow17_bands=8 \
    -shld_co -shld_ci -hdf5 --hdf5_path=$HOME/libs/hdf5_gnu "${ATHENA_GRAMMAR_FLAGS[@]}" \
    > $RUN/configure.out
  make clean > /dev/null
  make -j 64 > $RUN/build.out 2>&1
  cp -f bin/athena $RUN/athena
  cd $RUN
fi
export LD_LIBRARY_PATH=$HOME/libs/hdf5_gnu/lib:${LD_LIBRARY_PATH:-}
# At NX^3 (128 MeshBlocks either way), each case to t = 0.02: cfl10, cfl100, dt_user
# limits the step; nolimit, the hydrodynamic step alone, the chemistry substepping
# within it (photchem/flag_update_dt_main = false).
NX=${1:-128}
if [ $NX = 128 ]; then mb=(meshblock/nx1=32 meshblock/nx2=32 meshblock/nx3=16); else mb=(); fi
for case in cfl10 cfl100 nolimit; do
  mkdir -p $RUN/n$NX/$case
  cd $RUN/n$NX/$case
  case $case in
    cfl10) opts=(photchem/cfl_photchem=10) ;;
    cfl100) opts=(photchem/cfl_photchem=100) ;;
    nolimit) opts=(photchem/flag_update_dt_main=false) ;;
  esac
  start=$(date +%s)
  mpirun -np $SLURM_NTASKS --mca btl ^ofi --bind-to core --map-by core \
    -x LD_LIBRARY_PATH $RUN/athena -i $RUN/athinput.n1e3_256 time/tlim=0.02 \
    mesh/nx1=$NX mesh/nx2=$NX mesh/nx3=$NX "${mb[@]}" \
    output1/dt=1e-5 output2/dt=4e-3 "${opts[@]}" > run.out 2> run.err
  echo "real $(( $(date +%s) - start ))" > time.out
  cd $RUN
done
