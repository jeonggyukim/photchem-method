#!/bin/bash
#SBATCH -J hiidiag2
#SBATCH -p normal
#SBATCH -N 2
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=64
#SBATCH --mem-per-cpu=2G
#SBATCH -t 3:00:00
#SBATCH -o %x_%j.out
# Timestep rules at NX^3 (default 128), n_H = 1e3, dust, molecular start, 128 ranks,
# each case to t = 0.05:
#   all     dt_user from the cooling, H+, H2 and CO timescales (cfl_photchem = 10)
#   hponly  dt_user from the H+ timescale alone (photchem/dt_limit = Hp)
#   nolimit the chemistry limits only the first step (flag_update_dt_main = false)
# usage: sbatch diag2.sbatch [NX]   (run in the directory holding athinput.n1e3_256)
set -e
TREE=/gpfs/jeonggyukim/tigris-photchem/tree_7838f2cbd_wip2
RUN=$PWD
BUILD=$RUN/build2
set +u
source ~/ai-notes/scripts/grammar_modules.sh
set -u
if [ ! -x $BUILD/athena ]; then
  mkdir -p $BUILD
  cd $TREE
  python3 configure.py --prob=hii -gow17 --photchem_ions=O3,S3,N3 --gow17_bands=8 \
    -shld_co -shld_ci -hdf5 --hdf5_path=$HOME/libs/hdf5_gnu "${ATHENA_GRAMMAR_FLAGS[@]}" \
    > $BUILD/configure.out
  make clean > /dev/null
  make -j 64 > $BUILD/build.out 2>&1
  cp -f bin/athena $BUILD/athena
  cd $RUN
fi
export LD_LIBRARY_PATH=$HOME/libs/hdf5_gnu/lib:${LD_LIBRARY_PATH:-}
# dt_limit can be set on the command line only once the input holds it
awk '{print} /^<photchem>$/{print "dt_limit = cool,Hp,H2,CO  # timescales that limit dt_user"}' \
  $RUN/athinput.n1e3_256 > $RUN/athinput.n1e3_256_dtlimit
NX=${1:-128}
if [ $NX = 128 ]; then mb=(meshblock/nx1=32 meshblock/nx2=32 meshblock/nx3=16); else mb=(); fi
for case in all hponly nolimit; do
  mkdir -p $RUN/v2_n$NX/$case
  cd $RUN/v2_n$NX/$case
  case $case in
    all) opts=(photchem/cfl_photchem=10) ;;
    hponly) opts=(photchem/cfl_photchem=10 photchem/dt_limit=Hp) ;;
    nolimit) opts=(photchem/flag_update_dt_main=false) ;;
  esac
  start=$(date +%s)
  mpirun -np $SLURM_NTASKS --mca btl ^ofi --bind-to core --map-by core \
    -x LD_LIBRARY_PATH $BUILD/athena -i $RUN/athinput.n1e3_256_dtlimit time/tlim=0.05 \
    mesh/nx1=$NX mesh/nx2=$NX mesh/nx3=$NX "${mb[@]}" \
    output1/dt=1e-5 output2/dt=1e-2 "${opts[@]}" > run.out 2> run.err
  echo "real $(( $(date +%s) - start ))" > time.out
  cd $RUN
done
