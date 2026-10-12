#!/bin/bash
#SBATCH -J hii256hp
#SBATCH -p normal
#SBATCH -N 2
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=64
#SBATCH --mem-per-cpu=2G
#SBATCH -t 12:00:00
#SBATCH -o %x_%j.out
# Expanding H II region, 256^3, n_H = 1e3, dust, molecular start, 128 ranks, to t = 1,
# with dt_user from the H+ timescale alone (photchem/dt_limit = Hp; diag2.sbatch).
# usage: sbatch run_hp.sbatch [RESTART]   (in a run directory beside ../diag; RESTART,
# a restart file to continue from)
set -e
BUILD=/gpfs/jeonggyukim/tigris-photchem/hii_dtype_n1e3_256/diag/build2
INPUT=/gpfs/jeonggyukim/tigris-photchem/hii_dtype_n1e3_256/diag/athinput.n1e3_256_dtlimit
set +u
source ~/ai-notes/scripts/grammar_modules.sh
set -u
export LD_LIBRARY_PATH=$HOME/libs/hdf5_gnu/lib:${LD_LIBRARY_PATH:-}
if [ $# -gt 0 ]; then start_args=(-r "$1"); else start_args=(-i $INPUT); fi
start=$(date +%s)
mpirun -np $SLURM_NTASKS --mca btl ^ofi --bind-to core --map-by core \
  -x LD_LIBRARY_PATH $BUILD/athena "${start_args[@]}" photchem/dt_limit=Hp \
  >> run.out 2>> run.err
echo "real $(( $(date +%s) - start ))" >> time.out
