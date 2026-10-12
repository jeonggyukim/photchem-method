#!/bin/bash
#SBATCH -J hii256
#SBATCH -p normal
#SBATCH -N 2
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=64
#SBATCH --mem-per-cpu=2G
#SBATCH -t 12:00:00
#SBATCH -o %x_%j.out
# Expanding H II region, 256^3, n_H = 1e3, dust, molecular start, 128 ranks.
# usage: sbatch [-t TIME] run.sbatch [NLIM]   (NLIM > 0: a smoke run of NLIM cycles)
set -e
NLIM=${1:--1}
TREE=/gpfs/jeonggyukim/tigris-photchem/tree_f698fe8b9
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
# the parallel HDF5 is a shared library outside the module tree
export LD_LIBRARY_PATH=$HOME/libs/hdf5_gnu/lib:${LD_LIBRARY_PATH:-}
start=$(date +%s)
# --mca btl ^ofi: across nodes the OFI BTL segfaulted in MPI_Iprobe (libfabric 1.18,
# job 1176336); UCX carries the inter-node traffic instead
mpirun -np $SLURM_NTASKS --mca btl ^ofi --bind-to core --map-by core \
  -x LD_LIBRARY_PATH ./athena -i athinput.n1e3_256 time/nlim=$NLIM > run.out 2> run.err
echo "real $(( $(date +%s) - start ))" > time.out
