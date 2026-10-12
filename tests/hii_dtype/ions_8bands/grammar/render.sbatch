#!/bin/bash
#SBATCH -J hiirender
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH -t 2:00:00
#SBATCH -o %x_%j.out
# Slice frames of an hii_dtype run with slice_movie.py, on one node.
# usage: sbatch render.sbatch RUNDIR OUTDIR TITLE
#   (photchem-method copied to /gpfs/jeonggyukim/tigris-photchem/photchem-method)
set -e
PM=/gpfs/jeonggyukim/tigris-photchem/photchem-method
export TIGRIS_DIR=/gpfs/jeonggyukim/tigris-photchem/tree_7838f2cbd_wip2
export PHOTCHEM_RUNS=/gpfs/jeonggyukim/tigris-photchem MPLBACKEND=Agg
export PYTHONPATH=$PM
cd $PM/tests/hii_dtype/ions_8bands
~/miniconda3/envs/pyathena/bin/python slice_movie.py "$1" "$2" --nproc 16 \
  --molecular --radiation --temp-log --density-scale 10 --title "$3"
