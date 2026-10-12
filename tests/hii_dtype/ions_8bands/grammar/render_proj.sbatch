#!/bin/bash
#SBATCH -J hiiproj
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH -t 2:00:00
#SBATCH -o %x_%j.out
# Projection frames of an hii_dtype run with proj_movie.py, on one node.
# usage: sbatch render_proj.sbatch RUNDIR OUTDIR TITLE
set -e
PM=/gpfs/jeonggyukim/tigris-photchem/photchem-method
export TIGRIS_DIR=/gpfs/jeonggyukim/tigris-photchem/tree_7838f2cbd_wip2
export PHOTCHEM_RUNS=/gpfs/jeonggyukim/tigris-photchem MPLBACKEND=Agg
export PYTHONPATH=$PM
cd $PM/tests/hii_dtype/ions_8bands
~/miniconda3/envs/pyathena/bin/python proj_movie.py "$1" "$2" --nproc 12 --title "$3"
