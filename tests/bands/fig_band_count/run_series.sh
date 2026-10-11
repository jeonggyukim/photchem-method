#!/bin/zsh
# F10 band-count series: the postproc Stromgren sphere (64^3, 8 ranks) with 3, 5, 6, 7
# point-source bands, GOW17 + O3,S3,N3, from the current tigris-gow17 tree.
# usage: run_series.sh WORKDIR   (WORKDIR holds athinput.stromgren_postproc)
#   Each build goes to WORKDIR/build_b<n> (athena, build.log) and each run to
#   WORKDIR/run64_b<n> (run.out, run.err, time.out, run.log).
source "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)/tests/env.sh"
W=${1:A}
for nb in 3 5 6 7; do
  tigris_build $W/build_b$nb -mpi -gow17 --photchem_ions=O3,S3,N3 --gow17_bands=$nb \
     --prob=photchem_postproc -hdf5 || { echo "build b$nb failed"; exit 1; }
  run_logged $W/run64_b$nb $W/build_b$nb mpi_launch 8 ../build_b$nb/athena \
     -i ../athinput.stromgren_postproc
  echo "b$nb exit=$?"
done
