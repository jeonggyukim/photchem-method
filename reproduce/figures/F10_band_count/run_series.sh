#!/bin/zsh
# F10 band-count series: the postproc Stromgren sphere (64^3, 8 ranks) with 3, 5, 6, 7
# point-source bands, GOW17 + O3,S3,N3, from the current tigris-gow17 tree.
# usage: run_series.sh WORKDIR   (WORKDIR holds athinput.stromgren_postproc)
W=${1:A}
T=~/Projects/tigris-gow17
export OMPI_CXX=/opt/homebrew/bin/g++-16
for nb in 3 5 6 7; do
  (cd $T && python3 configure.py -mpi -gow17 --photchem_ions=O3,S3,N3 --gow17_bands=$nb \
     --prob=photchem_postproc -hdf5 --hdf5_path=/opt/homebrew/opt/hdf5-mpi \
     --include=/opt/homebrew/opt/boost/include --fftw_path=/opt/homebrew/opt/fftw >/dev/null \
   && make clean >/dev/null && make -j8 >$W/build_b$nb.out 2>&1) || { echo "build b$nb failed"; exit 1; }
  /bin/cp -f $T/bin/athena $W/athena_b$nb
  mkdir -p $W/run64_b$nb
  (cd $W/run64_b$nb && /usr/bin/time -p mpirun -np 8 ../athena_b$nb -i ../athinput.stromgren_postproc >run.out 2>time.out)
  echo "b$nb exit=$?"
done
