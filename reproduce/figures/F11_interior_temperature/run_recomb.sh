#!/bin/zsh
# F11 runs: the post-processed Stromgren sphere (GOW17 + O3,S3,N3, 7 bands, 64^3, 8 ranks)
# with <photchem_gow17>/recombination_ots = true and = false (the diffuse tracer on,
# 11 diffuse groups: 7 point bands and 4 recombination groups).
# usage: run_recomb.sh WORKDIR   (athinput.stromgren_recomb is copied there)
W=${1:A}
T=~/Projects/tigris-gow17
export OMPI_CXX=/opt/homebrew/bin/g++-16
(cd $T && python3 configure.py -mpi -gow17 --photchem_ions=O3,S3,N3 --gow17_bands=7 \
   --ngroup_rayt_diffuse=11 --prob=photchem_postproc -hdf5 \
   --hdf5_path=/opt/homebrew/opt/hdf5-mpi --include=/opt/homebrew/opt/boost/include \
   --fftw_path=/opt/homebrew/opt/fftw >/dev/null && make clean >/dev/null \
   && make -j8 >$W/build.out 2>&1) || { echo "build failed"; exit 1; }
/bin/cp -f $T/bin/athena $W/athena
/bin/cp -f ${0:A:h}/athinput.stromgren_recomb $W/
mkdir -p $W/ots_true $W/ots_false
(cd $W/ots_true && /usr/bin/time -p mpirun -np 8 ../athena -i ../athinput.stromgren_recomb \
   >run.out 2>time.out)
(cd $W/ots_false && /usr/bin/time -p mpirun -np 8 ../athena -i ../athinput.stromgren_recomb \
   rayt_diffuse/rayt_diffuse=true photchem_gow17/recombination_ots=false >run.out 2>time.out)
grep -h "converged\|did not" $W/ots_true/cloud.postproc.txt $W/ots_false/cloud.postproc.txt
