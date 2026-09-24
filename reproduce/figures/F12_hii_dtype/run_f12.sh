#!/bin/zsh
# F12 run: hii_dtype with GOW17 + O3,S3,N3, seven bands, the arguments of
# tst/regression/scripts/tests/rayt_point/hii_dtype.py (gow17_ions mode) plus HDF5 dumps
# every 0.02 code. usage: run_f12.sh WORKDIR
W=${1:A}
T=~/Projects/tigris-gow17
tables=$T/inputs/tables
export OMPI_CXX=/opt/homebrew/bin/g++-16
(cd $T && python3 configure.py -mpi -gow17 --photchem_ions=O3,S3,N3 --gow17_bands=7 --prob=hii \
   -hdf5 --hdf5_path=/opt/homebrew/opt/hdf5-mpi --include=/opt/homebrew/opt/boost/include \
   --fftw_path=/opt/homebrew/opt/fftw >/dev/null && make clean >/dev/null \
   && make -j8 >$W/build.out 2>&1) || { echo "build failed"; exit 1; }
/bin/cp -f $T/bin/athena $W/athena
args=(output2/dt=0.02 output3/dt=-1 time/ncycle_out=10
      mesh/nx1=64 mesh/nx2=64 mesh/nx3=64 mesh/x1min=-18 mesh/x2min=-18 mesh/x3min=-18
      mesh/x1max=18 mesh/x2max=18 mesh/x3max=18 meshblock/nx1=32 meshblock/nx2=32
      meshblock/nx3=32 hydro/fofc=false hydro/neighbor_flooring=false problem/nH0=100.0
      problem/nH0_amb=100.0 problem/Qi=1e+49 problem/t0_src=0.0 problem/turb_flag=0
      rayt_point/rayt_point=true rayt_point/rays_per_cell=4 rayt_point/healpix_lev_min=4
      rayt_point/tau_max=30 rayt_point/rotate_rays=true photchem/flag_update_dt_main=true
      time/tlim=1.0 photchem/cfl_photchem=10.0
      photchem_gow17/rates_dir=$tables/rates
      photchem_gow17/ion_cooling_dir=$tables/chianti_v11/ion_cooling
      photchem_gow17/sed_file=$tables/sed/sb99_Z014_GenevaV00_2Myr.txt
      photchem_gow17/photoion_xsec_file=$tables/rates/verner96_photx.dat
      job/problem_id=hii_dtype_gow17_ions)
mkdir -p $W/run
(cd $W/run && /usr/bin/time -p mpirun -np 4 ../athena \
   -i $T/inputs/rayt_point/athinput.hii_dtype_gow17_ions "${args[@]}" >run.out 2>time.out)
echo "run exit=$?"
