#!/bin/zsh
# F14 timing: hii_dtype to t = 1 code with Simple, NCR, GOW17 core and GOW17 + O3,S3,N3,
# the regression test's arguments (tst/.../rayt_point/hii_dtype.py), no HDF5, 4 ranks,
# one run at a time. GOW17 core uses the same seven bands as GOW17 + ions.
# usage: run_hii_cost.sh WORKDIR
source "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)/tests/env.sh"
W=${1:A}
T=${TIGRIS_DIR:?set TIGRIS_DIR to a Tigris checkout}
tables=$T/inputs/tables
build() {  # name configure-args...
  local name=$1; shift
  tigris_build $W/build_$name -mpi --prob=hii "$@" || exit 1
}
common=(output2/dt=-1 output3/dt=-1 time/ncycle_out=10
        mesh/nx1=64 mesh/nx2=64 mesh/nx3=64 mesh/x1min=-18 mesh/x2min=-18 mesh/x3min=-18
        mesh/x1max=18 mesh/x2max=18 mesh/x3max=18 meshblock/nx1=32 meshblock/nx2=32
        meshblock/nx3=32 hydro/fofc=false hydro/neighbor_flooring=false problem/nH0=100.0
        problem/nH0_amb=100.0 problem/Qi=1e+49 problem/t0_src=0.0 problem/turb_flag=0
        rayt_point/rayt_point=true rayt_point/rays_per_cell=4 rayt_point/healpix_lev_min=4
        rayt_point/tau_max=30 rayt_point/rotate_rays=true photchem/flag_update_dt_main=true
        time/tlim=1.0)
gow=(photchem/cfl_photchem=10.0 photchem_gow17/rates_dir=$tables/rates
     photchem_gow17/ion_cooling_dir=$tables/chianti_v11/ion_cooling
     photchem_gow17/sed_file=$tables/sed/sb99_Z014_GenevaV00_2Myr.txt
     photchem_gow17/photoion_xsec_file=$tables/rates/verner96_photx.dat)
run() {  # name input args...
  local name=$1 input=$2; shift 2
  mkdir -p $W/run_$name
  /bin/cp -f $tables/tigress_coolftn_ncr.txt $W/run_$name/
  run_logged $W/run_$name $W/build_$name mpi_launch 4 ../build_$name/athena -i $input \
    "${common[@]}" "$@" job/problem_id=hii
  grep -q "FATAL" $W/run_$name/run.out && echo "$name FATAL" || echo "$name done"
}
build simple --nspecies=1 --ngroup_rayt_point=1
run simple $T/inputs/rayt_point/athinput.hii_dtype photchem/mode=simple \
  photchem_simple/f_dt_rad=0.5 'photchem_simple/sigma_d[0]=0.0' \
  'photchem_simple/sigma_pi_HI[0]=3.0e-18' photchem_simple/tgas_HII=8000 \
  photchem_simple/tgas_HI=100
build ncr -ncr --ngroup_rayt_point=3
run ncr $T/inputs/rayt_point/athinput.hii_dtype photchem/mode=ncr photchem/cfl_photchem=10.0
build core -gow17 --gow17_bands=7
run core $T/inputs/rayt_point/athinput.hii_dtype_gow17_ions "${gow[@]}"
build ions -gow17 --photchem_ions=O3,S3,N3 --gow17_bands=7
run ions $T/inputs/rayt_point/athinput.hii_dtype_gow17_ions "${gow[@]}"
