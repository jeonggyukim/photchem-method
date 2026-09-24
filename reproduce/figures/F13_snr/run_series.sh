#!/bin/zsh
# F13 SNR series: classic cooling, NCR, GOW17 core, GOW17 + O3,S3,N3 with and without
# the returned grain C, O, Si in the hot-gas cooling. 64^3, 8 ranks, tlim 0.05, nghost 4
# as tst/regression/scripts/tests/feedback/rad_snr.py.
# usage: run_series.sh WORKDIR   (WORKDIR holds athinput.snr_gow17 and src_noreturn, a copy
# of the tigris-gow17 tree without the three returned-metal lines of SharedHighIonTable)
W=${1:A}
T=~/Projects/tigris-gow17
export OMPI_CXX=/opt/homebrew/bin/g++-16
common=(-mpi -fb --prob=radiative_snr --nghost=4 -hdf5 --hdf5_path=/opt/homebrew/opt/hdf5-mpi
        --include=/opt/homebrew/opt/boost/include --fftw_path=/opt/homebrew/opt/fftw)
reg=(time/ncycle_out=10 time/dt_diagnostics=-1 time/integrator=rk2 cooling/coolftn=tigress
     output2/dt=-1 output3/dt=-1 time/tlim=0.05 time/ops_task=true feedback/ischeme=Subcell)
build() {  # tree name extra-configure-args...
  local tree=$1 name=$2; shift 2
  (cd $tree && python3 configure.py "${common[@]}" "$@" >/dev/null && make clean >/dev/null \
     && make -j8 >$W/build_$name.out 2>&1) || { echo "build $name failed"; exit 1; }
  /bin/cp -f $tree/bin/athena $W/athena_$name
}
run() {  # name input args...
  local name=$1 input=$2; shift 2
  mkdir -p $W/run_$name
  /bin/cp -f $T/inputs/tables/tigress_coolftn.txt $T/inputs/tables/tigress_coolftn_ncr.txt $W/run_$name/
  (cd $W/run_$name && /usr/bin/time -p mpirun -np 8 ../athena_$name -i $input "$@" \
     job/problem_id=snr >run.out 2>time.out)
  echo "$name exit=$?"
}
build $T classic
run classic $T/inputs/feedback/athinput.radiative_snr "${reg[@]}" cooling/cooling=op_split \
  hydro/neighbor_flooring=true
build $T ncr -ncr
run ncr $T/inputs/feedback/athinput.radiative_snr "${reg[@]}" cooling/cooling=none \
  photchem/photchem=true photchem_ncr/cool_hyd_cie_flag=true
gow=(time/ncycle_out=10 time/dt_diagnostics=-1 output2/dt=0.005 output3/dt=-1 time/tlim=0.05
     cooling/cooling=none photchem/photchem=true photchem_gow17/cool_hyd_cie_flag=true)
build $T core -gow17
run core ../athinput.snr_gow17 "${gow[@]}"
build $T ions -gow17 --photchem_ions=O3,S3,N3
run ions ../athinput.snr_gow17 "${gow[@]}"
build $W/src_noreturn noreturn -gow17 --photchem_ions=O3,S3,N3
run noreturn ../athinput.snr_gow17 "${gow[@]}"
