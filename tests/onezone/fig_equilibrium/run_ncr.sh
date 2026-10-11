#!/bin/zsh
# F06 NCR reference: tst/regression/scripts/tests/photchem/equil.py, iprob = 0 (density
# sweep), z_gas = z_dust = 1, interp_mode = 1. usage: run_ncr.sh WORKDIR
W=${1:A}
T=${TIGRIS_DIR:?set TIGRIS_DIR to a Tigris checkout}
(cd $T && python3 configure.py -ncr --prob=photchem_equil \
   --include=/opt/homebrew/opt/boost/include --fftw_path=/opt/homebrew/opt/fftw >/dev/null \
   && make clean >/dev/null && make -j8 >$W/build.out 2>&1) || { echo "build failed"; exit 1; }
/bin/cp -f $T/bin/athena $W/athena
/bin/cp -f $T/inputs/tables/tigress_coolftn_ncr.txt $W/
args=(job/problem_id=equil_ncr output1/file_type=tab output1/variable=prim,uov time/ops_task=true
      hydro/active=fixed problem/iprob=0 photchem/mode=ncr photchem/bookkeeping=true
      photchem_ncr/cool_dust_flag=0 photchem_ncr/h2_diss_bg_flag=1
      photchem_ncr/hi_phot_bg_flag=0 photchem_ncr/zeta_hi_phot0=0.0 photchem_ncr/nsub_max=50
      photchem_ncr/cfl_cool_sub=0.1 photchem_ncr/z_gas=1.0 photchem_ncr/z_dust=1.0
      photchem_ncr/chi0=1.0 photchem_ncr/xi_cr0=2e-16 photchem_ncr/shld_n_hyd0=1e2
      photchem_ncr/shld_len0=0.0 photchem_ncr/shld_pow_idx=-0.7
      photchem_ncr/shld_column_hyd_cr0=9.35e20
      'photchem_ncr/hnu[0]=18.0' 'photchem_ncr/sigma_d[0]=1.0e-21'
      'photchem_ncr/sigma_pi_HI[0]=3.1e-18' 'photchem_ncr/sigma_pi_H2[0]=4.6e-18'
      'photchem_ncr/hnu[1]=12.2' 'photchem_ncr/sigma_d[1]=2.0e-21'
      'photchem_ncr/sigma_pi_HI[1]=0.0' 'photchem_ncr/sigma_pi_H2[1]=0.0'
      'photchem_ncr/hnu[2]=9.0' 'photchem_ncr/sigma_d[2]=1e-21'
      'photchem_ncr/sigma_pi_HI[2]=0.0' 'photchem_ncr/sigma_pi_H2[2]=0.0'
      photchem_ncr/coolftn_file=tigress_coolftn_ncr.txt photchem_ncr/interp_mode=1)
(cd $W && ./athena -i $T/inputs/photchem/athinput.photchem_equil "${args[@]}" >run.out 2>&1)
echo "run exit=$?"
