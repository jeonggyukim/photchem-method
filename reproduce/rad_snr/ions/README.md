Build: bash build.sh (tigris-gow17 78db305a3+, -mpi -fb -gow17 --gow17_ions=O2,S3,N2 --prob=radiative_snr).
Run (test NCR arguments, CIE hydrogen on): mpirun -np 4 athena -i athinput.radiative_snr_ions output2/dt=-1 output3/dt=-1 time/ncycle_out=10 time/dt_diagnostics=-1 time/integrator=rk2 cooling/coolftn=tigress time/tlim=0.05 cooling/cooling=none photchem/photchem=true photchem_gow17/cool_hyd_cie_flag=true time/ops_task=true feedback/ischeme=Subcell job/problem_id=snr_ions  (43 s; needs tigress_coolftn_ncr.txt in the run dir)
Measure: python ../snr_sf.py snr_ions.hst athinput.runtime
