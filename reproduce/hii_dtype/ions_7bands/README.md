Binary: tigris-gow17 78db305a3, configure as T6_multi_ion/bands7/build.sh (-mpi -gow17 --gow17_ions=O2,S3,N2 --gow17_bands=7 --prob=hii -hdf5).
Run (as tst/.../rayt_point/hii_dtype.py): mpirun -np 4 athena -i athinput.hii_dtype_ions time/tlim=1.0 photchem/cfl_photchem=10.0  (105 s, 436 cycles)
Analyze: python compare_ions.py <rundir> out.png
