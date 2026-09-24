Binary: tigris-gow17 78db305a3, configure as T6_multi_ion/bands7/build.sh (-mpi -gow17 --gow17_ions=O2,S3,N2 --gow17_bands=7 --prob=hii -hdf5).
Run: mpirun -np 4 athena -i athinput.hii_dtype_ions time/tlim=1.0  (422 s, 1900 cycles)
Analyze: python compare_ions.py <rundir> out.png
