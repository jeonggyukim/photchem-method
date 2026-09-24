cd ~/Projects/tigris-gow17
export OMPI_CXX=/opt/homebrew/bin/g++-16
python3 configure.py -mpi -gow17 --gow17_ions=O2,S3,N2 --gow17_bands=6 --prob=hii -hdf5 --hdf5_path=/opt/homebrew/opt/hdf5-mpi --include=/opt/homebrew/opt/boost/include --fftw_path=/opt/homebrew/opt/fftw >/dev/null
make clean >/dev/null; make -j8 >$OLDPWD/build.out 2>&1; echo "make exit=$?"
cp bin/athena $OLDPWD/athena
