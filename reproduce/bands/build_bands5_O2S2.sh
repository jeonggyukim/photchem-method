cd "${TIGRIS_DIR:?set TIGRIS_DIR to a Tigris checkout}"
export OMPI_CXX=/opt/homebrew/bin/g++-16
python3 configure.py -mpi -gow17 --gow17_ions=O2,S2 --gow17_bands=5 --prob=hii -hdf5 --hdf5_path=/opt/homebrew/opt/hdf5-mpi --include=/opt/homebrew/opt/boost/include --fftw_path=/opt/homebrew/opt/fftw >/dev/null
make clean >/dev/null; make -j8 >$OLDPWD/build.out 2>&1; echo "make exit=$?"
cp bin/athena $OLDPWD/athena
