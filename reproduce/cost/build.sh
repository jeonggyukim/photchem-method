# Builds the NCR and GOW17 (core, 3 bands) pdr_slab binaries into this directory.
set -e
D=$(cd "$(dirname "$0")" && pwd)
cd ~/Projects/tigris-gow17
export OMPI_CXX=/opt/homebrew/bin/g++-16
common=(-mpi -shld_ci -shld_cr -hdf5 -h5double --prob=photchem_postproc --hdf5_path=/opt/homebrew/opt/hdf5-mpi --include=/opt/homebrew/opt/boost/include --fftw_path=/opt/homebrew/opt/fftw)
python3 configure.py -ncr "${common[@]}" >/dev/null
make clean >/dev/null; make -j8 >$D/build_ncr.out 2>&1; echo "ncr make exit=$?"
cp bin/athena $D/athena_ncr
python3 configure.py -gow17 -shld_co "${common[@]}" >/dev/null
make clean >/dev/null; make -j8 >$D/build_gow17.out 2>&1; echo "gow17 make exit=$?"
cp bin/athena $D/athena_gow17
