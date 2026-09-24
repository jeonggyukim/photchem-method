set -e
D=$(cd "$(dirname "$0")" && pwd)
cd ~/Projects/tigris-gow17
export OMPI_CXX=/opt/homebrew/bin/g++-16
python3 configure.py -mpi -fb -gow17 --gow17_ions=O2,S3,N2 --prob=radiative_snr --include=/opt/homebrew/opt/boost/include --fftw_path=/opt/homebrew/opt/fftw >/dev/null
make clean >/dev/null; make -j8 >$D/build.out 2>&1; echo "make exit=$?"
cp bin/athena $D/athena
