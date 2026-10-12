source "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)/tests/env.sh"
tigris_build "$PWD" -mpi -gow17 --gow17_ions=O2,S3 --gow17_bands=5 --prob=hii -hdf5
