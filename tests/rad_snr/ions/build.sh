source "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)/tests/env.sh"
set -e
D=$(cd "$(dirname "$0")" && pwd)
tigris_build "$D" -mpi -fb -gow17 --gow17_ions=O2,S3,N2 --prob=radiative_snr
