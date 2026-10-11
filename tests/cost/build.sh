# Builds the NCR and GOW17 (core, 3 bands) pdr_slab binaries into build_ncr/athena and
# build_gow17/athena of this directory.
source "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)/tests/env.sh"
set -e
D=$(cd "$(dirname "$0")" && pwd)
common=(-mpi -shld_ci -shld_cr -hdf5 -h5double --prob=photchem_postproc)
tigris_build "$D/build_ncr" -ncr "${common[@]}"
tigris_build "$D/build_gow17" -gow17 -shld_co "${common[@]}"
