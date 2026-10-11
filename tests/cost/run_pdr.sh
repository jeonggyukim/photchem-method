# pdr_slab, 4 ranks, exactly 5 postproc iterations (convergence tests off), NCR then GOW17,
# with the binaries of build.sh.
source "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)/tests/env.sh"
D=$(cd "$(dirname "$0")" && pwd)
T=${TIGRIS_DIR:?set TIGRIS_DIR to a Tigris checkout}
IN=$T/inputs/photchem
for m in ncr gow17; do
  mkdir -p $D/pdr_$m && cd $D/pdr_$m && rm -f *
  cp $T/inputs/tables/tigress_coolftn_ncr.txt .
  inp=$IN/athinput.pdr_slab; [ $m = gow17 ] && inp=$IN/athinput.pdr_slab_gow17
  args=(postproc/niteration=5 postproc/threshold_abund_mean=-1 postproc/threshold_temp=-1 postproc/threshold_rad=-1)
  run_logged $D/pdr_$m $D/build_$m mpi_launch 4 $D/build_$m/athena -i $inp "${args[@]}"
  echo "$m exit=$? wall $(cut -d' ' -f2 $D/pdr_$m/time.out) s"
done
