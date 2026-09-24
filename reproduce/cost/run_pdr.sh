# pdr_slab, 4 ranks, exactly 5 postproc iterations (convergence tests off), NCR then GOW17.
D=$(cd "$(dirname "$0")" && pwd)
IN=~/Projects/tigris-gow17/inputs/photchem
for m in ncr gow17; do
  mkdir -p $D/pdr_$m && cd $D/pdr_$m && rm -f *
  cp ~/Projects/tigris-gow17/inputs/tables/tigress_coolftn_ncr.txt .
  inp=$IN/athinput.pdr_slab; [ $m = gow17 ] && inp=$IN/athinput.pdr_slab_gow17
  args=(postproc/niteration=5 postproc/threshold_abund_mean=-1 postproc/threshold_temp=-1 postproc/threshold_rad=-1)
  s=$(date +%s)
  mpirun -np 4 $D/athena_$m -i $inp "${args[@]}" >run.out 2>&1
  echo "$m exit=$? wall $(( $(date +%s)-s )) s"
done
