#!/bin/bash
# F15: a Bonnor-Ebert sphere (1000 H cm^-3 at the centre, 2 pc, xi_max 8.4) in 0.05 H cm^-3
# gas, lit by an isotropic 1-Draine FUV field through the diffuse solver (48 HEALPix
# directions), post-processed to convergence with NCR and with GOW17 + ions (O3,S3,N3).
# Same input (inputs/photchem/athinput.cloud_postproc), same bands (LyC off, LW, PE);
# only <photchem>/mode and the build differ.
#   bash run.sh [build|run|all]   (default all)
set -e
TIGRIS=${TIGRIS_DIR:?set TIGRIS_DIR to a Tigris checkout}
RUNS=${RUNS:-${PHOTCHEM_POSTPROC_RUNS:?set PHOTCHEM_POSTPROC_RUNS to the post-processing runs directory}/be_net}
INP=$TIGRIS/inputs/photchem/athinput.cloud_postproc
export OMPI_CXX=/opt/homebrew/bin/g++-16
MPIINC=$(mpicxx --showme:incdirs | tr ' ' '\n' | head -1)
common=(--prob=photchem_postproc -mpi -hdf5 --hdf5_path=/opt/homebrew/opt/hdf5-mpi
        --include=$MPIINC -shld_ci -shld_cr)
mkdir -p "$RUNS"

build() {
  local name=$1; shift
  (cd "$TIGRIS" && python3 configure.py "${common[@]}" "$@" > "$RUNS/configure_$name.log" \
     && make clean > /dev/null && make -j8 > "$RUNS/build_$name.log" 2>&1)
  mkdir -p "$RUNS/$name"
  /bin/cp -f "$TIGRIS/bin/athena" "$RUNS/$name/athena"
  (cd "$TIGRIS" && git rev-parse --short HEAD && git status --short src) > "$RUNS/$name/tree.txt"
}

run() {
  local name=$1 mode=$2; shift 2
  cd "$RUNS/$name"
  /bin/cp -f "$TIGRIS/inputs/tables/tigress_coolftn_ncr.txt" .
  /bin/cp -f "$INP" athinput.cloud_postproc
  if [ "$mode" = gow17 ]; then
    # the ion ladders read their rate tables from rates_dir (default: the run
    # directory); the command line cannot add a key, so write it into the copy
    awk -v d="$TIGRIS/inputs/tables/rates/" -v c="$TIGRIS/inputs/tables/chianti_v11/ion_cooling/" '{print} /^<photchem_gow17>/ {
      print "rates_dir = " d; print "photoion_xsec_file = " d "verner96_photx.dat";
      print "ion_cooling_dir = " c}' \
      "$INP" > athinput.cloud_postproc
  fi
  { time mpirun -np 8 ./athena -i athinput.cloud_postproc photchem/mode=$mode "$@" \
      > run.out 2> run.err ; } 2> time.txt
}

what=${1:-all}
if [ "$what" = gow17 ]; then
  build gow17 -gow17 --photchem_ions=O3,S3,N3 -shld_co
  run gow17 gow17
  exit 0
fi
if [ "$what" = run_gow17 ]; then
  run gow17 gow17
  exit 0
fi
if [ "$what" = build ] || [ "$what" = all ]; then
  build ncr -ncr
  build gow17 -gow17 --photchem_ions=O3,S3,N3 -shld_co
fi
if [ "$what" = run ] || [ "$what" = all ]; then
  run ncr ncr
  run gow17 gow17
fi
