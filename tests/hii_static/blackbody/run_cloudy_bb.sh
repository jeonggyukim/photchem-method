#!/bin/zsh
# Cloudy 25.00 static Stromgren spheres for ../fig_ion_stages, with the
# default transfer of the diffuse field (no_ots/) and with `diffuse OTS` (ots/), each
# about 1 min. The models go to $PHOTCHEM_RUNS/T6_multi_ion/cloudy_<spectrum>/; the
# figure reads their save files with `plot.py SPECTRUM --from-runs`.
#   sb99   the Starburst99 2 Myr models of ../hii.in and ../thermal_with_N_ots/hii.in
#   bb<T>  a blackbody of temperature T (e.g. bb5e4): make_bb_sed.py writes the same
#          spectrum for both codes into T<T>/, and the inputs are ../hii.in with it
# usage: run_cloudy_bb.sh SPECTRUM
set -e
source "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)/tests/env.sh"
HERE=${0:A:h}
CLOUDY=${CLOUDY_EXE:?set CLOUDY_EXE to the cloudy.exe of Cloudy 25.00}
DATA=${CLOUDY_DATA_DIR:?set CLOUDY_DATA_DIR to the Cloudy data directory}
S=${1:?give the spectrum: sb99 or bb<T>, e.g. bb5e4}
W=${PHOTCHEM_RUNS:?set PHOTCHEM_RUNS}/T6_multi_ion/cloudy_$S
mkdir -p $W/no_ots $W/ots
if [[ $S == sb99 ]]; then
  SEDDIR=$HERE/..
  /bin/cp -f $HERE/../hii.in $W/no_ots/hii.in
  /bin/cp -f $HERE/../thermal_with_N_ots/hii.in $W/ots/hii.in
else
  T=${S#bb}
  SEDDIR=$HERE/T$T
  mkdir -p $SEDDIR
  "$PYTHON" $HERE/make_bb_sed.py $T T$T
  SED="bb_T$T.sed"
  sed -e "s/^title .*/title static HII region, blackbody $T K, Q(H)=1e49, nH=100/" \
      -e "s/^table SED .*/table SED \"$SED\"/" $HERE/../hii.in > $W/no_ots/hii.in
  sed -e "s/^title .*/title static HII region, blackbody $T K, Q(H)=1e49, nH=100, OTS/" \
      -e "s/^table SED .*/table SED \"$SED\"/" $HERE/../thermal_with_N_ots/hii.in \
      > $W/ots/hii.in
fi
export CLOUDY_DATA_PATH=$DATA:$SEDDIR
for m in no_ots ots; do
  (cd $W/$m && $CLOUDY -r hii > run.log 2>&1 && "$PYTHON" $HERE/../make_profiles.py >> run.log)
done
