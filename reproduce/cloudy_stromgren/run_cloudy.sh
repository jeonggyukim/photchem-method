#!/bin/zsh
# Cloudy 25.00 static Stromgren sphere. Run in this directory; takes ~1 min.
#   1. make_sed.py converts the SB99 table to Cloudy's `table SED` file and checks
#      the photon ratios Q(>24.59 eV)/Q and Q(>35.12 eV)/Q against the table.
#   2. cloudy.exe reads hii.in, writes hii.out and the save files.
#      Cloudy looks for the SED file on the data path, so this directory is appended.
#   3. make_profiles.py reduces the save files to radial_profiles.txt and prints
#      the H+-zone averages and coolant shares.
set -e
cd "${0:A:h}"
PY=/opt/homebrew/Caskroom/miniforge/base/envs/pyathena/bin/python
CLOUDY=/Users/jgkim/Documents/cloudy-grain-hii/cloudy/source/cloudy.exe
export CLOUDY_DATA_PATH=/Users/jgkim/Dropbox/Projects/cloudy/data:$PWD
$PY make_sed.py
$CLOUDY -r hii
$PY make_profiles.py
