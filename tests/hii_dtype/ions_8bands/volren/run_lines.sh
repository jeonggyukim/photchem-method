#!/bin/zsh
# H-alpha, [O III] 5007 and [S II] 6716, 6731 images of hii_dtype snapshots (GOW17 +
# O3,S3,N3) from the Tigris volume renderer: volren_postproc reads the run's athdf dumps
# with the chemical network constructed and not advanced, and a parallel camera looks
# along -z over the whole box, one pixel per cell. Once without dust and once with the
# Weingartner & Draine Milky Way grains at Z_d = 1.
#
#   ./run_lines.sh RUNDIR NUM_START [NUM_END [NUM_STEP]]
#     RUNDIR: the run directory (athinput.runtime and HII.out2.NNNNN.athdf)
# The physics blocks come from RUNDIR/athinput.runtime, so the renderer sees the
# network settings the run had. Output goes to RUNDIR/../volren_lines/{thin,dust}.
source "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)/tests/env.sh"
set -e
RUN=${1:A}
N0=$2
N1=${3:-$2}
DN=${4:-1}
T=${TIGRIS_DIR:?set TIGRIS_DIR to a Tigris checkout}
OUT=${RUN:h}/volren_lines
BUILD=$OUT/build
tigris_build $BUILD -mpi -gow17 --photchem_ions=O3,S3,N3 --gow17_bands=8 \
  --prob=volren_postproc -hdf5 -shld_co -shld_ci -rayt_volren --nbands_rayt_volren=4
nx=$(awk '/^<mesh>/{m=1} m&&/^nx1/{print $3; exit}' $RUN/athinput.runtime)
lo=$(awk '/^<mesh>/{m=1} m&&/^x1min/{print $3; exit}' $RUN/athinput.runtime)
hi=$(awk '/^<mesh>/{m=1} m&&/^x1max/{print $3; exit}' $RUN/athinput.runtime)
width=$(( hi - lo ))
for mode in thin dust; do
  dir=$OUT/$mode
  mkdir -p $dir
  zd=$([[ $mode == dust ]] && echo 1.0 || echo 0.0)
  {
    cat <<EOF
<comment>
problem = volume-rendered line images of $RUN

<job>
problem_id = lines

<output1>
file_type        = volren
dcycle           = 1
file_id          = lines
ncameras         = 1
format           = h5
frames_per_output = 1
interp           = false
cam0_npx              = $nx
cam0_npy              = $nx
cam0_projection       = parallel
cam0_plane_center_x   = $(( 0.5*(lo + hi) ))
cam0_plane_center_y   = $(( 0.5*(lo + hi) ))
cam0_plane_center_z   = $(( hi + 1.0 ))
cam0_look_direction_x = 0.0
cam0_look_direction_y = 0.0
cam0_look_direction_z = -1.0
cam0_up_x             = 0.0
cam0_up_y             = 1.0
cam0_up_z             = 0.0
cam0_plane_width      = $width
cam0_plane_height     = $width
cam0_camera_rank      = 0

<time>
cfl_number = 0.3
tlim       = 1.0e30
integrator = rk2
xorder     = 2
ncycle_out = 1

EOF
    # the run's physics; no ray tracing, outflow boundaries
    awk '/^</{keep=($0~/^<(mesh|meshblock|hydro|photchem|photchem_gow17|units|dust)>/)} keep' \
      $RUN/athinput.runtime | sed -e 's/^\([io]x[123]_bc *=\).*/\1 outflow/'
    cat <<EOF

<rayt_point>
rayt_point = false

<rayt_diffuse>
rayt_diffuse = false

<problem>
basedir    = $RUN
snap_id    = HII
snap_outid = out2
num_start  = $N0
num_end    = $N1
num_step   = $DN
read_rad   = false

<rayt_volren>
scheme      = physical
emit[0]     = halpha
emit[1]     = oiii5007
emit[2]     = sii6716
emit[3]     = sii6731
z_dust      = $zd
dust_model  = wd01
dust_file   = $T/inputs/tables/dust/kext_albedo_WD_MW_3.1A_60_D03.dat
levels_dir  = $T/inputs/tables/chianti_v11/levels
line_method = table
EOF
  } > $dir/athinput.lines
  run_logged $dir $BUILD mpi_launch 8 $BUILD/athena -i athinput.lines photchem/evolve=false
done
