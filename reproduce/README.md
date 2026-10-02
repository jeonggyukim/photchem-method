# reproduce

Inputs, harnesses and scripts that regenerate each number and figure the paper
quotes, grouped by the section they serve. Code revision: Tigris branch `rayt-photchem-updates` (was `gow17`)
at the commit named in each subdirectory's entry below.

The C++ harnesses include Tigris headers directly and build with

    OMPI_CXX=/opt/homebrew/bin/g++-16 mpicxx -std=c++17 -O2 -ffp-contract=off \
        -I $TIGRIS_DIR/src <harness>.cpp -o <harness>

after `configure.py -gow17 [...]` has written `defs.hpp`. Bit-identity checks
need `-ffp-contract=off`. Python scripts read Tigris output with
`vis/python/athena_read.py` of the checkout in `TIGRIS_DIR`.

The harnesses find code checkouts and simulation output through the environment
variables below; none has a default, and a harness that needs an unset one stops
with a message naming it. The simulation output is not in the repository. The
figure scripts under `figures/` that plot simulation output also read a reduced
copy of the arrays they plot (`reduced*.txt` beside the script), so they remake
the figures without the raw runs.

## Environment variables

| variable | points to |
| --- | --- |
| `TIGRIS_DIR` | a Tigris checkout, branch `rayt-photchem-updates` (headers, `inputs/tables`, `tst/regression/data`, `vis/python`) |
| `PHOTCHEM_RUNS` | the directory holding the Tigris runs: `T3_onezone`, `T6_multi_ion`, `M5_hii_dtype_ions`, `M6_rad_snr`, `T7_cost` |
| `PHOTCHEM_POSTPROC_RUNS` | the directory holding the post-processing runs (`be_net` of F15) |
| `ATHENAPP_PDR_DIR` | an Athena++ checkout with `tst/regression/data/chem_pdr_static.vtk` (`pdr_slab/`) |
| `ATHENAK_DIR` | an AthenaK checkout with the chemistry module (`onezone/`) |
| `PYATHENA_DIR` | a pyathena checkout, for its `data/` directory (`ions/`) |
| `TIGRIS_DATA_DIR` | a checkout of jeonggyukim/tigris-data (dust models in `pdr_sphere/athinput.cloud_postproc_dust`) |
| `CLOUDY_EXE`, `CLOUDY_DATA_DIR` | the Cloudy 25.00 executable and its data directory (`cloudy_stromgren/run_cloudy.sh`) |
| `XUVTOP` | the CHIANTI 11 database (`ions/ion_cooling/build_ion_cooling.py`) |

Python reads them through `paths.py`, C++ through `paths.hpp`, and the shell
scripts as `${TIGRIS_DIR:?...}`. Athena++ does not expand variables in an input
file, so the `athinput.*` files here write `${TIGRIS_DIR}` where a table path
goes; expand them before a run, e.g.
`envsubst '$TIGRIS_DIR $TIGRIS_DATA_DIR' < athinput.X > athinput.local`.

## Running a harness

    pip install ".[reproduce]"          # from the repository root
    export TIGRIS_DIR=/path/to/tigris PHOTCHEM_RUNS=/path/to/runs   # as the harness needs
    python reproduce/figures/F12_hii_dtype/plot.py

Scripts are run in place, from any directory. A figure script with a `reduced*.txt`
reads it and needs no variable; `--from-runs` makes it read the raw runs and
rewrite the reduced file. The figures go to `figures/` at the repository root.

| directory | serves | what |
| --- | --- | --- |
| `onezone/` | 5.1 | B1: AthenaK GOW17_uniform in one zone, Tigris port vs AthenaK and Athena++ |
| `solver/` | 4, 5.1 | Riccati step check, ionization-front one-zone error vs substeps, change-based step control |
| `pdr_slab/` | 5.2 | PDR slab vs Athena++ chem_pdr_static; CO shielding; per-reaction dust |
| `ions/` | 2.4 | ion-ladder balance, band-averaged ion cross sections, CIE tables and pool |
| `cloudy_stromgren/` | 5.3 | Cloudy 25.00 static Stromgren sphere: input, SED conversion, reduced radial profiles |
| `hii_dtype/` | 5.4 | D-type H II region: input, shell radius vs NCR, diagnostic figure and movie |
| `rad_snr/` | 5.5 | radiative SNR: input, shell-formation mass and momentum vs NCR |

Code revisions: all copied at Tigris `gow17` 2026-09-24, c4e12eed1, after the
`--gow17_bands` per-band rates commit.

## cloudy_stromgren

`CLOUDY_DATA_PATH=<cloudy>/data:$PWD cloudy.exe -r hii` (the SED file must be on
the data path). `make_sed.py` converts the SB99 table to Cloudy's `table SED`
format and checks the photon ratios; `make_profiles.py` reduces the Cloudy
outputs to `radial_profiles.txt` and prints the H+-zone averages.

## bands

`build.sh` builds -gow17 --gow17_ions=O2,S2 --gow17_bands=5 (Tigris c4e12eed1);
`mpirun -np 1 ./athena -i athinput.hii_bands5 time/nlim=0` writes the computed
band averages into athinput.runtime; `check_photon_frac.py athinput.runtime`
compares the ionizing photon fractions with the SED ratios the Cloudy run uses
(section 3; agreement 6e-4 and 7e-4 relative).
`run_cloudy.sh` runs all three steps. The Tigris side is `athinput.static_bands5`
(the hii pgen with `hydro/active = fixed`, no dust, no ISRF, no CRs) run as
`mpirun -np 8 athena -i athinput.static_bands5 mesh/x1min=-6 mesh/x1max=6
mesh/x2min=-6 mesh/x2max=6 mesh/x3min=-6 mesh/x3max=6 time/tlim=0.05 time/nlim=500
output2/dt=0.005` with the 5-band build of `bands/build.sh`; `compare_cloudy.py
OUT.png "5 bands:RUNDIR"` draws the radial profiles and prints the H+-zone averages.

## Rule

Every script, harness or input that produced a number or a figure, including
debug harnesses, is copied here when it is made and committed. The working
folder (`$PHOTCHEM_RUNS`) holds runs and
scratch; this repository is the record. `ncr_physics/` holds the checks of the
NCR physics ported into GOW17 (hot gas, LyC path, equilibrium solver, limiter).
