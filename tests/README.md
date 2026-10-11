# tests

The paper's test problems: inputs, harnesses and scripts that regenerate each
number and figure the paper quotes, grouped by the section they serve. Code revision: Tigris branch `rayt-photchem-updates` (was `gow17`)
at the commit named in each subdirectory's entry below.

The C++ harnesses include Tigris headers directly and build with

    source tests/env.sh     # sets OMPI_CXX and CXX from tests/machine.sh
    mpicxx -std=c++17 -O2 -ffp-contract=off -I $TIGRIS_DIR/src <harness>.cpp -o <harness>

after `configure.py -gow17 [...]` has written `defs.hpp`. Bit-identity checks
need `-ffp-contract=off`. Python scripts read Tigris output with
`vis/python/athena_read.py` of the checkout in `TIGRIS_DIR`.

The harnesses find code checkouts and simulation output through the environment
variables below; none has a default, and a harness that needs an unset one stops
with a message naming it. The simulation output is not in the repository. Each
paper figure has a `fig_<name>/` folder inside the test it plots, holding its
`plot.py` and, where the figure comes from simulation output, a reduced copy of
the arrays it plots (`reduced*.txt`), so it remakes the figure without the raw
runs.

## Environment variables

| variable | points to |
| --- | --- |
| `TIGRIS_DIR` | a Tigris checkout, branch `rayt-photchem-updates` (headers, `inputs/tables`, `tst/regression/data`, `vis/python`) |
| `PHOTCHEM_RUNS` | the directory holding the Tigris runs: `T3_onezone`, `T6_multi_ion`, `M5_hii_dtype_ions`, `M6_rad_snr`, `T7_cost` |
| `PHOTCHEM_POSTPROC_RUNS` | the directory holding the post-processing runs (`be_net`, `be_sphere/`) |
| `ATHENAPP_PDR_DIR` | an Athena++ checkout with `tst/regression/data/chem_pdr_static.vtk` (`pdr_slab/`) |
| `ATHENAK_DIR` | an AthenaK checkout with the chemistry module (`onezone/`) |
| `PYATHENA_DIR` | a pyathena checkout, for its `data/` directory (`ions/`) |
| `TIGRIS_DATA_DIR` | a checkout of jeonggyukim/tigris-data (dust models in `pdr_sphere/athinput.cloud_postproc_dust`) |
| `CLOUDY_EXE`, `CLOUDY_DATA_DIR` | the Cloudy 25.00 executable and its data directory (`cloudy_stromgren/run_cloudy.sh`) |
| `XUVTOP` | the CHIANTI 11 database (`ions/ion_cooling/build_ion_cooling.py`, which also needs ChiantiPy 0.16; `pip install ChiantiPy`) |

Python reads them through `photchem_method/paths.py`, C++ through `paths.hpp`, and the shell
scripts as `${TIGRIS_DIR:?...}`. Athena++ does not expand variables in an input
file, so the `athinput.*` files here write `${TIGRIS_DIR}` where a table path
goes; expand them before a run, e.g.
`envsubst '$TIGRIS_DIR $TIGRIS_DATA_DIR' < athinput.X > athinput.local`.

## Running a harness

    pip install -e ".[tests]"           # from the repository root
    export TIGRIS_DIR=/path/to/tigris PHOTCHEM_RUNS=/path/to/runs   # as the harness needs
    python tests/hii_dtype/fig_hii_dtype/plot.py

Scripts are run in place, from any directory. A figure script with a `reduced*.txt`
reads it and needs no variable; `--from-runs` makes it read the raw runs and
rewrite the reduced file. The figures go to `figures/` at the repository root.

| directory | serves | what |
| --- | --- | --- |
| `onezone/` | 5.1 | B1: AthenaK GOW17_uniform in one zone, Tigris port vs AthenaK and Athena++ |
| `ode_solver/` | 4, 5.1 | time integration of the per-cell ODEs: semi-implicit error vs substeps, Gauss-Seidel order and blocks, Riccati step, step control |
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

## Machines

Build and run scripts take every machine-specific setting (compilers, MPI launcher,
library prefixes, Python) from `tests/machine.sh`, which is not tracked. Make it once
per machine from the preset in `tests/machines/`:

    cp tests/machines/grammar.sh tests/machine.sh     # or macos-homebrew, syntax, stellar

and edit it if your machine differs. A new machine gets a new preset in
`tests/machines/`, defining the names listed at the top of `macos-homebrew.sh`.
`stellar.sh` is incomplete: its module loads still have to be filled in.

## Conventions

These rules let anyone with a Tigris or AthenaK checkout rerun a result and know which
code produced it. `python -m photchem_method.check` tests them, and the pre-commit hook
runs that check, so a commit that breaks a rule is refused. Enable the hook once per
clone with `git config core.hooksPath .githooks`.

1. **One folder per result.** A test lives in `tests/<problem>/`; a paper figure from
   it in `tests/<problem>/fig_<name>/`, named by content, never by figure number.
   Every `plot*.py` under `tests/` is in a `fig_<name>/` folder.
2. **Builds and runs go through `tests/env.sh`.** A build script sources it and calls
   `tigris_build OUTDIR <configure flags>`, which writes `OUTDIR/build.log`: the code's
   repository, branch and commit, whether the tree was dirty, the configure line and
   the compiler. A dirty build also saves `OUTDIR/code.diff` and lists untracked
   files under `src/`, so the commit plus the diff is the source that was built. A run
   goes through `run_logged RUNDIR BUILDDIR mpi_launch N ...`, which copies
   `build.log` into `RUNDIR/run.log` and adds the command line, host, exit status and
   wall time.
3. **No machine paths in tracked files.** Only `tests/machines/*.sh` name a compiler,
   library prefix, Python or home directory.
4. **Every location comes from an environment variable** listed in
   `photchem_method/paths.py` (the table above). Python reads them through
   `photchem_method.paths`, C++ through `tests/paths.hpp`, shell as `${NAME:?...}`.
5. **A figure reads its reduced file by default**; `--from-runs` reads the runs and
   rewrites it. `reduced.save` requires the list of runs and files the arrays came
   from and writes them, with their commits, as `# provenance:` lines at the top of
   the file.
6. **Results of other codes are read, not copied.** A reference value from Athena++,
   AthenaK, Cloudy or CHIANTI is read from that code's checkout or output at
   `--from-runs` time and named in the provenance lines; it is never typed into a
   script.
7. **Everything that produced a number or a figure is committed here** when it is
   made, including checks written while debugging. `$PHOTCHEM_RUNS` holds the runs;
   this repository is the record.
