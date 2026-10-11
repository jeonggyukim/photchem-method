# tests

This directory holds tests of the TIGRESS++ photochemistry module in Tigris and
AthenaK, and the scripts that make figures from them. Each test has its inputs, its
check programs and the scripts that produce its results. Most tests run on Tigris, a
code based on Athena++. GOW17 is the chemical network of Gong, Ostriker & Wolfire
(2017). NCR is the photochemistry and cooling of TIGRESS-NCR (Kim J.-G. et al. 2023).

Every result records the code revision it came from. A build writes the code commit
to `build.log`. Every reduced data file names its sources and their commits in
`# provenance:` lines. The rules are under Conventions below.

| directory | contents |
| --- | --- |
| `bands/` | Tigris builds with 5 and 6 point-source bands, the check of the ionizing photon fractions, the band-count figure |
| `be_sphere/` | a Bonnor-Ebert sphere in an isotropic far-ultraviolet field, post-processed with NCR and with GOW17 and ions |
| `cloudy_stromgren/` | a static Stromgren sphere against Cloudy 25.00: Cloudy and Tigris inputs, the Starburst99 spectrum conversion, comparison figures |
| `cost/` | the cost per cycle of each thermal and chemical mode, split by task |
| `dust_temperature/` | the closed-form dust temperature against the grain equilibrium and the gas-dust balance |
| `hii_dtype/` | a D-type H II region: input, shell radius against NCR, diagnostic figure and movie |
| `ions/` | ion-ladder balance, band-averaged ion cross sections, tables of collisional ionization equilibrium (CIE) and the higher ions |
| `ncr_physics/` | NCR physics inside GOW17: NCR equilibrium method, hot-gas cooling and electrons, grain recombination in hot gas, substep limiter, Lyman-continuum photoionization |
| `ode_solver/` | time integration of the chemical equations of each cell: semi-implicit error against substeps, Gauss-Seidel order and blocks, Riccati step, step control |
| `onezone/` | test B1, the AthenaK GOW17_uniform problem in one zone: the Tigris copy against AthenaK and Athena++; the equilibrium figure |
| `pdr_slab/` | a photodissociation region (PDR) slab against Athena++ chem_pdr_static; CO shielding; dust in each reaction |
| `pdr_sphere/` | a uniform sphere lit from every side, post-processed with NCR and with GOW17 |
| `rad_snr/` | a radiative supernova remnant: shell-formation mass and momentum against NCR, and the ions of the shell |

## Building a harness

A harness is a small C++ program that includes Tigris headers and runs one check.
First run `configure.py -gow17 [...]` in the Tigris checkout. It writes `defs.hpp`.
Then build the harness:

    source tests/env.sh     # sets OMPI_CXX and CXX from tests/machine.sh
    mpicxx -std=c++17 -O2 -ffp-contract=off -I $TIGRIS_DIR/src <harness>.cpp -o <harness>

The flag `-ffp-contract=off` stops the compiler from fusing a multiply and an add
into one instruction. Checks that compare two codes bit for bit need it.

## Remaking a figure

    pip install -e .                    # from the repository root
    python tests/hii_dtype/fig_hii_dtype/plot.py

A script runs from any directory. Each figure has a `fig_<name>/` folder inside its
test. The folder holds `plot.py`. For a figure from simulation output it also holds
a reduced copy of the plotted arrays (`reduced*.txt`). The script reads the reduced
file and needs no environment variable. With `--from-runs` it reads the raw runs and
rewrites the reduced file. Figures go to `figures/` at the repository root. A script
writes a PNG only. Set `PHOTCHEM_FIGURE_FORMATS=png,pdf` to also write a PDF.

The raw simulation output is not in the repository.
`tests/pdr_sphere/fig_slices/plot.py` reads its reduced slices from
`$PHOTCHEM_DATA_DIR/pdr_sphere`, outside the repository.

`pip install -e ".[tests]"` adds pandas and pyathena. Only
`ions/pyathena_rate_check/compare_rates.py` needs them.

## Environment variables

The harnesses find code checkouts and simulation output through these variables.
Plotting a reduced file needs none of them. `--from-runs`, the harnesses and building
and running the codes need them. None has a default. A script that needs an unset
variable stops with a message naming it.

| variable | points to |
| --- | --- |
| `TIGRIS_DIR` | a Tigris checkout, branch `rayt-photchem-updates` (headers, `inputs/tables`, `tst/regression/data`, `vis/python`) |
| `PHOTCHEM_RUNS` | the directory of the Tigris runs: `T3_onezone`, `T6_multi_ion`, `M5_hii_dtype_ions`, `M6_rad_snr`, `T7_cost` |
| `PHOTCHEM_POSTPROC_RUNS` | the directory of the post-processing runs: `be_net` (`be_sphere/`) |
| `PHOTCHEM_DATA_DIR` | the directory holding the large reduced arrays that are not tracked: `pdr_sphere/` (`pdr_sphere/run.sh`, `pdr_sphere/fig_slices/`) |
| `ATHENAPP_PDR_DIR` | an Athena++ checkout with `tst/regression/data/chem_pdr_static.vtk` (`pdr_slab/`) |
| `ATHENAK_DIR` | an AthenaK checkout with the chemistry module (`onezone/`) |
| `PYATHENA_DIR` | a pyathena checkout, for its `data/` directory (`ions/`) |
| `TIGRIS_DATA_DIR` | a checkout of https://github.com/jeonggyukim/tigris-data, for the dust models (`pdr_sphere/athinput.cloud_postproc_dust`) |
| `CLOUDY_EXE`, `CLOUDY_DATA_DIR` | the Cloudy 25.00 executable and its data directory (`cloudy_stromgren/run_cloudy.sh`) |
| `XUVTOP` | the CHIANTI 11 database, read by `ions/ion_cooling/build_ion_cooling.py` through ChiantiPy 0.16 (`pip install ChiantiPy`) |

Athena++ does not expand variables in an input file. Where a table path goes, the
`athinput.*` files here write `${TIGRIS_DIR}`. Expand them before a run, for example
with `envsubst '$TIGRIS_DIR $TIGRIS_DATA_DIR' < athinput.X > athinput.local`.

## cloudy_stromgren

`run_cloudy.sh` runs the Cloudy side in three steps. `make_sed.py` converts the
Starburst99 spectrum of Tigris to the `table SED` format of Cloudy and checks the
photon ratios. Cloudy then reads `hii.in`. `make_profiles.py` reduces the Cloudy
output to `radial_profiles.txt` and prints the averages over the H+ zone.

The Tigris side is `athinput.static_bands5`, the hii problem generator with
`hydro/active = fixed`. Run it with a 5-band Tigris build, such as one of the
`bands/build_bands5_*.sh` builds:

    mpirun -np 8 athena -i athinput.static_bands5 mesh/x1min=-6 mesh/x1max=6 \
        mesh/x2min=-6 mesh/x2max=6 mesh/x3min=-6 mesh/x3max=6 \
        time/tlim=0.05 time/nlim=500 output2/dt=0.005

`compare_cloudy.py OUT.png "5 bands:RUNDIR"` draws the radial profiles against
Cloudy and prints the averages over the H+ zone.

## bands

Each `build_bands*.sh` script builds Tigris with GOW17, the ions in its name and 5
or 6 bands. For example, `build_bands5_O2S2.sh` uses `--gow17_ions=O2,S2
--gow17_bands=5`. The run

    mpirun -np 1 ./athena -i athinput.hii_bands5 time/nlim=0

writes the band averages that Tigris computes from the spectrum into
`athinput.runtime`. `check_photon_frac.py athinput.runtime` compares the photon
fractions above 24.59 eV and above 35.12 eV with the ratios of the spectrum that
Cloudy reads. The relative differences are 6e-4 and 7e-4.
`fig_band_count/run_series.sh` runs the post-processed sphere with 3, 5, 6 and 7
bands for the band-count figure.

## Machines

Build and run scripts take every machine-specific setting (compilers, MPI launcher,
library prefixes, Python) from `tests/machine.sh`. That file is not tracked. Make it
once per machine from a preset in `tests/machines/`:

    cp tests/machines/grammar.sh tests/machine.sh     # or macos-homebrew, syntax, stellar

Edit it if your machine differs. A new machine gets a new preset in
`tests/machines/` that defines the names listed at the top of `macos-homebrew.sh`.
`stellar.sh` is incomplete: its module loads are not filled in yet.

## Conventions

These rules let anyone with a Tigris or AthenaK checkout rerun a result and know which
code produced it. `python -m photchem_method.check` tests them. The pre-commit hook
runs that check and refuses a commit that breaks a rule. Enable the hook once per
clone with `git config core.hooksPath .githooks`.

1. **One folder per result.** A test lives in `tests/<problem>/`. A figure from it
   lives in `tests/<problem>/fig_<name>/`, named by content, never by figure number.
   Every `plot*.py` under `tests/` is in a `fig_<name>/` folder.
2. **Builds and runs go through `tests/env.sh`.** A build script sources it and calls
   `tigris_build OUTDIR <configure flags>`. That writes `OUTDIR/build.log` with the
   code's repository, branch and commit, whether the tree was dirty, the configure
   line and the compiler. A dirty build also saves `OUTDIR/code.diff` and lists the
   untracked files under `src/`, so the commit plus the diff is the source that was
   built. A run goes through `run_logged RUNDIR BUILDDIR mpi_launch N ...`. It copies
   `build.log` into `RUNDIR/run.log` and adds the command line, host, exit status and
   wall time.
3. **No machine paths in tracked files.** Only `tests/machines/*.sh` name a compiler,
   library prefix, Python or home directory.
4. **Every location comes from an environment variable** listed in
   `photchem_method/paths.py` (the table above). Python reads them through
   `photchem_method.paths`, C++ through `tests/paths.hpp`, and shell scripts as
   `${NAME:?...}`.
5. **A figure reads its reduced file by default.** `--from-runs` reads the runs and
   rewrites it. `reduced.save` requires the list of runs and files the arrays came
   from. It writes them, with their commits, as `# provenance:` lines at the top of
   the file. A figure script creates its output directory before it saves.
6. **Results of other codes are read, not copied.** A reference value from Athena++,
   AthenaK, Cloudy or CHIANTI is read from that code's checkout or output at
   `--from-runs` time and named in the provenance lines. It is never typed into a
   script.
7. **Everything that produced a number or a figure is committed here** when it is
   made, including checks written while debugging. `$PHOTCHEM_RUNS` holds the runs.
   This repository is the record.
