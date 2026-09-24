# reproduce

Inputs, harnesses and scripts that regenerate each number and figure the paper
quotes, grouped by the section they serve. Code revision: Tigris branch `gow17`
at the commit named in each subdirectory's entry below.

The C++ harnesses include Tigris headers directly and build with

    OMPI_CXX=/opt/homebrew/bin/g++-16 mpicxx -std=c++17 -O2 -ffp-contract=off \
        -I ~/Projects/tigris-gow17/src <harness>.cpp -o <harness>

after `configure.py -gow17 [...]` has written `defs.hpp`. Bit-identity checks
need `-ffp-contract=off`. Python scripts use the `pyathena` conda env and read
Tigris output with `vis/python/athena_read.py`. Paths inside the scripts are
absolute and point at the author's checkout; that is to be replaced by
arguments before the repository is shared.

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
