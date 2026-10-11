# Paper figures

Figures for the photochemistry method paper, and the scripts that make them.

This repository is meant to be shareable. Collaborators may get it, and it may
be published with the paper, so nothing private belongs here and nothing bulky
does either -- raw simulation output and Cloudy runs go in the parent
directory, which Dropbox backs up without committing anything.

| directory | holds |
| --- | --- |
| `scripts/` | one script per figure, named for the figure it makes |
| `reproduce/` | the inputs and harnesses that regenerate each quoted number |
| `figures/` | the output, tracked once a figure is near final |

A script has to run from something a reader can get, so where a figure comes
from a large simulation output, commit the reduced array it actually plots and
leave the raw output in the parent directory.

Copy a figure into the Overleaf `figures/` once it is going into the paper, so
that directory holds only what the manuscript includes.

## Requirement

Python **3.10** or higher, with numpy, matplotlib and healpy. LaTeX is not
needed.

## Installation

The scripts run in place from `scripts/`; nothing is installed as a package.
`pyproject.toml` declares the dependencies, so installing the repository
installs only those.

```sh
git clone https://github.com/jeonggyukim/photchem-method.git
cd photchem-method
```

Then pick one of the options below.

**Option 1 — conda** ([miniforge](https://github.com/conda-forge/miniforge)
recommended; swap `conda` for `mamba` if you have it)
```sh
conda create -n photchem-method -c conda-forge python=3.12 numpy matplotlib healpy
conda activate photchem-method
```

**Option 2 — [uv](https://docs.astral.sh/uv/) (no conda)**
```sh
uv venv --python 3.12
source .venv/bin/activate
uv pip install .
```

Plain `python -m venv .venv` followed by `pip install .` works the same way.

Two extras add what the other scripts read:

```sh
pip install ".[slides]"      # python-pptx, for make_cr_slides.py
pip install ".[reproduce]"   # pandas and pyathena, for the reproduce/ harnesses
```

Some harnesses also need `athena_read.py` from the Athena++ or Tigris
`vis/python` directory, which is not on PyPI. `reproduce/README.md` lists the
environment variables that tell them where the codes and runs are.

Sanity check any option with:
```sh
python scripts/rad_cr_chem_coupling.py --help
```

### Uninstall

```sh
conda deactivate && conda env remove -n photchem-method   # Option 1
deactivate && rm -rf .venv                                # Option 2
```

## Making the figures

Run each script from `scripts/`. Each figure is written as a `.png` and a
`.pdf` of the same name to `figures/`. `rad_cr_chem_coupling.py`,
`ncr_coupling.py`, `fig_network.py` and `fig_network_talk.py` take `--help`.

| script | output in `figures/` |
| --- | --- |
| `rad_cr_chem_coupling.py [talk\|paper] [--mhd] [--cr-left] [--crpic=A\|B]` | `rad_cr_chem_coupling_<style>[_mhd][_crleft][_crA\|_crB]` |
| `ncr_coupling.py [--no-cr-icon]` | `ncr_coupling`, or `ncr_coupling_nocricon` |
| `art_schematic_talk.py` | `art_schematic_talk` |
| `fig_network.py [--outdir DIR]` | `fig_network_core`, `_tigris_greycore`, `_2panel`, `_1panel_colour`, and the contact sheet `fig_network_versions.png` |
| `fig_network_talk.py [--outdir DIR]` | `fig_network_talk_all`, `fig_network_talk_ncr` |
| `adaptive_rays_3d.py [outdir]` | `adaptive_rays_3d`, `adaptive_rays_3d_panel` |
| `diffuse_rays_3d.py [outdir]` | `diffuse_rays_3d` |
| `make_cr_slides.py [output.pptx]` | `cr-coupling.pptx`, from the `rad_cr_chem_coupling` figures |
| `pdr_sphere_slices.py [npz ...]` | `pdr_sphere_<mode>`, from reduced arrays in `../data/pdr_sphere`, beside the repository |

For example:
```sh
cd scripts
python rad_cr_chem_coupling.py paper --mhd --cr-left --crpic=B
```

The paper figures F02–F14 are made by `reproduce/figures/F*/plot.py`.

matplotlib 3.10 and 3.11 place text slightly differently, so a figure is
reproduced pixel for pixel only with the same matplotlib version.
