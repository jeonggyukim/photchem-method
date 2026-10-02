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

## The TIGRESS++ coupling schematic

`scripts/rad_cr_chem_coupling.py` draws how radiation transfer, cosmic-ray
transport, photochemistry and gas dynamics exchange quantities in TIGRESS++. It
needs Python 3.10 or later with numpy, matplotlib and healpy (no LaTeX), and the
three modules beside it: `adaptive_rays_3d.py` (point-source ray trees, on
HEALPix directions), `diffuse_rays_3d.py` (diffuse-field rays) and
`clumpy_cloud.py` (the cloud both ray panels cross). `pyproject.toml` declares
these; the scripts are run in place, not imported as a package.

```
git clone https://github.com/jeonggyukim/photchem-method.git
cd photchem-method
python -m venv .venv && source .venv/bin/activate   # or any environment
pip install .
cd scripts
python rad_cr_chem_coupling.py 3d2 paper --mhd --cr-left --crpic=B
```

`pip install ".[slides]"` adds python-pptx for `make_cr_slides.py`, and
`pip install ".[reproduce]"` adds what the `reproduce/` harnesses read (pandas and
pyathena; some also need `athena_read.py` from the Athena++ or Tigris
`vis/python` directory).

The output, `rad_cr_chem_coupling_3d2_paper_mhd_crleft_crB.png` and `.pdf`, goes to
`figures/` (created if absent). The options:

| option | effect |
| --- | --- |
| `3d2` / `3d` / `2d` | both ray panels in 3D / point sources only in 3D / flat |
| `paper` / `talk` | symbols with a few words each / symbols in large type for slides |
| `--mhd` | adds the gas-dynamics box of a coupled simulation |
| `--cr-left` | cosmic-ray box on the left, radiation on the right |
| `--crpic=A` / `--crpic=B` | cosmic-ray picture: resonant scattering / the self-confinement loop |

The two variants used most:

```
python rad_cr_chem_coupling.py 3d2 paper --mhd --cr-left --crpic=B   # paper style
python rad_cr_chem_coupling.py 3d2 talk  --mhd --cr-left --crpic=B   # slides
```

Labels, colours and box positions are set near the top of the script (`LABELS`,
`POINT_SOURCES`, `make_cloud`).
