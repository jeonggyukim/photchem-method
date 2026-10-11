# photchem-method

This repository has the scripts behind the figures of the TIGRESS++
photochemistry method paper. It also has the test problems those figures come
from, set up to run with Tigris and AthenaK.

- `schematics/` draws the diagrams: how the modules couple, the ray tracer, and
  the chemical network.
- `tests/` has one folder per test problem. A paper figure lives in
  `tests/<problem>/fig_<name>/`, next to the data it plots.
- `photchem_method/` is a small Python package the scripts share.

You do not need the simulation runs to remake a figure. Each figure that comes
from a simulation keeps a small `reduced*.txt` file with the numbers it plots.

## Installation

You need Python 3.10 or newer. Clone the repository first:

```sh
git clone https://github.com/jeonggyukim/photchem-method.git
cd photchem-method
```

With conda:

```sh
conda create -n photchem-method -c conda-forge python=3.12 numpy matplotlib healpy
conda activate photchem-method
pip install -e . --no-deps
```

Or with [uv](https://docs.astral.sh/uv/):

```sh
uv venv --python 3.12
source .venv/bin/activate
uv pip install -e .
```

A few scripts need more packages. `pip install -e ".[slides]"` adds python-pptx
for `schematics/make_cr_slides.py`, and `pip install -e ".[tests]"` adds pandas
and pyathena for the scripts under `tests/` that need them.

## Making a figure

Run a script from anywhere. It writes its figure to `figures/`.

```sh
python schematics/rad_cr_chem_coupling.py paper --mhd --cr-left --crpic=B
python tests/hii_dtype/fig_hii_dtype/plot.py
```

The comment at the top of each script lists its options and the files it writes.

## Rerunning the tests

To rebuild Tigris or AthenaK, rerun a test, or redraw a figure from the raw runs
(`plot.py --from-runs`), see `tests/README.md`. It explains how to set up your
machine and which environment variables point to the codes and runs.

## Keeping things consistent

`python -m photchem_method.check` looks for common mistakes, such as a hard-coded
path to one person's computer or a data file that does not say where it came from.
To run it before every commit, turn on the hook once:

```sh
git config core.hooksPath .githooks
```

GitHub runs the same check, and every figure script, on each push to `main`.

Figures can differ by a few pixels between matplotlib versions, because text is
placed slightly differently. The data are the same.
