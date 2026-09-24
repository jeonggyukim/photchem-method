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
