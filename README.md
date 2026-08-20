# util-library

Routines for bfield calculation, profile fitting, EFIT/geqdsk processing, and
basic geometry utilities.

## MATLAB

Add the MATLAB bfield library to your path:

```matlab
addpath(genpath('/path/to/util-library/matlab/bfield_library_jdl'))
```

## Python

Install the Python utilities in editable mode:

```bash
cd /path/to/util-library
python -m pip install -e .
```

This installs the top-level modules:

```python
import EFITutils
import GeoUtils
```

MPEX-specific Python routines also need the private `MPEX-modeling-data` repo:

```bash
source /path/to/MPEX-modeling-data/setup_paths.sh
```

and command-line helpers such as:

```bash
gfile_vessel
gfile_update
gfile_vertical_mirror
gfile_symmetrize
plot_gfile_transform_tests
plot_ogr
order_ogr
refine_ogr_elements
convert_ogr_to_structure_dat
double_gfile_resolution
```

`gfile_update` updates selected GEQDSK fields from an EQU file, `rzpsi.dat`, or
an OGR vessel file. `gfile_vertical_mirror` flips an equilibrium about a
specified Z plane. `gfile_symmetrize` builds a vertically symmetric gfile from
the upper or lower side. `plot_gfile_transform_tests` writes before/after PNGs
for visual inspection of these transformations; add `--symmetrize-vessel` to
include the limiter/vessel in the symmetrization example plots.

Most commands only need numpy. Plotting commands, such as `plot_ogr` and
`plot_gfile_transform_tests`, also need matplotlib. If matplotlib is not
already available:

```bash
python -m pip install matplotlib
```

or install the optional plotting dependency with the repo:

```bash
python -m pip install -e '.[plot]'
```

## Attribution

Developed and maintained by J.D. Lore.
