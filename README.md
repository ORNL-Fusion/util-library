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
plot_ogr
order_ogr
refine_ogr_elements
convert_ogr_to_structure_dat
double_gfile_resolution
```

Most commands only need numpy.  Plotting commands, such as `plot_ogr`, also
need matplotlib.  If matplotlib is not already available:

```bash
python -m pip install matplotlib
```

or install the optional plotting dependency with the repo:

```bash
python -m pip install -e '.[plot]'
```

## Attribution

Developed and maintained by J.D. Lore.
