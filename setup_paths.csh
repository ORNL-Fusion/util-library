# Source this file to use util-library without installing a Python package.
#
# From the util-library root:
#   source setup_paths.csh
#
# Or from elsewhere:
#   source /path/to/util-library/setup_paths.csh /path/to/util-library

if ($#argv >= 1) then
  set util_root = "$argv[1]"
else
  set util_root = "$cwd"
endif

if ($?PYTHONPATH) then
  setenv PYTHONPATH "${util_root}/python:${PYTHONPATH}"
else
  setenv PYTHONPATH "${util_root}/python"
endif

setenv PATH "${util_root}/scripts:${PATH}"
setenv UTIL_LIBRARY_ROOT "${util_root}"
echo "util-library paths enabled: ${util_root}"
