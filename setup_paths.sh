#!/usr/bin/env bash

# Source this file to use util-library without installing a Python package:
#
#   source /path/to/util-library/setup_paths.sh

util_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

case ":${PYTHONPATH:-}:" in
  *":${util_root}/python:"*) ;;
  *) export PYTHONPATH="${util_root}/python${PYTHONPATH:+:${PYTHONPATH}}" ;;
esac

case ":${PATH}:" in
  *":${util_root}/scripts:"*) ;;
  *) export PATH="${util_root}/scripts:${PATH}" ;;
esac

export UTIL_LIBRARY_ROOT="${util_root}"
echo "util-library paths enabled: ${util_root}"
