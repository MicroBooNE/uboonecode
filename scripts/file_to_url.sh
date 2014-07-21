#! /bin/bash
#----------------------------------------------------------------------
#
# Name: file_to_url.sh
#
# Purpose: Convert a filesystem path to an xrootd url.
#
# Usage:
#
# file_to_url.sh <path>
#
# Notes:
#
# 1.  Pnfs paths are converted to an xrootd url as
#     /pnfs/... -> root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/...
#
# 2.  Other paths (including local and bluearc) are left unchanged.
#
#----------------------------------------------------------------------

# Parse arguments.

path=''
if [ $# -gt 0 ]; then
  path=$1
fi
if [ x$path = x ]; then
  echo "No file specified."
  exit 1
fi
echo $path | sed 's;^/pnfs/;root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/;'
