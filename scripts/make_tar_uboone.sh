#! /bin/bash
#------------------------------------------------------------------
#
# Name: make_tar_uboone.sh
#
# Purpose: Make a tarball for a larsoft test release, for purpose
#          of being shipped to grid worker.
#
# Usage:
#
# make_tar_uboone.sh [-h|--help] [-d dir] <tarball-name>
#
# Options:
#
# -h, --help - Print help message.
# -d dir     - Specify location of test release (default $SRT_PRIVATE_CONTEXT
#              or $MRB_INSTALL).
#
# Notes.
#
# 1.  Resulting tar file is always created in the current directory.
#
# 2.  This scripts attempts to exclude non-relevant files from the
#     tarball, including the following.
#     a) The tmp subdirectory.
#     b) .svn directories.
#     c) *.root files in the top directory.
#
# 3.  This script can be invoked in the test release directory, or in
#     another directory.
#
# Created 8-Nov-2013  H. Greenlee
#
#------------------------------------------------------------------

# Help function.

function dohelp {
  echo "Usage: make_tar_uboone.sh [-h|--help] [-d dir] <tarball-name>"
}

# Parse arguments.

if [ $# -eq 0 ]; then
  dohelp
  exit
fi

# Defaults.

tar=''
dir=''
if [ x$SRT_PRIVATE_CONTEXT != x ]; then
  dir=$SRT_PRIVATE_CONTEXT
fi
if [ x$MRB_INSTALL != x ]; then
  dir=$MRB_INSTALL
fi
if [ x$dir = x ]; then
  dir=`pwd`
fi

while [ $# -gt 0 ]; do
  case "$1" in

    # Help.
    -h|--help )
      dohelp
      exit
      ;;

    # Directory.
    -d )
      if [ $# -gt 1 ]; then
        dir=$2
        shift
      fi
      ;;

    # Other options.
    -* )
      echo "Unrecognized option $1"
      dohelp
      exit
      ;;

    # Positional.
    * )
      if [ x$tar = x ]; then
        tar=$1
      else
        echo "Too many arguments."
        dohelp
        exit 1
      fi

  esac
  shift
done

# Make sure source directory is defined and exists.

if [ x$dir = x ]; then
  echo "No source directory specified."
  dohelp
  exit 1
fi
if [ ! -d $dir ]; then
  echo "Directory $dir doesn't exist."
  exit 1
fi

# Make sure that a tarball was specified.

if [ x$tar = x ]; then
  echo "No tarball specified."
  dohelp
  exit 1
fi

# If the tarball already exists, delete it.

if [ -f $tar ]; then
  rm $tar
fi

ls -A $dir | egrep -v '.root$|tmp' | tar -C $dir -T- -czf $tar --exclude=.svn --exclude=\*.tar



