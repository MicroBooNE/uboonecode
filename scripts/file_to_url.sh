#! /bin/bash
#----------------------------------------------------------------------
#
# Name: file_to_url.sh
#
# Purpose: Convert a filesystem path to an xrootd url.
#
# Usage:
#
# file_to_url.sh [-h|--help]
# file_to_url.sh <file1> <file2> ...
# file_to_url.sh @<filelist>
# cat <filelist> | file_to_url.sh
#
# Notes:
#
# 1.  Pnfs paths are converted to an xrootd url as
#     /pnfs/... -> root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/...
#
# 2.  Other paths (including local and bluearc) are left unchanged.
#
# 3.  Converted files written to standard output, one per line.
#
#----------------------------------------------------------------------

# Parse arguments.

list=''
nfile=0
while [ $# -gt 0 ]; do
  case "$1" in

    # Help

    -h|--help )
      echo "Usage:"
      echo "file_to_url.sh [-h|--help]"
      echo "file_to_url.sh <file1> <file2> ..."
      echo "file_to_url.sh @<filelist>"
      echo "cat <filelist> | file_to_url.sh"
      exit
    ;;

    # Other unknown options.

    -* )
      echo "Unknown option $1"
      exit 1
    ;;

    # File list.

    @* )
      list=`echo $1 | cut -c2-`
      if [ ! -f $list ]; then
        echo "File list $list not found."
      else
        cat $list | while read file
        do
          echo $file | sed 's;^/pnfs/;root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/;'
          nfile=$(( $nfile + 1 ))
        done
      fi
    ;;

    # Interpret anything as a file name.

    * )
      file=$1
      echo $file | sed 's;^/pnfs/;root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/;'
      nfile=$(( $nfile + 1 ))
    ;;
  esac
  shift
done

# Maybe read files from standard input.

if [ $nfile -eq 0 ]; then
  while read file
  do
    echo $file | sed 's;^/pnfs/;root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/;'
    nfile=$(( $nfile + 1 ))
  done
fi
