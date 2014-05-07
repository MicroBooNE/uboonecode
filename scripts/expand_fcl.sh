#! /bin/bash
#------------------------------------------------------------------
#
# Purpose: Prints a canonicalized expansion of an fcl configuration
#          file.  
#
# Usage: expand_fcl.sh <fcl-file>
#
# Usage notes.
#
# 1.  The larsoft environment must be initialized, in particular:
#
#     a) The lar executable must be on the execution path.
#
#     b) The fcl file search path $FHICL_FILE_PATH must be
#        properly initialized.
# 
# Created: H. Greenlee, 30-Jul-2013
#
#------------------------------------------------------------------

# Parse arguments.

fcl=''
if [ $# -gt 0 ]; then
  fcl=$1
fi

if [ x$fcl = x ]; then
  echo "Usage: expand_fcl.sh <fcl-file>"
  exit 1
fi

# Use lar executable to expand the fcl script.
# Get rid of art boilerplate and redirect to standard output.

export ART_DEBUG_CONFIG=1
lar -c $fcl 2>&1 >/dev/null | sed -n '2,$p'

