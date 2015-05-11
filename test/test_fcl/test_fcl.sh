#! /bin/bash

# Loop over all installed fcl files.

find $MRB_BUILDDIR/uboonecode/job -name \*.fcl -print | while read fcl
do
  echo "Testing fcl file $fcl"

  # Parse this fcl file.

  out=`basename ${fcl}`.out
  lar -c $fcl --debug-config $out

  # Exit status 1 counts as success

  stat=$?
  if [ $stat -ne 0 -a $stat -ne 1 ]; then
    exit $stat
  fi
  
done
