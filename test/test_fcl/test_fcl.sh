#! /bin/bash

# Loop over all installed fcl files.

find $MRB_BUILDDIR/uboonecode/job -name \*.fcl -print | while read fcl
do
  echo "Testing fcl file $fcl"

  # Parse this fcl file.

  fclout=`basename ${fcl}`.out
  larout=`basename ${fcl}`.lar.out
  larerr=`basename ${fcl}`.lar.err
  lar -c $fcl --debug-config $fclout > $larout 2> $larerr

  # Exit status 1 counts as success.
  # Any other exit status exit immediately.

  stat=$?
  if [ $stat -ne 0 -a $stat -ne 1 ]; then
    echo "Error parsing ${fcl}."
    exit $stat
  fi

  # Flag files that have services.user blocks.

  if grep -q user: $fclout; then
    echo "Deprecated services.user found in ${fcl}."
    exit 1
  fi

  # Check for certain kinds of diagnostic output.

  if egrep -iq 'deprecated|no longer supported' $larerr; then
    echo "Deprecated fcl construct found in ${fcl}."
    exit 1
  fi

  # We consider it an error if the diagnostic output from lar has more than two lines.

  if [ `cat $larerr | wc -l` -gt 2 ]; then
    echo "Excess diagnostic output while parsing ${fcl}."
    exit 1
  fi
  
done
