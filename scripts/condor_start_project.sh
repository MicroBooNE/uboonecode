#! /bin/bash
#------------------------------------------------------------------
#
# Purpose: A batch worker script for starting a sam project.
#
# Usage:
#
# condor_start_project.sh [options]
#
# --sam_user <arg>    - Specify sam user (default $GRID_USER).
# --sam_group <arg>   - Specify sam group (default "uboone").
# --sam_station <arg> - Specify sam station (default "uboone").
# --sam_defname <arg> - Sam dataset definition name (required).
# --sam_project <arg> - Sam project name (required).
# --outdir <arg>      - Specify output directory (optional). 
# -g, --grid          - Be grid-friendly.
#
# End options.
#
# Created: H. Greenlee, 29-Aug-2012
#
#------------------------------------------------------------------

# Parse arguments.

SAM_USER=$GRID_USER
SAM_GROUP="uboone"
SAM_STATION="uboone"
SAM_DEFNAME=""
SAM_PROJECT=""
OUTDIR=""
GRID=0
IFDH_OPT=""

while [ $# -gt 0 ]; do
  case "$1" in

    # Help.
    -h|--help )
      awk '/^# Usage:/,/^# End options/{print $0}' $0 | cut -c3- | head -n -2
      exit
      ;;

    # Sam user.
    --sam_user )
      if [ $# -gt 1 ]; then
        SAM_USER=$2
        shift
      fi
      ;;

    # Sam group.
    --sam_group )
      if [ $# -gt 1 ]; then
        SAM_GROUP=$2
        shift
      fi
      ;;

    # Sam station.
    --sam_station )
      if [ $# -gt 1 ]; then
        SAM_STATION=$2
        shift
      fi
      ;;

    # Sam dataset definition name.
    --sam_defname )
      if [ $# -gt 1 ]; then
        SAM_DEFNAME=$2
        shift
      fi
      ;;

    # Sam project name.
    --sam_project )
      if [ $# -gt 1 ]; then
        SAM_PROJECT=$2
        shift
      fi
      ;;

    # Output directory.
    --outdir )
      if [ $# -gt 1 ]; then
        OUTDIR=$2
        shift
      fi
      ;;

    # Grid flag.
    -g|--grid )
      GRID=1
      ;;

    # Other.
    * )
      echo "Unknown option $1"
      exit 1
  esac
  shift
done

# Done with arguments.

echo "Nodename: `hostname`"

# Check and print configuraiton options.

echo "Sam user: $SAM_USER"
echo "Sam group: $SAM_GROUP"
echo "Sam station: $SAM_STATION"
echo "Sam dataset definition: $SAM_DEFNAME"
echo "Sam project name: $SAM_PROJECT"

# Complain if SAM_DEFNAME is not defined.

if [ x$SAM_DEFNAME = x ]; then
  echo "Sam dataset was not specified (use option --sam_defname)."
  exit 1
fi

# Complain if SAM_PROJECT is not defined.

if [ x$SAM_PROJECT = x ]; then
  echo "Sam project name was not specified (use option --sam_project)."
  exit 1
fi

# Set options for ifdh.

if [ $GRID -ne 0 ]; then
  echo "X509_USER_PROXY = $X509_USER_PROXY"
  if ! echo $X509_USER_PROXY | grep -q Production; then
    IFDH_OPT="--force=expgridftp"
  fi
fi
echo "IFDH_OPT=$IFDH_OPT"

# Create the scratch directory in the condor scratch diretory.
# Copied from condor_lBdetMC.sh.
# Scratch directory path is stored in $TMP.
# Scratch directory is automatically deleted when shell exits.

# Do not change this section.
# It creates a temporary working directory that automatically cleans up all
# leftover files at the end.
TMP=`mktemp -d ${_CONDOR_SCRATCH_DIR}/working_dir.XXXXXXXXXX`
TMP=${TMP:-${_CONDOR_SCRATCH_DIR}/working_dir.$$}

{ [[ -n "$TMP" ]] && mkdir -p "$TMP"; } || \
  { echo "ERROR: unable to create temporary directory!" 1>&2; exit 1; }
trap "[[ -n \"$TMP\" ]] && { cd ; rm -rf \"$TMP\"; }" 0
cd $TMP
# End of the section you should not change.

echo "Scratch directory: $TMP"

# See if we need to set umask for group write.

if [ $GRID -eq 0 -a x$OUTDIR != x ]; then
  OUTUSER=`stat -c %U $OUTDIR`
  CURUSER=`whoami`
  if [ $OUTUSER != $CURUSER ]; then
    echo "Setting umask for group write."
    umask 002
  fi
fi

# Save the project name in a file.

echo $SAM_PROJECT > sam_project.txt

# Start the project.

echo "Starting project."
ifdh startProject $SAM_PROJECT $SAM_STATION $SAM_DEFNAME $SAM_USER $SAM_GROUP
if [ $? -eq 0 ]; then
  echo "Project successfully started."
else
  echo "Start project error status $?"
fi

# Stash all of the files we want to save in a local
# directory with a unique name.  Then copy this directory
# and its contents recursively.

if [ x$OUTDIR != x ]; then
  OUTPUT_SUBDIR=${CLUSTER}_start
  mkdir $OUTPUT_SUBDIR
  for outfile in *; do
    if [ $outfile != $OUTPUT_SUBDIR ]; then
      mv $outfile $OUTPUT_SUBDIR
    fi
  done
  echo "ifdh cp -r $IFDH_OPT $OUTPUT_SUBDIR ${OUTDIR}/$OUTPUT_SUBDIR"
  ifdh cp -r $IFDH_OPT $OUTPUT_SUBDIR ${OUTDIR}/$OUTPUT_SUBDIR
  stat=$?
  if [ $stat -ne 0 ]; then
    echo "ifdh cp failed with status ${stat}."
    exit $stat
  fi 
fi
