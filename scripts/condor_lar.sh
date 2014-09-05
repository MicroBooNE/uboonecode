#! /bin/bash
#------------------------------------------------------------------
#
# Purpose: A general purpose larsoft batch worker script.
#
# Adapted from condor_lBdetMC.sh by E. Church.
#
# Usage:
#
# condor_lar.sh [options]
#
# Lar options:
#
# -c, --config <arg>      - Configuration (fcl) file (required).
# -s, --source <arg>      - Input file (full path).
# -S, --source-list <arg> - Input file list (full path, one per line).
# -o, --output <arg>      - Output file name.
# -n, --nevts <arg>       - Number of events to process.
# --nskip <arg>           - Number of events to skip.
# --nfile <arg>           - Number of files to process per worker.
# --nfile_skip <arg>      - Number of files to skip (use with option -S).
# --args <args...>        - Arguments for lar command line (place at end).
#
# Sam and parallel project options.
#
# --sam_user <arg>        - Specify sam user (default $GRID_USER).
# --sam_group <arg>       - Specify sam group (default "uboone").
# --sam_station <arg>     - Specify sam station (default "uboone").
# --sam_defname <arg>     - Sam dataset definition name.
# --sam_project <arg>     - Sam project name.
# --njobs <arg>           - Parallel project with specified number of jobs (default one).
#
# Larsoft options.
#
# -r, --release <arg>     - Release tag.
# -q, -b, --build <arg>   - Release build qualifier (default "debug", or "prof").
# --localdir <arg>        - Larsoft local test release directory (default none).
# --localtar <arg>        - Tarball of local test release.
# --mrb                   - Ignored (for compatibility).
# --srt                   - Exit with error status (SRT run time no longer supported).
#
# Microboone options.
#
# --ubfcl <arg>           - Ubfcl version (default none, depracated).
#
# Other options.
#
# -h, --help              - Print help.
# -i, --interactive       - For interactive use.
# -g, --grid              - Be grid-friendly.
# --group <arg>           - Group or experiment (required).
# --workdir <arg>         - Work directory (required).
# --outdir <arg>          - Output directory (required).
# --scratch <arg>         - Scratch directory (only for interactive).
# --cluster <arg>         - Job cluster (override $CLUSTER)
# --process <arg>         - Process within cluster (override $PROCESS).
# --procmap <arg>         - Name of process map file (override $PROCESS).
# --init-script <arg>     - User initialization script execute.
# --init-source <arg>     - User initialization script to source (bash).
# --end-script <arg>      - User end-of-job script to execute.
#
# End options.
#
# Run time environment setup.
#
# MRB run-time environmental setup is controlled by four options:
#  --release (-r), --build (-b, -q), --localdir, and --localtar.  
#
# a) Use option --release or -r to specify uboonecode version.  Note that
#    larsoft is setup as a dependent product of uboonecode.
# b) Use option --build or -b to specify build full qualifiers (e.g. 
#    "debug:e5" or "e5:prof").
# c) Options --localdir or --localtar are used to specify your local
#    test release.  Use one or the other (not both).
#
#    Use --localdir to specify the location of your local install
#    directory ($MRB_INSTALL).
#
#    Use --localtar to specify thye location of a tarball of your
#    install directory (made relative to $MRB_INSTALL).
#
#    Note that --localdir is not grid-friendly.
#
# Notes.
#
# 1.  Each batch worker is uniquely identified by two numbers stored
#     in environment variables $CLUSTER and $PROCESS (the latter is 
#     a small integer that starts from zero and varies for different
#     jobs in a parallel job group).  These environment variables are
#     normally set by the batch system, but can be overridden by options 
#     --cluster, --process, and --procmap (e.g. to rerun failed jobs).
#
# 2.  The work directory must be set to an existing directory owned
#     by the submitter and readable by the batch worker.  Files from the 
#     work directory are copied to the batch worker scratch directory at
#     the start of the job.
#
# 3.  The job configuration file (-c option), initialization and end-of-job
#     scripts (optins --init-script, --init-source, --end-script) may
#     be stored in the work directory specified by option --workdir, or they
#     may be specified as absolute paths visible on the worker node.
#
# 4.  A local test release may be specified as an absolute path using
#     --localdir, or a tarball using --localtar.  The location of the tarball
#     may be specified as an absolute path visible on the worker, or a 
#     relative path relative to the work directory.
#
# 5.  The output directory must exist and be writable by the batch
#     worker (i.e. be group-writable for grid jobs).  The worker
#     makes a new subdirectory called ${CLUSTER}_${PROCESS} in the output
#     directory and copies all files in the batch scratch directory there 
#     at the end of the job.  If the output directory is not specified, the
#     default is /grid/data/<group>/outstage/<user> (user is defined as 
#     owner of work directory).
#
# 6.  Parallel projects are specified whenever --njobs is specified to
#     be greater than one.  Parallel projects are supported for single file,
#     file list, and sam project input.
#
#     In all cases, each worker processes some number of complete files.
#     If the number of jobs is greater than the number of input files, some
#     workers will not have any input files to process.
#
#     In any case, options --nfile and --nevts can be used to limit the
#     number of files or events that are processed by a single worker, 
#     regardless of the way files are divided among the workers.
#
#     Option --njobs is incompatible with options --nskip, and --nfile_skip.
#
#     a) Non-sam (single file or file list) input.
#
#     In this case, input files are preassigned to workers such that all input
#     files are approximately evenly divided among the workers.  All files
#     preassigned to this worker are copied to the scratch directory at the 
#     start of the job.
#
#     b) Sam project input.
#
#     In this case, files are assigned to workers in a non-deterministic
#     manner by the sam system.  The sam system fetches input files to the
#     scratch directory and deletes processed input files during job execution.
#
#
# 7.  Using option -n or --nevts to limit number of events processed:
#
#     a) If no input files are specified (e.g. mc generation), --nevts
#        specifies total number of events among all workers.
#
#     b) If input files are specified, --nevts specifies total number of 
#        events processed by each worker or from each input file, whichever
#        is less.
#
# 8.  The interactive option (-i or --interactive) allows this script
#     to be run interactively by overriding some settings that are normally
#     obtained from the batch system, including $CLUSTER, $PROCESS, and
#     the scratch directory.  Interactive jobs always set PROCESS=0 (unless
#     overridden by --process).
#
# 9.  The grid option (-g or --grid) instructs this script to use grid-
#     friendly tools.  This means that there must be no direct access to
#     bluearc disks.  File transfers are done using gridftp or other
#     grid-friendly protocol.  Local test releases are not allowed to 
#     be specified as directories (--localdir), but may be specified as
#     tarballs (--localtar).
#
#
# Created: H. Greenlee, 29-Aug-2012
#
#------------------------------------------------------------------

cd

# Parse arguments.

FCL=""
INFILE=""
INLIST=""
OUTFILE=""
NEVT=0
NSKIP=0
FIRST_EVENT=0
SUBRUN=1
NFILE=0
NFILE_SKIP=0
NJOBS=1
ARGS=""
REL=""
QUAL=""
LOCALDIR=""
LOCALTAR=""
UBFCL=""
INTERACTIVE=0
GRP=""
WORKDIR=""
OUTDIR=""
SCRATCH=""
CLUS=""
PROC=""
PROCMAP=""
INITSCRIPT=""
INITSOURCE=""
ENDSCRIPT=""
SAM_USER=$GRID_USER
SAM_GROUP="uboone"
SAM_STATION="uboone"
SAM_DEFNAME=""
SAM_PROJECT=""
USE_SAM=0
GRID=0
IFDH_OPT=""

while [ $# -gt 0 ]; do
  case "$1" in

    # Help.
    -h|--help )
      awk '/^# Usage:/,/^# End options/{print $0}' $0 | cut -c3- | head -n -2
      exit
      ;;

    # Config file.
    -c|--config )
      if [ $# -gt 1 ]; then
        FCL=$2
        shift
      fi
      ;;

    # Input file.
    -s|--source )
      if [ $# -gt 1 ]; then
        INFILE=$2
        shift
      fi
      ;;

    # Input file list.
    -S|--source-list )
      if [ $# -gt 1 ]; then
        INLIST=$2
        shift
      fi
      ;;

    # Output file.
    -o|--output )
      if [ $# -gt 1 ]; then
        OUTFILE=$2
        shift
      fi
      ;;

    # Number of events.
    -n|--nevts )
      if [ $# -gt 1 ]; then
        NEVT=$2
        shift
      fi
      ;;

    # Number of events to skip.
    --nskip )
      if [ $# -gt 1 ]; then
        NSKIP=$2
        shift
      fi
      ;;

    # Number of files to process.
    --nfile )
      if [ $# -gt 1 ]; then
        NFILE=$2
        shift
      fi
      ;;

    # Number of files to skip.
    --nfile_skip )
      if [ $# -gt 1 ]; then
        NFILE_SKIP=$2
        shift
      fi
      ;;

    # Number of parallel jobs.
    --njobs )
      if [ $# -gt 1 ]; then
        NJOBS=$2
        shift
      fi
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
        USE_SAM=1
        shift
      fi
      ;;

    # Sam project name.
    --sam_project )
      if [ $# -gt 1 ]; then
        SAM_PROJECT=$2
        USE_SAM=1
        shift
      fi
      ;;

    # General arguments for lar command line.
    --args )
      if [ $# -gt 1 ]; then
        shift
        ARGS=$@
        break
      fi
      ;;

    # Release tag.
    -r|--release )
      if [ $# -gt 1 ]; then
        REL=$2
        shift
      fi
      ;;

    # Release build qualifier.
    -q|-b|--build )
      if [ $# -gt 1 ]; then
        QUAL=$2
        shift
      fi
      ;;

    # Local test release directory.
    --localdir )
      if [ $# -gt 1 ]; then
        LOCALDIR=$2
        shift
      fi
      ;;

    # Local test release tarball.
    --localtar )
      if [ $# -gt 1 ]; then
        LOCALTAR=$2
        shift
      fi
      ;;

    # MRB flag.
    --mrb )
      ;;

    # SRT flag.
    --srt )
      echo "SRT run time environment is no longer supported."
      exit 1
      ;;

    # Ubfcl version.
    --ubfcl )
      if [ $# -gt 1 ]; then
        UBFCL=$2
        shift
      fi
      ;;

    # Interactive flag.
    -i|--interactive )
      INTERACTIVE=1
      ;;

    # Grid flag.
    -g|--grid )
      GRID=1
      ;;

    # Group.
    --group )
      if [ $# -gt 1 ]; then
        GRP=$2
        shift
      fi
      ;;

    # Work directory.
    --workdir )
      if [ $# -gt 1 ]; then
        WORKDIR=$2
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

    # Scratch directory.
    --scratch )
      if [ $# -gt 1 ]; then
        SCRATCH=$2
        shift
      fi
      ;;

    # Job cluster.
    --cluster )
      if [ $# -gt 1 ]; then
        CLUS=$2
        shift
      fi
      ;;

    # Process within cluster.
    --process )
      if [ $# -gt 1 ]; then
        PROC=$2
        shift
      fi
      ;;

    # Process map.
    --procmap )
      if [ $# -gt 1 ]; then
        PROCMAP=$2
        shift
      fi
      ;;

    # User initialization script.
    --init-script )
      if [ $# -gt 1 ]; then
        INITSCRIPT=$2
        shift
      fi
      ;;

    # User source initialization script.
    --init-source )
      if [ $# -gt 1 ]; then
        INITSOURCE=$2
        shift
      fi
      ;;

    # User end-of-job script.
    --end-script )
      if [ $# -gt 1 ]; then
        ENDSCRIPT=$2
        shift
      fi
      ;;

    # Other.
    * )
      echo "Unknown option $1"
      exit 1
  esac
  shift
done

#echo "FCL=$FCL"
#echo "INFILE=$INFILE"
#echo "INLIST=$INLIST"
#echo "OUTFILE=$OUTFILE"
#echo "NEVT=$NEVT"
#echo "NSKIP=$NSKIP"
#echo "NFILE=$NFILE"
#echo "NFILE_SKIP=$NFILE_SKIP"
#echo "NJOBS=$NJOBS"
#echo "ARGS=$ARGS"
#echo "REL=$REL"
#echo "QUAL=$QUAL"
#echo "LOCALDIR=$LOCALDIR"
#echo "LOCALTAR=$LOCALTAR"
#echo "INTERACTIVE=$INTERACTIVE"
#echo "GRP=$GRP"
#echo "WORKDIR=$WORKDIR"
#echo "OUTDIR=$OUTDIR"
#echo "SCRATCH=$SCRATCH"
#echo "CLUS=$CLUS"
#echo "PROC=$PROC"
#echo "INITSCRIPT=$INITSCRIPT"
#echo "INITSOURCE=$INITSOURCE"
#echo "ENDSCRIPT=$ENDSCRIPT"

# Done with arguments.

echo "Nodename: `hostname`"
id

# Set defaults.

if [ x$QUAL = x ]; then
  QUAL="debug:e5"
fi

# Initialize microboone ups products and mrb.

OASIS_DIR="/cvmfs/oasis.opensciencegrid.org/microboone/products/"
FERMIAPP_DIR="/grid/fermiapp/products/uboone/"

echo "Initializing ups and mrb."
  
if [[ -d "${FERMIAPP_DIR}" ]]; then
  echo "Sourcing ${FERMIAPP_DIR}setup_uboone.sh file"
  source ${FERMIAPP_DIR}/setup_uboone.sh

elif [[ -d "${OASIS_DIR}" ]]; then
  echo "Sourcing the ${OASIS_DIR}setup_uboone.sh file"
  source ${OASIS_DIR}/setup_uboone.sh

else
  echo "Could not find MRB initialization script setup_uboone.sh"
  exit 1
fi

# Ifdh may already be setup by jobsub wrapper.
# If not, set it up here.

echo "IFDHC_DIR=$IFDHC_DIR"
if [ x$IFDHC_DIR = x ]; then
  echo "Setting up ifdhc, because jobsub did not set it up."
  setup ifdhc
fi
echo "IFDHC_DIR=$IFDHC_DIR"

# Set GROUP environment variable.

unset GROUP
if [ x$GRP != x ]; then
  GROUP=$GRP
else
  echo "GROUP not specified."
  exit 1  
fi
export GROUP
echo "Group: $GROUP"

# Set options for ifdh.

if [ $GRID -ne 0 ]; then

  # Figure out if this is a production job.
  # This option is only used when copying back output.
  # It affects the ownership of copied back files.

  echo "X509_USER_PROXY = $X509_USER_PROXY"
  if ! echo $X509_USER_PROXY | grep -q Production; then
    IFDH_OPT="--force=expgridftp"
  fi
fi
echo "IFDH_OPT=$IFDH_OPT"

# Make sure work directory is defined and exists.

if [ x$WORKDIR = x ]; then
  echo "Work directory not specified."
  exit 1
fi
if [ $GRID -eq 0 -a ! -d $WORKDIR ]; then
  echo "Work directory $WORKDIR does not exist."
  exit 1
fi
echo "Work directory: $WORKDIR"

# Make sure fcl file argument was specified.

if [ x$FCL = x ]; then
  echo "No configuration option (-c|--config) was specified."
  exit 1
fi

# Make sure output directory exists and is writable.

if [ x$OUTDIR = x ]; then
  echo "Output directory not specified."
  exit 1
fi
if [ $GRID -eq 0 -a \( ! -d $OUTDIR -o ! -w $OUTDIR \) ]; then
  echo "Output directory $OUTDIR does not exist or is not writable."
  exit 1
fi
echo "Output directory: $OUTDIR"

# See if we need to set umask for group write.

if [ $GRID -eq 0 ]; then
  OUTUSER=`stat -c %U $OUTDIR`
  CURUSER=`whoami`
  if [ $OUTUSER != $CURUSER ]; then
    echo "Setting umask for group write."
    umask 002
  fi
fi

# Make sure scratch directory is defined.
# For batch, the scratch directory is always $_CONDOR_SCRATCH_DIR
# For interactive, the scratch directory is specified by option 
# --scratch or --outdir.

if [ $INTERACTIVE -eq 0 ]; then
  SCRATCH=$_CONDOR_SCRATCH_DIR
else
  if [ x$SCRATCH = x ]; then
    SCRATCH=$OUTDIR
  fi
fi
if [ x$SCRATCH = x -o ! -d "$SCRATCH" -o ! -w "$SCRATCH" ]; then
  echo "Local scratch directory not defined or not writable."
  exit 1
fi

# Create the scratch directory in the condor scratch diretory.
# Copied from condor_lBdetMC.sh.
# Scratch directory path is stored in $TMP.
# Scratch directory is automatically deleted when shell exits.

# Do not change this section.
# It creates a temporary working directory that automatically cleans up all
# leftover files at the end.
TMP=`mktemp -d ${SCRATCH}/working_dir.XXXXXXXXXX`
TMP=${TMP:-${SCRATCH}/working_dir.$$}

{ [[ -n "$TMP" ]] && mkdir -p "$TMP"; } || \
  { echo "ERROR: unable to create temporary directory!" 1>&2; exit 1; }
trap "[[ -n \"$TMP\" ]] && { cd ; rm -rf \"$TMP\"; }" 0
cd $TMP
# End of the section you should not change.

echo "Scratch directory: $TMP"

# Copy files from work directory to scratch directory.

echo "Fetching files from work directory."
ifdh cp -r $WORKDIR work
stat=$?
if [ $stat -ne 0 ]; then
  echo "ifdh cp failed with status ${stat}."
  exit $stat
fi 
cd work
echo "Local working directoroy:"
pwd
ls
echo

# Save the hostname and condor job id.

hostname > hostname.txt
echo ${CLUSTER}.${PROCESS} > jobid.txt

# Set default CLUSTER and PROCESS environment variables for interactive jobs.

if [ $INTERACTIVE -ne 0 ]; then
  CLUSTER=`date +%s`   # From time stamp.
  PROCESS=0            # Default zero for interactive.
fi

# Override CLUSTER and PROCESS from command line options.

if [ x$CLUS != x ]; then
  CLUSTER=$CLUS
fi
if [ x$PROC != x ]; then
  PROCESS=$PROC
fi
if [ x$PROCMAP != x ]; then
  if [ -f $PROCMAP ]; then
    PROCESS=`sed -n $(( $PROCESS + 1 ))p $PROCMAP`
  else
    echo "Process map file $PROCMAP not found."
    exit 1
  fi
fi
if [ x$CLUSTER = x ]; then
  echo "CLUSTER not specified."
  exit 1
fi
if [ x$PROCESS = x ]; then
  echo "PROCESS not specified."
  exit 1
fi
echo "Procmap: $PROCMAP"
echo "Cluster: $CLUSTER"
echo "Process: $PROCESS"

# Construct name of output subdirectory.

OUTPUT_SUBDIR=${CLUSTER}_${PROCESS}
echo "Output subdirectory: $OUTPUT_SUBDIR"

# Copy stuff back to output directory before we actually start running lar.

#ifdh mkdir $IFDH_OPT ${OUTDIR}/${OUTPUT_SUBDIR}/hostname.txt
#ifdh cp $IFDH_OPT hostname.txt ${OUTDIR}/${OUTPUT_SUBDIR}/hostname.txt
#ifdh cp $IFDH_OPT jobid.txt ${OUTDIR}/${OUTPUT_SUBDIR}/jobid.txt

# Make sure fcl file exists.

if [ ! -f $FCL ]; then
  echo "Configuration file $FCL does not exist."
  exit 1
fi

# Make sure init script exists and is executable (if specified).

if [ x$INITSCRIPT != x -a ! -x "$INITSCRIPT" ]; then
  echo "Initialization script $INITSCRIPT does not exist."
  exit 1
fi

# Make sure init source script exists (if specified).

if [ x$INITSOURCE != x -a ! -f "$INITSOURCE" ]; then
  echo "Initialization script $INITSOURCE does not exist."
  exit 1
fi

# Make sure end-of-job script exists and is executable (if specified).

if [ x$ENDSCRIPT != x -a ! -x "$ENDSCRIPT" ]; then
  echo "Initialization script $ENDSCRIPT does not exist."
  exit 1
fi

# MRB run time environment setup goes here.

# Setup local test release, if any.

if [ x$LOCALDIR != x ]; then
  #if [ $GRID -ne 0 ]; then
  #  echo "Local test release may not be specified in grid-friendly mode."
  #  exit 1
  #fi
  echo "Setting up local test release in directory ${LOCALDIR}."
  if [ ! -d $LOCALDIR ]; then
    echo "Local test release directory $LOCALDIR does not exist."
    exit 1
  fi
  if [ ! -f ${LOCALDIR}/setup ]; then
    echo "Local test release directory $LOCALDIR does not contain a setup script."
    exit
  fi
  echo "Initializing localProducts from ${LOCALDIR}."
  . $LOCALDIR/setup
  #echo "MRB_INSTALL=${MRB_INSTALL}."
  #echo "MRB_QUALS=${MRB_QUALS}."
  echo "Setting up all localProducts."
  if [ x$IFDHC_DIR != x ]; then
    unsetup ifdhc
  fi
  mrbslp
fi
cd $TMP/work

# Setup local larsoft test release from tarball.

if [ x$LOCALTAR != x ]; then
  mkdir $TMP/local
  cd $TMP/local

  # Fetch the tarball.

  echo "Fetching test release tarball ${LOCALTAR}."

  # Make sure ifdhc is setup.

  if [ x$IFDHC_DIR = x ]; then
    echo "Setting up ifdhc before fetching tarball."
    setup ifdhc
  fi
  echo "IFDHC_DIR=$IFDHC_DIR"
  ifdh cp $LOCALTAR local.tar
  stat=$?
  if [ $stat -ne 0 ]; then
    echo "ifdh cp failed with status ${stat}."
    exit $stat
  fi 

  # Extract the tarball.

  tar -xf local.tar

  # Setup the environment.

  cd $TMP/work
  echo "Initializing localProducts from tarball ${LOCALTAR}."
  sed "s@setenv MRB_INSTALL.*@setenv MRB_INSTALL ${TMP}/local@" $TMP/local/setup | \
  sed "s@setenv MRB_TOP.*@setenv MRB_TOP ${TMP}@" > $TMP/local/setup.local
  . $TMP/local/setup.local
  #echo "MRB_INSTALL=${MRB_INSTALL}."
  #echo "MRB_QUALS=${MRB_QUALS}."
  echo "Setting up all localProducts."
  if [ x$IFDHC_DIR != x ]; then
    unsetup ifdhc
  fi
  mrbslp
fi

# Setup specified version of uboonecode (if specified, and if local
# test release did not set it up).

if [ x$UBOONECODE_DIR == x -a x$REL != x ]; then
  echo "Setting up uboonecode $REL -q ${QUAL}."
  if [ x$IFDHC_DIR != x ]; then
    unsetup ifdhc
  fi
  setup uboonecode $REL -q $QUAL
fi

cd $TMP/work

# In case mrb setup didn't setup a version of ifdhc, set up ifdhc again.

if [ x$IFDHC_DIR = x ]; then
  echo "Setting up ifdhc again, because larsoft did not set it up."
  setup ifdhc
fi
echo "IFDH_ART_DIR=$IFDH_ART_DIR"
echo "IFDHC_DIR=$IFDHC_DIR"

# Set up ubfcl product version, if any.

if [ x$UBFCL != x ]; then
  setup ubfcl $UBFCL
fi 
echo "UBFCL_DIR=${UBFCL_DIR}"

# Get input files to process, either single file, file list, or sam.
#
# For non-sam input, copy all files local using ifdh cp, and make a 
# local file list called input.list.  Save the remote file names (uri's)
# in another file called transferred_uris.list

rm -f input.list
rm -f transferred_uris.list
NFILE_TOTAL=0

if [ $USE_SAM -eq 0 -a x$INFILE != x ]; then

  # Single file case.

  # Don't allow any list-related options in single file case:
  # -S, --source-list, --nfile, --nfile_skip

  if [ x$INLIST != x -o $NFILE -ne 0 -o $NFILE_SKIP -ne 0 ]; then
    echo "File list options specified with single input file."
    exit 1
  fi

  # Copy input file to scratch directoroy.

  NFILE_TOTAL=1
  LOCAL_INFILE=`basename $INFILE`
  echo "Copying $INFILE"
  ifdh cp $INFILE $LOCAL_INFILE
  stat=$?
  if [ $stat -ne 0 ]; then
    echo "ifdh cp failed with status ${stat}."
    exit $stat
  fi 
  if [ -f $LOCAL_INFILE -a $stat -eq 0 ]; then
    echo $INFILE > transferred_uris.list
    echo $LOCAL_INFILE > input.list
  else
    echo "Error fetching input file ${INFILE}."
    exit 1
  fi

elif [ $USE_SAM -eq 0 -a x$INLIST != x ]; then

  # Input list case.

  # Make sure input file list exists.

  if [ ! -f $INLIST ]; then
    echo "Input file list $INLIST does not exist."
    exit 1
  fi

  # Remember how many files are in the input file list.

  NFILE_TOTAL=`cat $INLIST | wc -l`
  echo "Input file list contains $NFILE_TOTAL total files."

  # If --njobs was specified, calculate how many files
  # to skip and process in this worker.

  if [ $NJOBS -ne 0 ]; then

    # Don't allow option --nfile_skip in this case.

    if [ $NFILE_SKIP -ne 0 ]; then
      echo "Illegal options specified with --njobs."
      exit 1
    fi

    # Clamp NJOBS to be a maximum of $NFILE_TOTAL.
    # This means that workers with $PROCESS >= $NFILE_TOTAL will not have 
    # any input files to process.

    MYNJOBS=$NJOBS
    if [ $MYNJOBS -gt $NFILE_TOTAL ]; then
      MYNJOBS=$NFILE_TOTAL
    fi

    # Calculate number of files to skip and number of files to process.

    NFILE_SKIP=$(( $PROCESS * $NFILE_TOTAL / $MYNJOBS ))
    MYNFILE=$(( ( $PROCESS + 1 ) * $NFILE_TOTAL / $MYNJOBS - $NFILE_SKIP ))
    if [ $MYNFILE -eq 0 -o $NFILE_SKIP -ge $NFILE_TOTAL ]; then
      echo "This worker did not get any input files."
      exit 1
    fi
    if [ $MYNFILE -lt $NFILE -o $NFILE -eq 0 ]; then
      NFILE=$MYNFILE
    fi
  fi

  # Report number of files to skip and process.

  echo "Skipping $NFILE_SKIP files."
  if [ $NFILE -eq 0 ]; then
    echo "Processing all remaining files."
  else
    echo "Processing $NFILE files."
  fi

  # Copy input files and construct local input file list.

  nfile=0
  nfskip=$NFILE_SKIP
  nmax=$NFILE
  while read infile; do
    if [ $nfskip -gt 0 ]; then
      nfskip=$(( $nfskip - 1 ))
    else

      # Retain the original file name as the local file name, if possible.
      # Otherwise, generate a new (hopefully) unique name.

      if [ ! -f input.list ]; then
        touch input.list
      fi
      LOCAL_INFILE=`basename $infile`
      if grep -q $LOCAL_INFILE input.list; then
        LOCAL_INFILE=input${nfile}.root
      fi
      echo "Copying $infile"
      ifdh cp $infile $LOCAL_INFILE
      stat=$?
      if [ $stat -ne 0 ]; then
        echo "ifdh cp failed with status ${stat}."
        exit $stat
      fi 
      if [ -f $LOCAL_INFILE -a $stat -eq 0 ]; then
        echo $infile >> transferred_uris.list
        echo $LOCAL_INFILE >> input.list
      else
        echo "Error fetching input file ${infile}."
        exit 1
      fi
      nmax=$(( $nmax - 1 ))
      if [ $nmax -eq 0 ]; then
        break
      fi
    fi
    nfile=$(( $nfile + 1 ))
  done < $INLIST
fi

NFILE_LOCAL=0
if [ $USE_SAM -eq 0 ]; then
  if [ -f input.list ]; then
    echo "Local input file list:"
    cat input.list
    NFILE_LOCAL=`cat input.list | wc -l`
  else
    echo "No local input files."
  fi
  echo "Local input list has $NFILE_LOCAL files."
fi

# In case no input files were specified, and we are not getting input
# from sam (i.e. mc generation), recalculate the first event number,
# the subrun number, and the number of events to generate in this worker.

if [ $USE_SAM -eq 0 -a $NFILE_TOTAL -eq 0 ]; then

  # Don't allow --nskip.

  if [ $NSKIP -gt 0 ]; then
    echo "Illegal option --nskip specified with no input."
    exit 1
  fi

  # Do calculation.

  NSKIP=$(( $PROCESS * $NEVT / $NJOBS ))
  NEV=$(( ( $PROCESS + 1 ) * $NEVT / $NJOBS - $NSKIP ))
  FIRST_EVENT=$(( $NSKIP + 1 ))
  NSKIP=0
  NEVT=$NEV

  # Set subrun=$PROCESS+1 in a wrapper fcl file.

  SUBRUN=$(( $PROCESS + 1))
  cat <<EOF > subrun_wrapper.fcl
#include "$FCL"

source.firstSubRun: $SUBRUN

EOF
  FCL=subrun_wrapper.fcl
  
  echo "First MC event: $FIRST_EVENT"
  echo "MC subrun: $SUBRUN"
  echo "Number of MC events: $NEVT"
fi

# For sam input, start project (if necessary), and consumer process.

PURL=''
CPID=''
if [ $USE_SAM -ne 0 ]; then

  # Make sure a project name has been specified.

  if [ x$SAM_PROJECT = x ]; then
    echo "No sam project was specified."
    exit 1
  fi
  echo "Sam project: $SAM_PROJECT"

  # If a sam definition name was specified, attempt to start a project.
  # Issue the start project command, even though the project might
  # already be running.  In the case of multiple-worker projects, 
  # the first worker to issue this command will start the project.

  #if [ x$SAM_DEFNAME != x ]; then
  #  echo "Sam dataset: $SAM_DEFNAME"
  #  echo "Attempting to start project."
  #  ifdh startProject $SAM_PROJECT $SAM_STATION $SAM_DEFNAME $SAM_USER $SAM_GROUP
  #  echo "Start project status $?"
  #fi

  # Get the project url of a running project (maybe the one we just started,
  # or maybe started externally).  This command has to succeed, or we can't
  # continue.

  PURL=`ifdh findProject $SAM_PROJECT $SAM_STATION`
  if [ x$PURL = x ]; then
    echo "Unable to find url for project ${SAM_PROJECT}."
    exit 1
  else
    echo "Project url: $PURL"
  fi

  # Start the consumer process.  This command also has to succeed.

  NODE=`hostname`
  APPFAMILY=art

  # Use lar to parse fcl file to extract process_name, and use that
  # as the application name for starting the consumer process.

  export ART_DEBUG_CONFIG=1
  APPNAME=`lar -c $FCL 2>&1 > /dev/null | grep process_name: | awk '{print $2}'`
  unset ART_DEBUG_CONFIG
  if [ x$APPNAME = x ]; then
    echo "Trouble determining application name."
    exit 1
  fi

  CPID=`ifdh establishProcess $PURL $APPNAME $REL $NODE $SAM_USER $APPFAMILY $FCL $NFILE`
  if [ x$CPID = x ]; then
    echo "Unable to start consumer process for project url ${PURL}."
    exit 1
  else
    echo "Consumer process id $CPID"
  fi

  # Stash away the project name and consumer process id in case we need them
  # later for bookkeeping.

  echo $SAM_PROJECT > sam_project.txt
  echo $CPID > cpid.txt

  # Generate an fcl wrapper for all sam-related fcl parameters.

  cat <<EOF > sam_wrapper.fcl
#include "$FCL"

services.user.IFDH:
{
  IFDH_BASE_URI: "http://samweb.fnal.gov:8480/sam/uboone/api"
}

services.user.CatalogInterface:
{
  service_provider: "IFCatalogInterface"
  webURI: "$PURL"
}

services.user.FileTransfer:
{
  service_provider: "IFFileTransfer"
}

source.fileNames: [ "$CPID" ]

EOF
  FCL=sam_wrapper.fcl
fi

# Construct options for lar command line.

LAROPT="-c $FCL"
if [ -f input.list ]; then
  LAROPT="$LAROPT -S input.list"
fi
if [ x$OUTFILE != x ]; then
  LAROPT="$LAROPT -o $OUTFILE"
fi
if [ $NEVT -ne 0 ]; then
  LAROPT="$LAROPT -n $NEVT"  
fi
if [ $NSKIP -ne 0 ]; then
  LAROPT="$LAROPT --nskip $NSKIP"
fi
if [ $FIRST_EVENT -ne 0 ]; then
  LAROPT="$LAROPT -e $FIRST_EVENT"
fi
if [ -n "$ARGS" ]; then
  LAROPT="$LAROPT $ARGS"  
fi

# Run/source optional initialization scripts.

if [ x$INITSCRIPT != x ]; then
  echo "Running initialization script ${INITSCRIPT}."
  if ! ./${INITSCRIPT}; then
    exit $?
  fi
fi
if [ x$INITSOURCE != x ]; then
  echo "Sourceing initialization script ${INITSOURCE}."
  . $INITSOURCE
  status=$?
  if [ $status -ne 0 ]; then
    exit $status
  fi
fi

# Save a copy of the environment, which can be helpful for debugging.

env > env.txt

# Save a canonicalized version of the fcl configuration.

ART_DEBUG_CONFIG=cfg.fcl lar -c $FCL

# Run lar.

echo "lar $LAROPT"
echo "lar $LAROPT" > command.txt
lar $LAROPT > lar.out 2> lar.err
stat=$?
echo $stat > lar.stat
echo "lar completed with exit status ${stat}."

# Setup up current version of ifdhc (may be different than version setup by larsoft).

#echo "Setting up current version of ifdhc."
#if [ x$IFDHC_DIR != x ]; then
#  unsetup ifdhc
#fi
#setup ifdhc v1_3_2
echo "IFDHC_DIR=$IFDHC_DIR"

# Sam cleanups.

if [ $USE_SAM -ne 0 ]; then

  # Get list of consumed files.

  ifdh translateConstraints "consumer_process_id $CPID and consumed_status consumed" > consumed_files.list

  # End consumer process.

  ifdh endProcess $PURL $CPID
fi

# Delete input files.

if [ $USE_SAM -eq 0 -a -f input.list ]; then
  while read file; do
    rm -f $file
  done < input.list
fi

# Run optional end-of-job script.

if [ x$ENDSCRIPT != x ]; then
  echo "Running end-of-job script ${ENDSCRIPT}."
  if ! ./${ENDSCRIPT}; then
    exit $?
  fi
fi

# Do root file checks.

for root in *.root; do
  root_metadata.py $root > ${root}.json
done

# Make local files group write, if appropriate.

if [ $GRID -eq 0 -a $OUTUSER != $CURUSER ]; then
  chmod -R g+rw .
fi

# Stash all of the files we want to save in a local
# directory with a unique name.  Then copy this directory
# and its contents recursively.

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

