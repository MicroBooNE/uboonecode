#! /bin/bash
#------------------------------------------------------------------
#
# Purpose: A general purpose larsoft batch submit script.
#
# Usage:
#
# submit_lar.sh [options]
#
# Lar options:
#
# -c, --config <arg>      - Configuration (fcl) file (required).
# -s, --source <arg>      - Input file (full path).
# -S, --source-list <arg> - Input file list (full path, one per line).
# -o, --output <arg>      - Output file name.
# -n, --nevts <arg>       - Number of events to process.
# --nskip <arg>           - Number of events to skip.
# --nfile <arg>           - Number of files to process (use with option -S).
# --nfile_skip <arg>      - Number of files to skip (use with option -S).
# --args <args...>        - Arguments for lar command line (place at end).
#
# Larsoft options.
#
# -r, --release <arg>     - Release tag (default "development").          
# -q, -b, --build <arg>   - Release build qualifier (default "debug", or "prof").
# --localdir <arg>        - Larsoft local test release directory (default none).
# --localtar <arg>        - Tarball of local test release.
# --mrb                   - Use mrb-style environment initialization (see below).
# --srt                   - Use srt-style environment initialization (see below).
#
# Jobsub options.
#
# -N, --numjobs           - Number of jobs to submit (default 1).
# --opportunistic         - Use opportunistic grid.
# --OS <arg>              - Specify OS (comma-separated list with no spaces, e.g. SL5,SL6).
#
# Microboone options.
#
# --ubfcl <arg>           - Ubfcl version (default none).
#
# Other options.
#
# -h, --help              - Print help.
# --group <arg>           - Group or experiment (default $GROUP).
# --workdir <arg>         - Work directory (default current directory).
# --outdir <arg>          - Output directory.
# --cluster <arg>         - Job cluster.
# --process <arg>         - Process within cluster.
# --procmap <arg>         - Name of process map file (override $PROCESS).
# --init-script <arg>     - User initialization script execute.
# --init-source <arg>     - User initialization script to source (bash).
# --end-script <arg>      - User end-of-job script to execute.
#
# End options.
#
# MRB vs. SRT environment setup.
#
# Environmental setup inside the worker batch script is controlled by 
# options --release (-r), --build (-b, -q), --localdir, and --localtar.  These options are 
# interpreted differently depending on whether one is using an MRB (--mrb) 
# or SRT (--srt) environment.  Options --srt and --mrb are mutually
# exclusive.  If neither is specified, --mrb is the default.
#
# For SRT-style environment:
#
# a) Use option --release or -r to specify larsoft frozen release.
# b) Use option --build, -b, or -q to specify build type ("debug" or "prof")
# c) If you have a local test release, specify the (bluearc) location of
#    the test release using --localdir, or use --localtar to specify the 
#    absolute or relative path of the local release tarball.
#
# For MRB-style environment:
#
# a) Use option --release or -r to specify uboonecode version.  Note that
#    larsoft is setup as a dependent product of uboonecode.
# b) Use option --build or -b to specify build full qualifiers (e.g. 
#    "debug:e4" or "e4:prof").
# c) Options --localdir and --localtar are used in a similar way as for
#    SRT.  Here, --localdir should point to, or your tarball should be 
#    made relative to, the mrb localProducts directory ($MRB_INSTALL).
#    If your local test release includes uboonecode, then you do not
#    need options --release or --build.
#
# Notes.
#
# 1.  Each batch worker is uniquely identified by two numbers stored
#     in environment variables $CLUSTER and $PROCESS (the latter is 
#     a small integer that starts from zero and varies for different
#     jobs in a parallel job group).  These environment variables are
#     normally set by the batch system, but can be overridden by options 
#     --cluster and --process (e.g. to rerun failed jobs).
#
#     It is not allowed to specify multiple jobs (-N option with >1 job) 
#     and to specify --process at the same time.
#
# 2.  The work directory must be set to an existing directory owned
#     by the submitter and readable by the batch worker.  Files from the 
#     work directory are copied to the batch worker scratch directory at
#     the start of the job.
#
# 3.  The job configuration file (-c option), initialization and end-of-job
#     scripts (optins --init-script, --init-source, --end-script) may
#     be stored in the work directory specified by option --workdir, or they
#     may be specified as absolute paths.
#
# 4.  The output directory must exist and be writable by the batch
#     worker (i.e. be group-writable for grid jobs).  The worker
#     makes a new subdirectory called ${CLUSTER}_${PROCESS} in the output
#     directory and copies all files in the batch scratch directory there 
#     at the end of the job.  If the output directory is not specified, the
#     default is /grid/data/<group>/outstage/<user> (user is defined as 
#     owner of work directory).
#
# 5.  The number of jobs option (-N or --numjobs) is used in two ways.
#     It is passed to jobsub to generate the requested number of batch jobs,
#     and it is passed to condor_lar.sh (as option --njobs) to control how
#     condor_lar.sh delivers files and events to individual workers.
#
# Created: H. Greenlee, 29-Aug-2012
#
#------------------------------------------------------------------

# Parse arguments.

FCL=""
INFILE=""
INLIST=""
OUTFILE=""
NEVT=0
NSKIP=0
NFILE=0
NFILE_SKIP=0
ARGS=""
REL=""
QUAL=""
LOCALDIR=""
LOCALTAR=""
MRB=0
SRT=0
UBFCL=""
NUMJOBS=1
OPPORTUNISTIC=0
GRP=$GROUP
OS=""
WORKDIR=`pwd`
OUTDIR=""
CLUS=""
PROC=""
PROCMAP=""
INITSCRIPT=""
INITSOURCE=""
ENDSCRIPT=""

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
      MRB=1
      ;;

    # SRT flag.
    --srt )
      SRT=1
      ;;

    # Ubfcl version.
    --ubfcl )
      if [ $# -gt 1 ]; then
        UBFCL=$2
        shift
      fi
      ;;

    # Number of batch jobs.
    -N|--numjobs )
      if [ $# -gt 1 ]; then
        NUMJOBS=$2
        shift
      fi
      ;;

    # Opportunistic grid flag.
    --opportunistic )
      OPPORTUNISTIC=1
      ;;

    # Group.
    -g|--group )
      if [ $# -gt 1 ]; then
        GRP=$2
        shift
      fi
      ;;

    # OS.
    --OS )
      if [ $# -gt 1 ]; then
        OS=$2
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

# Do MRB/SRT checks, and set some environment-specific defaults.

# Make sure at lease one of MRB and SRT is specified.
# If neither is specified, set MRB by default.

if [ $MRB -eq 0 -a $SRT -eq 0 ]; then
  MRB=1
fi

# Make sure both MRB and SRT are not set.

if [ $MRB -ne 0 -a $SRT -ne 0 ]; then
  echo "Both --mrb and --srt were specified."
  exit 1
fi

# Set defaults for MRB.

if [ $MRB -ne 0 ]; then
  if [ x$QUAL = x ]; then
    QUAL="debug:e4"
  fi
fi

# Set defaults for SRT.

if [ $SRT -ne 0 ]; then
  if [ x$REL = x ]; then
    REL="development"
  fi
  if [ x$QUAL = x ]; then
    QUAL="debug"
  fi
fi

#echo "FCL=$FCL"
#echo "INFILE=$INFILE"
#echo "INLIST=$INLIST"
#echo "OUTFILE=$OUTFILE"
#echo "NEVT=$NEVT"
#echo "NSKIP=$NSKIP"
#echo "NFILE=$NFILE"
#echo "NFILE_SKIP=$NFILE_SKIP"
#echo "ARGS=$ARGS"
#echo "REL=$REL"
#echo "QUAL=$QUAL"
#echo "LOCALDIR=$LOCALDIR"
#echo "LOCALTAR=$LOCALTAR"
#echo "MRB=$MRB"
#echo "SRT=$SRT"
#echo "UBFCL=$UBFCL"
#echo "NUMJOBS=$NUMJOBS"
#echo "OPPORTUNISTIC=$OPPORTUNISTIC"
#echo "GRP=$GRP"
#echo "WORKDIR=$WORKDIR"
#echo "OUTDIR=$OUTDIR"
#echo "CLUS=$CLUS"
#echo "PROC=$PROC"
#echo "PROCMAP=$PROCMAP"
#echo "INITSCRIPT=$INITSCRIPT"
#echo "INITSOURCE=$INITSOURCE"
#echo "ENDSCRIPT=$ENDSCRIPT"

# Done with arguments.

# Make sure all required options have been specified.

if [ x$GROUP = x ]; then
  echo "Group not specified."
  exit 1
fi
if [ x$OUTDIR = x ]; then
  echo "Output directory not specified."
  exit 1
fi
if [ x$WORKDIR = x ]; then
  echo "Work directory not specified."
  exit 1
fi
if [ x$FCL = x ]; then
  echo "Configuration file not specified."
  exit 1
fi

# Construct options for jobsub

JOBSUB_OPTS=''
if [ x$GRP != x ]; then
  JOBSUB_OPTS="$JOBSUB_OPTS --group=$GRP"
fi
if [ $NUMJOBS -gt 1 ]; then
  JOBSUB_OPTS="$JOBSUB_OPTS -N $NUMJOBS"
fi
JOBSUB_OPTS="$JOBSUB_OPTS --grid"
if [ $OPPORTUNISTIC -ne 0 ]; then
  JOBSUB_OPTS="$JOBSUB_OPTS --opportunistic"
fi
if [ x$OS != x ]; then
  JOBSUB_OPTS="$JOBSUB_OPTS --OS=$OS"
fi

# Construct condor_lar.sh options.

BATCH_OPTS=''
if [ x$FCL != x ]; then
  BATCH_OPTS="$BATCH_OPTS -c $FCL"
fi
if [ x$INFILE != x ]; then
  BATCH_OPTS="$BATCH_OPTS -s $INFILE"
fi
if [ x$INLIST != x ]; then
  BATCH_OPTS="$BATCH_OPTS -S $INLIST"
fi
if [ x$OUTFILE != x ]; then
  BATCH_OPTS="$BATCH_OPTS -o $OUTFILE"
fi
if [ $NEVT -gt 0 ]; then
  BATCH_OPTS="$BATCH_OPTS -n $NEVT"
fi
if [ $NSKIP -gt 0 ]; then
  BATCH_OPTS="$BATCH_OPTS --nskip $NSKIP"
fi
if [ $NFILE -gt 0 ]; then
  BATCH_OPTS="$BATCH_OPTS --nfile $NFILE"
fi
if [ $NFILE_SKIP -gt 0 ]; then
  BATCH_OPTS="$BATCH_OPTS --nfile_skip $NFILE_SKIP"
fi
if [ $NUMJOBS -gt 1 ]; then
  BATCH_OPTS="$BATCH_OPTS --njobs $NUMJOBS"
fi
if [ x$REL != x ]; then
  BATCH_OPTS="$BATCH_OPTS -r $REL"
fi
if [ x$QUAL != x ]; then
  BATCH_OPTS="$BATCH_OPTS -q $QUAL"
fi
if [ x$LOCALDIR != x ]; then
  BATCH_OPTS="$BATCH_OPTS --localdir $LOCALDIR"
fi
if [ x$LOCALTAR != x ]; then
  BATCH_OPTS="$BATCH_OPTS --localtar $LOCALTAR"
fi
if [ $MRB -ne 0 ]; then
  BATCH_OPTS="$BATCH_OPTS --mrb"
fi
if [ $SRT -ne 0 ]; then
  BATCH_OPTS="$BATCH_OPTS --srt"
fi
if [ x$UBFCL != x ]; then
  BATCH_OPTS="$BATCH_OPTS --ubfcl $UBFCL"
fi
if [ x$GRP != x ]; then
  BATCH_OPTS="$BATCH_OPTS --group $GRP"
fi
if [ x$WORKDIR != x ]; then
  BATCH_OPTS="$BATCH_OPTS --workdir $WORKDIR"
fi
if [ x$OUTDIR != x ]; then
  BATCH_OPTS="$BATCH_OPTS --outdir $OUTDIR"
fi
if [ x$CLUS != x ]; then
  BATCH_OPTS="$BATCH_OPTS --cluster $CLUS"
fi
if [ x$PROC != x ]; then
  BATCH_OPTS="$BATCH_OPTS --process $PROC"
fi
if [ x$PROCMAP != x ]; then
  BATCH_OPTS="$BATCH_OPTS --procmap $PROCMAP"
fi
if [ x$INITSCRIPT != x ]; then
  BATCH_OPTS="$BATCH_OPTS --init-script $INITSCRIPT"
fi
if [ x$INITSOURCE != x ]; then
  BATCH_OPTS="$BATCH_OPTS --init-source $INITSOURCE"
fi
if [ x$ENDSCRIPT != x ]; then
  BATCH_OPTS="$BATCH_OPTS --end-script $ENDSCRIPT"
fi

if [ -n "$ARGS" ]; then
  BATCH_OPTS="$BATCH_OPTS $ARGS"
fi

# Make a fresh copy of batch worker script in work directory and go there.

cp `which condor_lar.sh` $WORKDIR
cd $WORKDIR

# Submit jobs.

CMD="jobsub $JOBSUB_OPTS condor_lar.sh $BATCH_OPTS"
echo $CMD
eval $CMD
