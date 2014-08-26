#! /bin/bash
#----------------------------------------------------------------------
#
# Name: make_xml_mcc5.0.sh
#
# Purpose: Make xml files for mcc 5.0.  This script loops over all
#          generator-level fcl files in the source area of the currently 
#          setup version of uboonecode (that is, under 
#          $UBOONECODE_DIR/source/fcl/gen), and makes a corresponding xml
#          project file in the local directory.
#
# Usage:
#
# make_xml_mcc5.0.sh [-h|--help] [-r <release>] [-u|--user <user>] [--local <dir|tar>] [--nev <n>] [--nevjob <n>] [--nevgjob <n>]
#
# Options:
#
# -h|--help     - Print help.
# -r <release>  - Use the specified larsoft/uboonecode release.
# -u|--user <user> - Use users/<user> as working and output directories
#                    (default is to use uboonepro).
# --local <dir|tar> - Specify larsoft local directory or tarball (xml 
#                     tag <local>...</local>).
# --nev <n>     - Specify number of events for all samples.  Otherwise
#                 use hardwired defaults.
# --nevjob <n>  - Specify the default number of events per job.
# --nevgjob <n> - Specify the maximum number of events per gen/g4 job.
#
#----------------------------------------------------------------------

# Parse arguments.

rel=v02_05_01
userdir=uboonepro
userbase=$userdir
nevarg=0
nevjob=100
nevgjobarg=0
local=''

while [ $# -gt 0 ]; do
  case "$1" in

    # User directory.

    -h|--help )
      echo "Usage: make_xml_mcc5.0.sh [-h|--help] [-r <release>] [-u|--user <user>] [--local <dir|tar>] [--nev <n>] [--nevjob <n>] [--nevgjob <n>]"
      exit
    ;;

    # Release.

    -r )
    if [ $# -gt 1 ]; then
      rel=$2
      shift
    fi
    ;;

    # User.

    -u|--user )
    if [ $# -gt 1 ]; then
      userdir=users/$2
      userbase=$2
      shift
    fi
    ;;

    # Local release.

    --local )
    if [ $# -gt 1 ]; then
      local=$2
      shift
    fi
    ;;

    # Total number of events.

    --nev )
    if [ $# -gt 1 ]; then
      nevarg=$2
      shift
    fi
    ;;

    # Number of events per job.

    --nevjob )
    if [ $# -gt 1 ]; then
      nevjob=$2
      shift
    fi
    ;;

    # Number of events per gen/g4 job.

    --nevgjob )
    if [ $# -gt 1 ]; then
      nevgjobarg=$2
      shift
    fi
    ;;

  esac
  shift
done

# Delete existing xml files.

rm -f *.xml

find $UBOONECODE_DIR/source/fcl/gen -name \*.fcl | while read fcl
do
  if ! echo $fcl | grep -q common; then
    newprj=`basename $fcl .fcl`
    newxml=${newprj}.xml
    filt=1

    # Make xml file.

    echo "Making ${newprj}.xml"

    # Generator

    genfcl=`basename $fcl`

    # G4

    g4fcl=standard_g4_uboone.fcl

    # Detsim (optical + tpc).

    detsimfcl=standard_detsim_uboone.fcl
    if echo $newprj | grep -q dirt; then
      detsimfcl=standard_detsim_uboone_tpcfilt.fcl
      if echo $newprj | grep -q bnb; then
        filt=5
      else
        filt=20
      fi
    fi

    # Reco 2D

    reco2dfcl=standard_reco_uboone_2D.fcl

    # Reco 3D

    reco3dfcl=standard_reco_uboone_3D.fcl

    # Merge/Analysis

    mergefcl=standard_ana_uboone.fcl

    # Set number of gen/g4 events per job.

    nevgjob=$nevgjobarg
    if [ $nevgjob -eq 0 ]; then
      if echo $newprj | grep -q dirt; then
        if echo $newprj | grep -q cosmic; then
          nevgjob=200
        else
          nevgjob=2000
        fi
      else
        nevgjob=nevjob
      fi
    fi

    # Set number of events.

    nev=$nevarg
    if [ $nev -eq 0 ]; then
      if [ $newprj = prodgenie_bnb_nu_cosmic_uboone ]; then
        nev=50000
      elif [ $newprj = prodgenie_bnb_nu_uboone ]; then
        nev=20000
      elif [ $newprj = prodgenie_bnb_nue_cosmic_uboone ]; then
        nev=20000
      elif [ $newprj = prodgenie_bnb_nue_uboone ]; then
        nev=20000
      elif [ $newprj = prodcosmics_uboone ]; then
        nev=20000
      else
        nev=10000
      fi
    fi
    nev=$(( $nev * $filt ))

    # Calculate the number of worker jobs.

    njob1=$(( $nev / $nevgjob ))         # Pre-filter (gen, g4)
    njob2=$(( $nev / $nevjob / $filt ))  # Post-filter (detsim and later)
    if [ $njob1 -lt $njob2 ]; then
      njob1=$njob2
    fi

  cat <<EOF > $newxml
<?xml version="1.0"?>

<!-- Production Project -->

<!DOCTYPE project [
<!ENTITY release "$rel">
<!ENTITY file_type "mc">
<!ENTITY run_type "physics">
<!ENTITY name "$newprj">
<!ENTITY tag "mcc5.0">
]>

<project name="&name;">

  <!-- Group -->
  <group>uboone</group>

  <!-- Project size -->
  <numevents>$nev</numevents>

  <!-- Operating System -->
  <os>SL5,SL6</os>

  <!-- Batch resources -->
  <resource>DEDICATED,OPPORTUNISTIC</resource>

  <!-- Larsoft information -->
  <larsoft>
    <tag>&release;</tag>
    <qual>e5:prof</qual>
EOF
  echo "local=$local"
  if [ x$local != x ]; then
    echo "    <local>${local}</local>" >> $newxml
  fi
  cat <<EOF >> $newxml
  </larsoft>

  <!-- Project stages -->

  <stage name="gen">
    <fcl>$genfcl</fcl>
    <outdir>/pnfs/uboone/scratch/${userdir}/&release;/gen/&name;</outdir>
    <workdir>/uboone/app/users/${userbase}/&release;/gen/&name;</workdir>
    <numjobs>$njob1</numjobs>
    <datatier>generated</datatier>
    <defname>&name;_&tag;_gen</defname>
  </stage>

  <stage name="g4">
    <fcl>$g4fcl</fcl>
    <outdir>/pnfs/uboone/scratch/${userdir}/&release;/g4/&name;</outdir>
    <workdir>/uboone/app/users/${userbase}/&release;/g4/&name;</workdir>
    <numjobs>$njob1</numjobs>
    <datatier>simulated</datatier>
    <defname>&name;_&tag;_g4</defname>
  </stage>

EOF
  if [ x$detsimfcl != x ]; then
    cat <<EOF >> $newxml
  <stage name="detsim">
    <fcl>$detsimfcl</fcl>
    <outdir>/pnfs/uboone/scratch/${userdir}/&release;/detsim/&name;</outdir>
    <workdir>/uboone/app/users/${userbase}/&release;/detsim/&name;</workdir>
    <numjobs>$njob2</numjobs>
    <datatier>detector-simulated</datatier>
    <defname>&name;_&tag;_detsim</defname>
  </stage>

EOF
  fi
  if [ x$optsimfcl != x ]; then
    cat <<EOF >> $newxml
  <stage name="optsim">
    <fcl>$optsimfcl</fcl>
    <outdir>/pnfs/uboone/scratch/${userdir}/&release;/optsim/&name;</outdir>
    <workdir>/uboone/app/users/${userbase}/&release;/optsim/&name;</workdir>
    <numjobs>$njob2</numjobs>
    <datatier>optical-simulated</datatier>
    <defname>&name;_&tag;_optsim</defname>
  </stage>

EOF
  fi
  if [ x$tpcsimfcl != x ]; then
    cat <<EOF >> $newxml
  <stage name="tpcsim">
    <fcl>$tpcsimfcl</fcl>
    <outdir>/pnfs/uboone/scratch/${userdir}/&release;/tpcsim/&name;</outdir>
    <workdir>/uboone/app/users/${userbase}/&release;/tpcsim/&name;</workdir>
    <numjobs>$njob2</numjobs>
    <datatier>tpc-simulated</datatier>
    <defname>&name;_&tag;_tpcsim</defname>
  </stage>

EOF
  fi
  cat <<EOF >> $newxml
  <stage name="reco2D">
    <fcl>$reco2dfcl</fcl>
    <outdir>/pnfs/uboone/scratch/${userdir}/&release;/reco2D/&name;</outdir>
    <workdir>/uboone/app/users/${userbase}/&release;/reco2D/&name;</workdir>
    <numjobs>$njob2</numjobs>
    <datatier>reconstructed-2d</datatier>
    <defname>&name;_&tag;_reco2D</defname>
  </stage>

  <stage name="reco3D">
    <fcl>$reco3dfcl</fcl>
    <outdir>/pnfs/uboone/scratch/${userdir}/&release;/reco3D/&name;</outdir>
    <workdir>/uboone/app/users/${userbase}/&release;/reco3D/&name;</workdir>
    <numjobs>$njob2</numjobs>
    <datatier>reconstructed-3d</datatier>
    <defname>&name;_&tag;_reco3D</defname>
  </stage>

  <stage name="mergeana">
    <fcl>$mergefcl</fcl>
    <outdir>/pnfs/uboone/scratch/${userdir}/&release;/mergeana/&name;</outdir>
    <workdir>/uboone/app/users/${userbase}/&release;/mergeana/&name;</workdir>
    <numjobs>$njob2</numjobs>
    <targetsize>8000000000</targetsize>
    <datatier>reconstructed-3d</datatier>
    <defname>&name;_&tag;</defname>
  </stage>

  <!-- file type -->
  <filetype>&file_type;</filetype>

  <!-- run type -->
  <runtype>&run_type;</runtype>

</project>
EOF

  fi

done
