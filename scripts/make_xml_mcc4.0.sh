#! /bin/bash
#----------------------------------------------------------------------
#
# Name: make_xml_mcc4.0.sh
#
# Purpose: Make xml files for mcc 4.0.  This script loops over all
#          generator-level fcl files in the currently setup version
#          if ubfcl (that is, undert $UBFCL_DIR/gen), and makes a 
#          corresponding xml project file in $UBXML_DIR/mcc4.0.
#
# Usage:
#
# make_xml_mcc4.0.sh [--user <user>] [--nev <n>] [--njob <n>]
#
# Options:
#
# --user <user> - Use users/<user> as working and output directories
#                 (default is to use uboonepro).
# --nev <n>     - Specify number of events for all samples.  Otherwise
#                 use hardwired defaults.
# --njob <n>    - Specify the number of worker jobs.
#
#----------------------------------------------------------------------

# Parse arguments.

userdir=uboonepro
nevarg=0
njobarg=0

while [ $# -gt 0 ]; do
  case "$1" in

    # User directory.

    --user )
    if [ $# -gt 1 ]; then
      userdir=users/$2
      shift
    fi
    ;;

    # Number of events.

    --nev )
    if [ $# -gt 1 ]; then
      nevarg=$2
      shift
    fi
    ;;

    # Number of worker jobs.

    --njob )
    if [ $# -gt 1 ]; then
      njobarg=$2
      shift
    fi
    ;;

  esac
  shift
done

# Make sure ubfcl and ubxml are setup.

if [ x$UBFCL_DIR = x ]; then
  echo "Please setup ubfcl."
  exit 1
fi
if [ x$UBXML_DIR = x ]; then
  echo "Please setup ubxml."
  exit 1
fi

newxmldir=$UBXML_DIR/mcc4.0
mkdir -p $newxmldir
rm -f $newxmldir/*.xml

find $UBFCL_DIR/gen -name \*.fcl | while read fcl
do
  if ! echo $fcl | grep -q common; then
    newprj=`basename $fcl .fcl`
    newxml=$newxmldir/${newprj}.xml

    # Make xml file.

    echo "Making ${newprj}.xml"

    # Generator

    genfcl=`echo $fcl | sed "s;$UBFCL_DIR/;;"`

    # G4

    g4fcl=g4/standard_g4_uboone.fcl

    # Detsim (optical + tpc).

    detsimfcl=''
    #if echo $newprj | grep -q 3window; then
    #  detsimfcl=detsim/standard_detsim_3window_uboone.fcl
    #fi
    #echo "Using ${detsimfcl}."

    # Detsim optical

    optsimfcl=''
    #if echo $newprj | grep -q 3window; then
    #  optsimfcl=detsim/standard_detsim_3window_uboone_optical.fcl
    #fi
    #echo "Using ${optsimfcl}."

    # Detsim tpc

    tpcsimfcl=detsim/standard_detsim_uboone.fcl
    if echo $newprj | grep -q 3window; then
      tpcsimfcl=detsim/standard_detsim_3window_uboone_tpc.fcl
    fi
    echo "Using ${tpcsimfcl}."

    # Reco 2D

    reco2dfcl=reco/standard_reco_uboone_2D_noopt_nowires.fcl
    if echo $newprj | grep -q cosmic; then
      reco2dfcl=reco/standard_reco_uboone_2D_noopt_nowires_cosmic.fcl
    fi
    echo "Using ${reco2dfcl}."

    # Reco 3D

    reco3dfcl=reco/standard_reco_uboone_3D_noopt.fcl
    if echo $newprj | grep -q cosmic; then
      reco3dfcl=reco/standard_reco_uboone_3D_noopt_cosmic.fcl
    fi
    echo "Using ${reco3dfcl}."

    # Merge

    mergefcl=utility/copy.fcl

    nev=$nevarg
    njob=$njobarg

    # Set default number of events.

    if [ $nev -eq 0 ]; then
      if [ $newprj = prodgenie_bnb_nu_cosmic_3window_uboone ]; then
        nev=50000
      elif [ $newprj = prodgenie_bnb_nu_3window_uboone ]; then
        nev=20000
      else
        nev=10000
      fi
    fi

    # Set default number of workers, assuming 100 events/worker.

    if [ $njob -eq 0 ]; then
      njob=$(( $nev / 100 ))
    fi

  cat <<EOF > $newxml
<?xml version="1.0"?>

<!-- Production Project -->

<!DOCTYPE project [
<!ENTITY release "v1_00_04_ub01">
<!ENTITY ubfcl_version "v2_0_6">
<!ENTITY ubxml_version "v2_0_7">
<!ENTITY file_type "mc">
<!ENTITY run_type "physics">
<!ENTITY name "$newprj">
<!ENTITY tag "mcc4.0">
]>

<project name="&name;">

  <!-- Group -->
  <group>uboone</group>

  <!-- Project size -->
  <numevents>$nev</numevents>

  <!-- Operating System -->
  <os>SL5,SL6</os>

  <!-- Larsoft information -->
  <larsoft>
    <tag>&release;</tag>
    <qual>e4:prof</qual>
  </larsoft>

  <!-- Project stages -->

  <stage name="gen">
    <fcl>$genfcl</fcl>
    <outdir>/uboone/data/${userdir}/gen/&release;/&name;</outdir>
    <workdir>/uboone/app/users/${userdir}/gen/&release;/&name;</workdir>
    <numjobs>$njob</numjobs>
    <datatier>generated</datatier>
    <defname>&name;_&tag;_gen</defname>
  </stage>

  <stage name="g4">
    <fcl>$g4fcl</fcl>
    <outdir>/uboone/data/${userdir}/g4/&release;/&name;</outdir>
    <workdir>/uboone/app/users/${userdir}/g4/&release;/&name;</workdir>
    <numjobs>$njob</numjobs>
    <datatier>simulated</datatier>
    <defname>&name;_&tag;_g4</defname>
  </stage>

EOF
  if [ x$detsimfcl != x ]; then
    cat <<EOF >> $newxml
  <stage name="detsim">
    <fcl>$detsimfcl</fcl>
    <outdir>/uboone/data/${userdir}/detsim/&release;/&name;</outdir>
    <workdir>/uboone/app/users/${userdir}/detsim/&release;/&name;</workdir>
    <numjobs>$njob</numjobs>
    <datatier>optical-simulated</datatier>
    <defname>&name;_&tag;_detsim</defname>
  </stage>

EOF
  fi
  if [ x$optsimfcl != x ]; then
    cat <<EOF >> $newxml
  <stage name="optsim">
    <fcl>$optsimfcl</fcl>
    <outdir>/uboone/data/${userdir}/optsim/&release;/&name;</outdir>
    <workdir>/uboone/app/users/${userdir}/optsim/&release;/&name;</workdir>
    <numjobs>$njob</numjobs>
    <datatier>optical-simulated</datatier>
    <defname>&name;_&tag;_optsim</defname>
  </stage>

EOF
  fi
  if [ x$tpcsimfcl != x ]; then
    cat <<EOF >> $newxml
  <stage name="tpcsim">
    <fcl>$tpcsimfcl</fcl>
    <outdir>/uboone/data/${userdir}/tpcsim/&release;/&name;</outdir>
    <workdir>/uboone/app/users/${userdir}/tpcsim/&release;/&name;</workdir>
    <numjobs>$njob</numjobs>
    <datatier>tpc-simulated</datatier>
    <defname>&name;_&tag;_tpcsim</defname>
  </stage>

EOF
  fi
  cat <<EOF >> $newxml
  <stage name="reco2D">
    <fcl>$reco2dfcl</fcl>
    <outdir>/uboone/data/${userdir}/reco2D/&release;/&name;</outdir>
    <workdir>/uboone/app/users/${userdir}/reco2D/&release;/&name;</workdir>
    <numjobs>$njob</numjobs>
    <datatier>reconstructed-2d</datatier>
    <defname>&name;_&tag;_reco2D</defname>
  </stage>

  <stage name="reco3D">
    <fcl>$reco3dfcl</fcl>
    <outdir>/uboone/data/${userdir}/reco3D/&release;/&name;</outdir>
    <workdir>/uboone/app/users/${userdir}/reco3D/&release;/&name;</workdir>
    <numjobs>$njob</numjobs>
    <datatier>reconstructed-3d</datatier>
    <defname>&name;_&tag;_reco3D</defname>
  </stage>

  <stage name="merge">
    <fcl>$mergefcl</fcl>
    <outdir>/uboone/data/${userdir}/reco/&release;/&name;</outdir>
    <workdir>/uboone/app/users/${userdir}/reco/&release;/&name;</workdir>
    <numjobs>$njob</numjobs>
    <targetsize>2000000000</targetsize>
    <datatier>reconstructed</datatier>
    <defname>&name;_&tag;</defname>
  </stage>

  <!-- ubfcl version -->
  <ubfcl>&ubfcl_version;</ubfcl>

  <!-- ubxml version -->
  <ubxml>&ubxml_version;</ubxml>

  <!-- file type -->
  <filetype>&file_type;</filetype>

  <!-- run type -->
  <runtype>&run_type;</runtype>

</project>
EOF

  fi

done
