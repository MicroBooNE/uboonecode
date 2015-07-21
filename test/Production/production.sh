#! /bin/bash

# Save a copy of the environment (for debugging).

env > env.txt

# Manually add larpandora to FW_SEARCH_PATH to work around lar-ci bug.

export FW_SEARCH_PATH=${LARPANDORA_DIR}/scripts:$FW_SEARCH_PATH

# This script runs the full mc+reco chain using standard released fcl files.

input=''
for fcl in prod_muminus_0.5-5.0GeV_25degf_uboone.fcl standard_g4_uboone.fcl standard_detsim_uboone.fcl reco_uboone_stage_1.fcl reco_uboone_stage_2_w_cluster3d.fcl standard_ana_uboone.fcl
do
  output=`basename $fcl .fcl`.root
  out=`basename $fcl .fcl`.out
  err=`basename $fcl .fcl`.err
  if [ x$input = x ]; then
    cmd="lar --rethrow-all -c $fcl -o $output -n 5"
  else
    cmd="lar --rethrow-all -c $fcl -s $input -o $output -n 5"
  fi
  echo $cmd
  $cmd > $out 2> $err
  stat=$?
  echo "Command finished with status $stat"
  if [ $stat -ne 0 ]; then
    exit $stat
  fi
  input=$output
done
