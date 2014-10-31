#!/bin/bash
#Test LArSoft code with "prodsingle_uboone.fcl".
# echurch@fnal.gov


lar -c ${UBOONECODE_DIR}/job/hitfinder_ana_uboone.fcl -s $1 -n 1000 -T hitana_uboone_canonical_hist.root
