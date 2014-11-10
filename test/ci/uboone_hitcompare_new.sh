#!/bin/bash
#Test LArSoft code with "prodsingle_uboone.fcl".
# echurch@fnal.gov



lar -c ${UBOONECODE_DIR}/job/hitfinder_ana_uboone.fcl -s ../lar_ci_hitana_reco2D_uboonecode/hitana_uboone_reco2D.root  -n -1 -T hitana_uboone_new_hist.root
