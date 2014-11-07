#!/bin/bash
#Test LArSoft code with "prodsingle_uboone.fcl".
# echurch@fnal.gov


cp  ${UBOONECODE_DIR}/job/standard_detsim_uboone.fcl .
#echo "services.user.FileCatalogMetadataExtras.RenameTemplate: '' " >> ./standard_detsim_uboone.fcl
echo "outputs.out1.fileName: 'hitana_uboone_detsim.root'" >> ./standard_detsim_uboone.fcl .
lar -c ./standard_detsim_uboone.fcl -s ../lar_ci_hitana_g4/hitana_uboone_g4.root -n -1 -T hitana_uboone_detsim_hist.root -o hitana_uboone_detsim.root 
