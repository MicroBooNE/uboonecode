#!/bin/bash
#Test LArSoft code with "prodsingle_uboone.fcl".
# echurch@fnal.gov

#cp  ${UBOONECODE_DIR}/job/standard_detsim_uboone.fcl .
#echo "services.FileCatalogMetadataExtras.RenameTemplate: '' " >> ./standard_detsim_uboone.fc
#echo "outputs.out1.fileName: 'openclose_detsim_uboone.root'" >> ./standard_detsim_uboone.fcl

lar --process-name citest-detsim -c standard_detsim_uboone.fcl -s $1 -n 1 -o openclose_detsim_uboone.root -T openclose_detsim_hist_uboone.root
