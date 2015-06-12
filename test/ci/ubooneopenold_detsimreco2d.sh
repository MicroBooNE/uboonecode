#!/bin/bash
#Test LArSoft code with "prodsingle_uboone.fcl".
# echurch@fnal.gov

cp  ${UBOONECODE_DIR}/job/standard_reco_uboone_2D.fcl .
#echo "services.user.FileCatalogMetadataExtras.RenameTemplate: '' " >> ./standard_reco_uboone_2D.fcl
echo "outputs.out1.fileName: 'openclose_reco2D_uboone.root'" >> ./standard_reco_uboone_2D.fcl
lar --process-name citest-reco2D -c reco_uboone_stage_1.fcl -s ../lar_ci_openold_detsim_uboonecode/openclose_detsim_uboone.root -n 1 -o openclose_reco2D_uboone.root -T openclose_reco2D_hist_uboone.root
