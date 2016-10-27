#!/bin/bash
#Test LArSoft code with "prodsingle_uboone.fcl".
# echurch@fnal.gov

#cp  ${UBOONECODE_DIR}/job/standard_reco_uboone_3D.fcl .
#echo "services.FileCatalogMetadataExtras.RenameTemplate: '' " >> ./standard_reco_uboone_3D.fcl
#echo "outputs.out1.fileName: 'openclose_reco3D_uboone.root'" >> ./standard_reco_uboone_3D.fcl
lar --process-name citest-reco3D -c reco_uboone_stage_2.fcl -s ../lar_ci_openold_detsim2d_uboonecode/openclose_reco2D_uboone.root -n 1 -o openclose_reco3D_uboone.root -T openclose_reco3D_hist_uboone.root
