#!/bin/bash
#Test LArSoft code with "prodsingle_uboone.fcl".
# echurch@fnal.gov



cp  ${UBOONECODE_DIR}/job/standard_reco_uboone_2D.fcl .
#echo "services.user.FileCatalogMetadataExtras.RenameTemplate: '' " >> ./standard_reco_uboone_2D.fcl
echo "outputs.out1.fileName: 'hitana_uboone_reco2D.root'" >> ./standard_reco_uboone_2D.fcl

# comment out below modules, as they are slow and we don't need them
sed -e '/fuzzy/s/^/#/g' -i standard_reco_uboone_2D.fcl 
sed -e '/ccluster/s/^/#/g' -i standard_reco_uboone_2D.fcl 
sed -e '/pandora/s/^/#/g' -i standard_reco_uboone_2D.fcl 
# This gives us back our needed terminating ] for reco module list
sed -e '/pandora \]/a \]' -i  standard_reco_uboone_2D.fcl
# Get rid of now-extraneous trailing comma
sed -e '/gaushit,/s//gaushit/g' -i standard_reco_uboone_2D.fcl 

lar -c ./standard_reco_uboone_2D.fcl -s ../lar_ci_hitana_detsim_uboonecode/hitana_uboone_detsim.root  -n -1 -T hitana_uboone_reco2D_hist.root -o hitana_uboone_reco2D.root 
