#!/bin/bash
#Test LArSoft code with "prodsingle_uboone.fcl".
# echurch@fnal.gov


cp  ${UBOONECODE_DIR}/job/prod_muminus_0.1-2.0GeV_isotropic_uboone.fcl .
echo "services.user.FileCatalogMetadataExtras.RenameTemplate: '' " >> ./prod_muminus_0.1-2.0GeV_isotropic_uboone.fcl

lar -c ./prod_muminus_0.1-2.0GeV_isotropic_uboone.fcl -n 100 -T hitana_uboone_prod_hist.root -o hitana_uboone_prod.root
