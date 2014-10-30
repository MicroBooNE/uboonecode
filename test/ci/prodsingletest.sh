#!/bin/bash
#Test LArSoft code with "prodsingle.fcl".
#By Mark Dykstra, updated 6/13/2014

cp  ${UBOONECODE_DIR}/job/prodsingle_uboone.fcl .
echo "services.user.FileCatalogMetadataExtras.RenameTemplate: '' " >> ./prodsingle_uboone.fcl
strace -o lar.strace lar -c ${UBOONECODE_DIR}/job/prodsingle_uboone.fcl -n 1 -o single_gen.root -T single_hist.root

