#!/bin/bash
#Test LArSoft code with "prodgenie.fcl".


cp  ${UBOONECODE_DIR}/job/prodgenie_uboone.fcl .
echo "services.user.FileCatalogMetadataExtras.RenameTemplate: '' " >> ./prodgenie_uboone.fcl
export FW_SEARCH_PATH=${FW_SEARCH_PATH}:/uboone/data

strace -o lar.strace lar -c ${UBOONECODE_DIR}/job/prodgenie_uboone.fcl -n 1 -o genie_gen.root -T genie_hist.root

