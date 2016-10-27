#!/bin/bash
#Test LArSoft code with "prodsingle.fcl".

# only try to strace if we have it..
strace() {
   if [ -x /usr/bin/strace ]
   then
       /usr/bin/strace "$@"
   else
       if [ "$1" = "-o" ]
       then
           shift
           shift
       fi
       "$@"
   fi
}

#cp  ${UBOONECODE_DIR}/job/prodsingle_uboone.fcl .
#echo "services.FileCatalogMetadataExtras.RenameTemplate: '' " >> ./prodsingle_uboone.fcl
strace -o lar.strace lar -c ${UBOONECODE_DIR}/job/prodsingle_uboone.fcl -n 1 -o single_gen.root -T single_hist.root

