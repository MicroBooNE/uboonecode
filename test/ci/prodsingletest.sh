#!/bin/bash
#Test LArSoft code with "prodsingle.fcl".
#By Mark Dykstra, updated 6/13/2014


strace -o lar.strace lar -c ${UBOONECODE_DIR}/job/prodsingle_uboone.fcl 

