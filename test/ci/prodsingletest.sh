#!/bin/bash
#Test LArSoft code with "prodsingle.fcl".
#By Mark Dykstra, updated 6/13/2014


strace -o lar.strace lar -c ${LARSIM_DIR}/job/prodsingle.fcl 

