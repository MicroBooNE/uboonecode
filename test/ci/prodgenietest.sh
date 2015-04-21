#!/bin/bash
#Test LArSoft code with "prodgenie.fcl".


cp  ${UBOONECODE_DIR}/job/prodgenie_bnb_nu_uboone.fcl .

# This is not working, for whatever reason. So, instead use the normal mechanism to pull in inputfiles via test_runner's built-in ifdh'ing in ci_tests.cfg.
#source `ups setup $SETUP_UPS`
#source /grid/fermiapp/products/uboone/setup_uboone.sh;
#setup ifdhc 
# Just pick a few bnb flux files and put 'em here.
#ifdh cp /uboone/data/uboonebeam/bnb_gsimple_fluxes_02.28.2014_470/gsimple_microboone-470-onaxis_mc_nu_dummy_ntrd_4216_numintp_00002.root gsimple_microboone-470-onaxis_mc_nu_dummy_ntrd_4216_numintp_00002.root

export FW_SEARCH_PATH=${FW_SEARCH_PATH}:.
echo "physics.producers.generator.FluxFiles: ['$1'] " >> ./prodgenie_bnb_nu_uboone.fcl 
strace -o lar.strace lar -c ./prodgenie_bnb_nu_uboone.fcl -n 1 -o genie_gen.root -T genie_hist.root

