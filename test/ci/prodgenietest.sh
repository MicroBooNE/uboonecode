#!/bin/bash
#Test LArSoft code with "prodgenie.fcl".


cp  ${UBOONECODE_DIR}/job/prodgenie_uboone.fcl .
echo "services.user.FileCatalogMetadataExtras.RenameTemplate: '' " >> ./prodgenie_uboone.fcl

setup ifdhc
# Just pick a few bnb flux files and put 'em here.
ifdh cp /uboone/data/uboonebeam/bnb_gsimple_fluxes_02.28.2014_470/gsimple_microboone-470-onaxis_mc_nu_dummy_ntrd_4216_numintp_00002.root gsimple_microboone-470-onaxis_mc_nu_dummy_ntrd_4216_numintp_00002.root
ifdh cp /uboone/data/uboonebeam/bnb_gsimple_fluxes_02.28.2014_470/gsimple_microboone-470-onaxis_mc_nu_dummy_ntrd_4216_numintp_00003.root gsimple_microboone-470-onaxis_mc_nu_dummy_ntrd_4216_numintp_00003.root
ifdh cp /uboone/data/uboonebeam/bnb_gsimple_fluxes_02.28.2014_470/gsimple_microboone-470-onaxis_mc_nu_dummy_ntrd_4216_numintp_00004.root gsimple_microboone-470-onaxis_mc_nu_dummy_ntrd_4216_numintp_00004.root

strace -o lar.strace lar -c ${UBOONECODE_DIR}/job/prodgenie_uboone.fcl -n 1 -o genie_gen.root -T genie_hist.root

