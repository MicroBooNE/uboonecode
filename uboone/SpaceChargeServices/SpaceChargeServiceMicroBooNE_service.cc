////////////////////////////////////////////////////////////////////////
// \file SpaceChargeMicroBooNE.cxx
//
// \brief implementation of class for storing/accessing space charge distortions for MicroBooNE
//
// \author mrmooney@bnl.gov
// 
////////////////////////////////////////////////////////////////////////

// C++ language includes
#include <iostream>

// LArSoft includes
#include "uboone/SpaceChargeServices/SpaceChargeServiceMicroBooNE.h"

// ROOT includes
#include "TMath.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

//-----------------------------------------------
spacecharge::SpaceChargeServiceMicroBooNE::SpaceChargeServiceMicroBooNE(fhicl::ParameterSet const& pset, art::ActivityRegistry &reg)
{
  fProp.reset(new spacecharge::SpaceChargeMicroBooNE(pset));

  reg.sPreBeginRun.watch(this, &SpaceChargeServiceMicroBooNE::preBeginRun);
}

//----------------------------------------------
void spacecharge::SpaceChargeServiceMicroBooNE::preBeginRun(const art::Run& run)
{
  fProp->Update(run.id().run());
}

//------------------------------------------------
void spacecharge::SpaceChargeServiceMicroBooNE::reconfigure(fhicl::ParameterSet const& pset)
{
  fProp->Configure(pset);  
  return;
}

//------------------------------------------------
DEFINE_ART_SERVICE_INTERFACE_IMPL(spacecharge::SpaceChargeServiceMicroBooNE, spacecharge::SpaceChargeService)
