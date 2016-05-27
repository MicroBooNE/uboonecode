////////////////////////////////////////////////////////////////////////
// \file SpaceChargeServiceMicroBooNE.h
//
// \brief header of service for storing/accessing space charge distortions for MicroBooNE
//
// \author mrmooney@bnl.gov
// 
////////////////////////////////////////////////////////////////////////
#ifndef SPACECHARGESERVICEMICROBOONE_H
#define SPACECHARGESERVICEMICROBOONE_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Principal/Run.h"
#include "uboone/SpaceCharge/SpaceChargeMicroBooNE.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"


namespace spacecharge{
  class SpaceChargeServiceMicroBooNE : public SpaceChargeService {
    public:
      
      // this enables art to print the configuration help:
      //using Parameters = art::ServiceTable<spacecharge::SpaceChargeMicroBooNE::ConfigurationParameters_t>;
      
      SpaceChargeServiceMicroBooNE(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

      virtual void   reconfigure(fhicl::ParameterSet const& pset) override;
      void   preBeginRun(const art::Run& run);


      virtual const  provider_type* provider() const override { return fProp.get();}

    private:

      std::unique_ptr<spacecharge::SpaceChargeMicroBooNE> fProp;

    }; // class SpaceChargeServiceMicroBooNE
} //namespace spacecharge
DECLARE_ART_SERVICE_INTERFACE_IMPL(spacecharge::SpaceChargeServiceMicroBooNE, spacecharge::SpaceChargeService, LEGACY)
#endif // SPACECHARGESERVICEMICROBOONE_H
