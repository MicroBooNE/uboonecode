#ifndef UBOONEPMTGAINSERVICE_CC
#define UBOONEPMTGAINSERVICE_CC

#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "larevt/CalibrationDBI/Interface/PmtGainService.h"
#include "larevt/CalibrationDBI/Providers/SIOVPmtGainProvider.h"
#include "UbooneCalibrationServiceHelper.h"

namespace lariov{

  /**
     \class UboonePmtGainService
     art service implementation of PmtGainService.  Implements 
     a pmt gain retrieval service for database scheme in which 
     all elements in a database folder share a common interval of validity
  */
  class UboonePmtGainService : public PmtGainService {
  
    public:
    
      UboonePmtGainService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
      ~UboonePmtGainService(){}
      
      void PreProcessEvent(const art::Event& evt);
     
    private:
    
      PmtGainProvider const& DoGetProvider() const override {
        return fProvider;
      }   
      
      PmtGainProvider const* DoGetProviderPtr() const override {
        return &fProvider; 
      }
    
      SIOVPmtGainProvider fProvider;
      UbooneCalibrationServiceHelper fHelper;
  };
}//end namespace lariov
      
DECLARE_ART_SERVICE_INTERFACE_IMPL(lariov::UboonePmtGainService, lariov::PmtGainService, LEGACY)
      

namespace lariov{

  UboonePmtGainService::UboonePmtGainService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg) 
  : fProvider(pset.get<fhicl::ParameterSet>("PmtGainProvider")),
    fHelper(pset.get<fhicl::ParameterSet>("CalibrationHelper"))
  {
    //register callback to update local database cache before each event is processed
    reg.sPreProcessEvent.watch(this, &UboonePmtGainService::PreProcessEvent);
  }
  
  
  void UboonePmtGainService::PreProcessEvent(const art::Event& evt) {
    
    fProvider.Update( fHelper.GetTimeStamp(evt, "PMT Gain") );
  } 
  
}//end namespace lariov

DEFINE_ART_SERVICE_INTERFACE_IMPL(lariov::UboonePmtGainService, lariov::PmtGainService)

#endif
