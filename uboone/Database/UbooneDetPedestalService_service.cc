#ifndef UBOONEDETPEDESTALSERVICE_CC
#define UBOONEDETPEDESTALSERVICE_CC

#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Providers/DetPedestalRetrievalAlg.h"
#include "UbooneCalibrationServiceHelper.h"

namespace lariov{

  /**
     \class UbooneDetPedestalService
     art service implementation of DetPedestalService.  Implements 
     a detector pedestal retrieval service for database scheme in which 
     all elements in a database folder share a common interval of validity
  */
  class UbooneDetPedestalService : public DetPedestalService {
  
    public:
    
      UbooneDetPedestalService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
      ~UbooneDetPedestalService(){}
      
      void PreProcessEvent(const art::Event& evt);
     
    private:
    
      const DetPedestalProvider& DoGetPedestalProvider() const override {
        return fProvider;
      }   
      
      UbooneCalibrationServiceHelper fHelper; 
    
      DetPedestalRetrievalAlg fProvider;
  };
}//end namespace lariov
      
DECLARE_ART_SERVICE_INTERFACE_IMPL(lariov::UbooneDetPedestalService, lariov::DetPedestalService, LEGACY)
      

namespace lariov{

  UbooneDetPedestalService::UbooneDetPedestalService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg) 
  : fHelper(pset.get<fhicl::ParameterSet>("CalibrationHelper")),
    fProvider(pset.get<fhicl::ParameterSet>("DetPedestalRetrievalAlg"))
  {
    //register callback to update local database cache before each event is processed
    reg.sPreProcessEvent.watch(this, &UbooneDetPedestalService::PreProcessEvent);
        
  }

  void UbooneDetPedestalService::PreProcessEvent(const art::Event& evt) {
    
    if (evt.isRealData() && evt.run() < 183) {
      fProvider.Update(1430000000000000000);
    }
    else fProvider.Update(fHelper.GetTimeStamp(evt, "Detector Pedestals"));
  }
    
}//end namespace lariov

DEFINE_ART_SERVICE_INTERFACE_IMPL(lariov::UbooneDetPedestalService, lariov::DetPedestalService)

#endif
