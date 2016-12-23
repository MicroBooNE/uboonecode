#ifndef UBOONEDETPEDESTALSERVICE_CC
#define UBOONEDETPEDESTALSERVICE_CC

#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Providers/DetPedestalRetrievalAlg.h"
#include "uboone/DataOverlay/DataOverlayProducts/EventMixingSummary.h"

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
      
      std::string fMixingModuleLabel; 
    
      DetPedestalRetrievalAlg fProvider;
  };
}//end namespace lariov
      
DECLARE_ART_SERVICE_INTERFACE_IMPL(lariov::UbooneDetPedestalService, lariov::DetPedestalService, LEGACY)
      

namespace lariov{

  UbooneDetPedestalService::UbooneDetPedestalService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg) 
  : fMixingModuleLabel(pset.get<std::string>("EventMixingModuleLabel")),
    fProvider(pset.get<fhicl::ParameterSet>("DetPedestalRetrievalAlg"))
  {
    //register callback to update local database cache before each event is processed
    reg.sPreProcessEvent.watch(this, &UbooneDetPedestalService::PreProcessEvent);
        
  }

  void UbooneDetPedestalService::PreProcessEvent(const art::Event& evt) {
    
    art::Handle< std::vector<mix::EventMixingSummary> > eventMixingSummary;
    evt.getByLabel(fMixingModuleLabel, eventMixingSummary);
    if (eventMixingSummary.isValid() && eventMixingSummary->size()>0) {
      if (eventMixingSummary->size() > 1) {
        std::cout<<"  INFO: "<<eventMixingSummary->size()<<" EventMixingSummary objects"<<std::endl;
      }
      art::Timestamp time_stamp = eventMixingSummary->front().Timestamp();
      std::cout<<"Using EventMixingSummary timestamp to query pedestal database: "<<time_stamp.value()<<std::endl;
      fProvider.Update(time_stamp.value());
    }
    else {  
      // This is a temporary kludge to allow microboone to analyze early data 
      // which did not have a proper timestamp in the daq header.
      if (evt.isRealData() && evt.run() < 183) {
        std::uint64_t kludge_stamp = 1430000000000000000; //yes, there really needs to be 16 zeroes
        std::cout<<"Using kludged timestamp to query pedestal database: "<<kludge_stamp<<std::endl;
	fProvider.Update(kludge_stamp);
      }
      else {
        std::cout<<"Using art::Event timestamp to query pedestal database: "<<evt.time().value()<<std::endl;
        fProvider.Update(evt.time().value());
      }
    }
  }
    
}//end namespace lariov

DEFINE_ART_SERVICE_INTERFACE_IMPL(lariov::UbooneDetPedestalService, lariov::DetPedestalService)

#endif
