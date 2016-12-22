#ifndef UBOONEPMTGAINSERVICE_CC
#define UBOONEPMTGAINSERVICE_CC

#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "larevt/CalibrationDBI/Interface/PmtGainService.h"
#include "larevt/CalibrationDBI/Providers/SIOVPmtGainProvider.h"
#include "uboone/DataOverlay/DataOverlayProducts/EventMixingSummary.h"

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
      std::string         fMixingModuleLabel;
  };
}//end namespace lariov
      
DECLARE_ART_SERVICE_INTERFACE_IMPL(lariov::UboonePmtGainService, lariov::PmtGainService, LEGACY)
      

namespace lariov{

  UboonePmtGainService::UboonePmtGainService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg) 
  : fProvider(pset.get<fhicl::ParameterSet>("PmtGainProvider")),
    fMixingModuleLabel(pset.get<std::string>("EventMixingModuleLabel"))
  {
    //register callback to update local database cache before each event is processed
    reg.sPreProcessEvent.watch(this, &UboonePmtGainService::PreProcessEvent);
  }
  
  
  void UboonePmtGainService::PreProcessEvent(const art::Event& evt) {
    
    art::Handle< std::vector<mix::EventMixingSummary> > eventMixingSummary;
    evt.getByLabel(fMixingModuleLabel, eventMixingSummary);
    if (eventMixingSummary.isValid() && eventMixingSummary->size()>0) {
      if (eventMixingSummary->size() > 1) {
        std::cout<<"  INFO: "<<eventMixingSummary->size()<<" EventMixingSummary objects"<<std::endl;
      }
      art::Timestamp time_stamp = eventMixingSummary->front().Timestamp();
      std::cout<<"Using EventMixingSummary timestamp to query PMT gain database: "<<time_stamp.value()<<std::endl;
      fProvider.Update(time_stamp.value());
    }
    else {
      std::cout<<"Using art::Event timestamp to query PMT Gain database: "<<evt.time().value()<<std::endl;
      //First grab an update from the database
      fProvider.Update(evt.time().value());
    }
  } 
  
}//end namespace lariov

DEFINE_ART_SERVICE_INTERFACE_IMPL(lariov::UboonePmtGainService, lariov::PmtGainService)

#endif
