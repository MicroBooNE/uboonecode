#ifndef UBOONEELECTRONICSCALIBSERVICE_CC
#define UBOONEELECTRONICSCALIBSERVICE_CC

#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "larevt/CalibrationDBI/Interface/ElectronicsCalibService.h"
#include "UbooneElectronicsCalibProvider.h"
#include "uboone/DataOverlay/DataOverlayProducts/EventMixingSummary.h"

namespace lariov{

  /**
     \class UbooneElectronicsCalibService
     art service implementation of ElectronicsCalibService.  Implements 
     an electronics calibration retrieval service for database scheme in which 
     all elements in a database folder share a common interval of validity
  */
  class UbooneElectronicsCalibService : public ElectronicsCalibService {
  
    public:
    
      UbooneElectronicsCalibService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
      ~UbooneElectronicsCalibService(){}
      
      void PreProcessEvent(const art::Event& evt);
     
    private:
    
      ElectronicsCalibProvider const& DoGetProvider() const override {
        return fProvider;
      }   
      
      ElectronicsCalibProvider const* DoGetProviderPtr() const override {
        return &fProvider; 
      }
    
      UbooneElectronicsCalibProvider fProvider;
      std::string                    fMixingModuleLabel;
  };
}//end namespace lariov
      
DECLARE_ART_SERVICE_INTERFACE_IMPL(lariov::UbooneElectronicsCalibService, lariov::ElectronicsCalibService, LEGACY)
      

namespace lariov{

  UbooneElectronicsCalibService::UbooneElectronicsCalibService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg) 
  : fProvider(pset.get<fhicl::ParameterSet>("ElectronicsCalibProvider")),
    fMixingModuleLabel(pset.get<std::string>("EventMixingModuleLabel"))
  {
    //register callback to update local database cache before each event is processed
    reg.sPreProcessEvent.watch(this, &UbooneElectronicsCalibService::PreProcessEvent);
  }
  
  void UbooneElectronicsCalibService::PreProcessEvent(const art::Event& evt) {
    
    art::Handle< std::vector<mix::EventMixingSummary> > eventMixingSummary;
    evt.getByLabel(fMixingModuleLabel, eventMixingSummary);
    if (eventMixingSummary.isValid() && eventMixingSummary->size()>0) {
      if (eventMixingSummary->size() > 1) {
        std::cout<<"  INFO: "<<eventMixingSummary->size()<<" EventMixingSummary objects"<<std::endl;
      }
      art::Timestamp time_stamp = eventMixingSummary->front().Timestamp();
      std::cout<<"Using EventMixingSummary timestamp to query ASICs calibration database: "<<time_stamp.value()<<std::endl;
      fProvider.Update(time_stamp.value());
    }
    else {
      std::cout<<"Using art::Event timestamp to query ASICs calibration database: "<<evt.time().value()<<std::endl;
      //First grab an update from the database
      fProvider.Update(evt.time().value());
    }
  } 

}//end namespace lariov

DEFINE_ART_SERVICE_INTERFACE_IMPL(lariov::UbooneElectronicsCalibService, lariov::ElectronicsCalibService)

#endif
