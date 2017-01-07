#ifndef UBOONEELECTRONICSCALIBSERVICE_CC
#define UBOONEELECTRONICSCALIBSERVICE_CC

#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "larevt/CalibrationDBI/Interface/ElectronicsCalibService.h"
#include "UbooneElectronicsCalibProvider.h"
#include "UbooneCalibrationServiceHelper.h"

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
      UbooneCalibrationServiceHelper fHelper;
  };
}//end namespace lariov
      
DECLARE_ART_SERVICE_INTERFACE_IMPL(lariov::UbooneElectronicsCalibService, lariov::ElectronicsCalibService, LEGACY)
      

namespace lariov{

  UbooneElectronicsCalibService::UbooneElectronicsCalibService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg) 
  : fProvider(pset.get<fhicl::ParameterSet>("ElectronicsCalibProvider")),
    fHelper(pset.get<fhicl::ParameterSet>("CalibrationHelper"))
  {
    //register callback to update local database cache before each event is processed
    reg.sPreProcessEvent.watch(this, &UbooneElectronicsCalibService::PreProcessEvent);
  }
  
  void UbooneElectronicsCalibService::PreProcessEvent(const art::Event& evt) {
    
    fProvider.Update( fHelper.GetTimeStamp(evt, "ASIC Calibrations") );
  } 

}//end namespace lariov

DEFINE_ART_SERVICE_INTERFACE_IMPL(lariov::UbooneElectronicsCalibService, lariov::ElectronicsCalibService)

#endif
