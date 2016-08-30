#ifndef UBOONEELECTRONLIFETIMESERVICE_CC
#define UBOONEELECTRONLIFETIMESERVICE_CC

#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "UbooneElectronLifetimeProvider.h"

namespace lariov{

  /**
     \class UbooneElectronLifetimeService
     art service provider for electron lifetime.  Implements 
     an electron lifetime retrieval service for database scheme in which 
     all elements in a database folder share a common interval of validity
  */
  class UbooneElectronLifetimeService {
  
    public:
    
      UbooneElectronLifetimeService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
      ~UbooneElectronLifetimeService(){}
      
      void PreProcessEvent(const art::Event& evt) {
        fProvider.Update( evt.run() );
      }
      
      const UbooneElectronLifetimeProvider& GetElectronLifetimeProvider() const {
        return fProvider;
      }
     
    private:

      UbooneElectronLifetimeProvider fProvider;
  };
}//end namespace lariov
      
DECLARE_ART_SERVICE(lariov::UbooneElectronLifetimeService, LEGACY)
      

namespace lariov{

  UbooneElectronLifetimeService::UbooneElectronLifetimeService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg) 
  : fProvider(pset.get<fhicl::ParameterSet>("ElectronLifetimeProvider"))
  {
    //register callback to update local database cache before each event is processed
    //reg.sPreProcessEvent.watch(&UbooneElectronLifetimeService::PreProcessEvent, *this);
    reg.sPreProcessEvent.watch(this, &UbooneElectronLifetimeService::PreProcessEvent);
  }

}//end namespace lariov

DEFINE_ART_SERVICE(lariov::UbooneElectronLifetimeService)

#endif
