#ifndef UBOONEELECTRONLIFETIMESERVICE_H
#define UBOONEELECTRONLIFETIMESERVICE_H

#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
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
        fProvider.Update( (DBTimeStamp_t)evt.run() );
      }
      
      const UbooneElectronLifetimeProvider& GetProvider() const {
        return fProvider;
      }
     
    private:

      UbooneElectronLifetimeProvider fProvider;
  };
}//end namespace lariov
      
DECLARE_ART_SERVICE(lariov::UbooneElectronLifetimeService, LEGACY)

#endif
