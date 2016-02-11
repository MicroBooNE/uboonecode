#ifndef UBOONEDETPEDESTALSERVICE_CC
#define UBOONEDETPEDESTALSERVICE_CC

#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "larevt/CalibrationDBI/Interface/IDetPedestalService.h"
#include "larevt/CalibrationDBI/Providers/DetPedestalRetrievalAlg.h"

namespace lariov{

  /**
     \class UbooneDetPedestalService
     art service implementation of IDetPedestalService.  Implements 
     a detector pedestal retrieval service for database scheme in which 
     all elements in a database folder share a common interval of validity
  */
  class UbooneDetPedestalService : public IDetPedestalService {
  
    public:
    
      UbooneDetPedestalService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
      ~UbooneDetPedestalService(){}
      
      void PreProcessEvent(const art::Event& evt) {
        // This is a temporary kludge to allow microboone to analyze early data 
	// which did not have a proper timestamp in the daq header.
        if (evt.isRealData() && evt.run() < 183) {
	  std::uint64_t kludge_stamp = 1430000000000000000; //yes, there really needs to be 16 zeroes
	  fProvider.Update(kludge_stamp);
        }
	else fProvider.Update(evt.time().value());
      }
     
    private:
    
      const IDetPedestalProvider& DoGetPedestalProvider() const override {
        return fProvider;
      }    
    
      DetPedestalRetrievalAlg fProvider;
  };
}//end namespace lariov
      
DECLARE_ART_SERVICE_INTERFACE_IMPL(lariov::UbooneDetPedestalService, lariov::IDetPedestalService, LEGACY)
      

namespace lariov{

  UbooneDetPedestalService::UbooneDetPedestalService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg) 
  : fProvider(pset.get<fhicl::ParameterSet>("DetPedestalRetrievalAlg"))
  {
    //register callback to update local database cache before each event is processed
    reg.sPreProcessEvent.watch(this, &UbooneDetPedestalService::PreProcessEvent);
  }

}//end namespace lariov

DEFINE_ART_SERVICE_INTERFACE_IMPL(lariov::UbooneDetPedestalService, lariov::IDetPedestalService)

#endif
