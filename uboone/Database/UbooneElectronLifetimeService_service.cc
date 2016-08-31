#ifndef UBOONEELECTRONLIFETIMESERVICE_CC
#define UBOONEELECTRONLIFETIMESERVICE_CC

#include "UbooneElectronLifetimeService.h"

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
