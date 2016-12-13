#ifndef CRTDetSim_HH_
#define CRTDetSim_HH_


#include "art/Framework/Core/ModuleMacros.h" // For DEFINE_ART_MODULE
#include "fhiclcpp/ParameterSet.h" // for ParameterSet

#include "art/Framework/Core/EDProducer.h" // Base Class
#include "art/Framework/Principal/Event.h" // For produce

#include <string>

namespace crt{
  class CRTDetSim :  public art:: EDProducer{
    uint32_t fThreshold;
    std::string fProducerName; //Where to pull te IDEs from
  public:

    /// Default ctor
    CRTDetSim(const fhicl::ParameterSet&);

    /// Default dtor
    ~CRTDetSim();

    /// art::EDProducer::produce implementation
    virtual void produce (art::Event&);

  };
}


#endif  //CRTDetSim_HH_
