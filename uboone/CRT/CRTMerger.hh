#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"



namespace crt{
  class CRTMerger : public art::EDProducer {
  public:
    explicit CRTMerger( fhicl::ParameterSet const &pset );
    virtual ~CRTMerger();
    
    void reconfigure( fhicl::ParameterSet const &pset ); // override;
    void produce( art::Event &evt ) override;
    void beginRun(art::Run& run);
    
  private:
    
    std::string fSwizzlerProducerLabel; 

  };
  DEFINE_ART_MODULE( CRTMerger )
}