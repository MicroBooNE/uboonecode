
#include "uboone/CRT/CRTMerger.hh"

namespace crt{

  void CRTMerger::reconfigure( fhicl::ParameterSet const &p ){
    fSwizzlerProducerLabel = p.get< std::string >( "SwizzlerProducerLabel" );
    return;
  }

  CRTMerger::CRTMerger( fhicl::ParameterSet const &pset ){
    this->reconfigure( pset );
    produces< std::vector<crt::CRTData> >();  
  }

  CRTMerger::~CRTMerger() {
  }

  void CRTMerger::beginRun(art::Run& run){
    
    //TODO: Fill
  }

  void CRTMerger::produce( art::Event &evt ){
  	//TODO: Fill
    return;    
  }


}
