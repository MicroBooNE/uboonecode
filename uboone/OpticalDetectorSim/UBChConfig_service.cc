#ifndef UBCHCONFIG_CXX
#define UBCHCONFIG_CXX

#include "UBChConfig.h"

namespace opdet {

  //---------------------------------------------------------------------------------
  UBChConfig::UBChConfig(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
  //---------------------------------------------------------------------------------
  {
    this->reconfigure(pset);
  }
  
  //-----------------------------------------------------------
  void UBChConfig::reconfigure(fhicl::ParameterSet const& pset)
  //-----------------------------------------------------------
  {
    
    fParams.at(kPedestalMean)   = pset.get<std::vector<float> >("PedestalMean");
    fParams.at(kPedestalSpread) = pset.get<std::vector<float> >("PedestalSpread");
    fParams.at(kQE)             = pset.get<std::vector<float> >("QE");
    fParams.at(kHighGain)       = pset.get<std::vector<float> >("HighGain");
    fParams.at(kLowGain)        = pset.get<std::vector<float> >("LowGain");
    fParams.at(kGainSpread)     = pset.get<std::vector<float> >("GainSpread");
    fParams.at(kT0)             = pset.get<std::vector<float> >("T0");
    fParams.at(kT0Spread)       = pset.get<std::vector<float> >("T0Spread");
    fParams.at(kDarkRate)       = pset.get<std::vector<float> >("DarkRate");

    art::ServiceHandle<geo::Geometry> geom;
    for(size_t i=0; i<kChConfigTypeMax; ++i)
      if(fParams.at(i).size() != geom->NOpChannels())
	throw UBOpticalException(Form("ChConfigType_t %zu # values (%zu) != # channels (%d)!",
				      i,
				      fParams.at(i).size(),
				      geom->NOpChannels()));
    
  }


  DEFINE_ART_SERVICE(UBChConfig)
    
}
#endif
