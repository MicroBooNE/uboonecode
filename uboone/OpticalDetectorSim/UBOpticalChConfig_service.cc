#ifndef UBOPTICALCHCONFIG_CXX
#define UBOPTICALCHCONFIG_CXX

#include "UBOpticalChConfig.h"
#include "Utilities/LArProperties.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace opdet {

  //---------------------------------------------------------------------------------
  UBOpticalChConfig::UBOpticalChConfig(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
  //---------------------------------------------------------------------------------
  {
    this->reconfigure(pset);
  }
  
  //-----------------------------------------------------------
  void UBOpticalChConfig::reconfigure(fhicl::ParameterSet const& pset)
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

    
    // Correct QE by prescaling set in LArProperties
    art::ServiceHandle<util::LArProperties>   LarProp;
    for (unsigned int i = 0; i < fParams.at(kQE).size(); i++) {

      if ( LarProp->ScintPreScale() > fParams.at(kQE)[i]) {
        fParams.at(kQE)[i] /= LarProp->ScintPreScale();
      }
      else {
        mf::LogError("UBOpticalChConfig_service") << "Quantum efficiency set in UBOpticalChConfig_service, "
                                                  << fParams.at(kQE)[i]
                                                  << " is too large.  It is larger than the prescaling applied during simulation, "
                                                  << LarProp->ScintPreScale()
                                                  << ".  Final QE must be equalt to or smaller than the QE applied at simulation time.";
        assert(false);
      }
    }

    art::ServiceHandle<geo::Geometry> geom;
    for(size_t i=0; i<kChConfigTypeMax; ++i)
      if(fParams.at(i).size() != geom->NOpChannels())
	throw UBOpticalException(Form("ChConfigType_t %zu # values (%zu) != # channels (%d)!",
				      i,
				      fParams.at(i).size(),
				      geom->NOpChannels()));
    
  }


  DEFINE_ART_SERVICE(UBOpticalChConfig)
    
}
#endif
