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
    _pset = pset;
    //this->reconfigure(pset);
    //fInitialized;
  }
  
  //-----------------------------------------------------------
  void UBOpticalChConfig::doInitialization() {
    reconfigure( _pset );
  }


  //-----------------------------------------------------------
  void UBOpticalChConfig::reconfigure(fhicl::ParameterSet const& pset)
  //-----------------------------------------------------------
  {
    
    std::vector< std::vector<float> > tmp_params;
    tmp_params.resize( kChConfigTypeMax );
    tmp_params.at( kPedestalMean )   = pset.get<std::vector<float> >("PedestalMean");
    tmp_params.at( kPedestalSpread ) = pset.get<std::vector<float> >("PedestalSpread");
    tmp_params.at( kQE )             = pset.get<std::vector<float> >("QE");
    tmp_params.at( kPMTGain )        = pset.get<std::vector<float> >("PMTGain");
    tmp_params.at( kSplitterGain )   = pset.get<std::vector<float> >("PMTGain");
    tmp_params.at( kGainSpread )     = pset.get<std::vector<float> >("GainSpread");
    tmp_params.at( kT0 )             = pset.get<std::vector<float> >("T0");
    tmp_params.at( kT0Spread )       = pset.get<std::vector<float> >("T0Spread");
    tmp_params.at( kDarkRate )       = pset.get<std::vector<float> >("DarkRate");

    
    // Correct QE by prescaling set in LArProperties
    art::ServiceHandle<util::LArProperties>   LarProp;
    auto tmp_QE = tmp_params.at( kQE );
    for (unsigned int i = 0; i < tmp_QE.size(); i++) {
      
      if ( LarProp->ScintPreScale() > tmp_QE.at(i) ) {
        tmp_QE[i] /= LarProp->ScintPreScale();
      }
      else {
        mf::LogError("UBOpticalChConfig_service") << "Quantum efficiency set in UBOpticalChConfig_service, "
                                                  << tmp_QE[i]
                                                  << " is too large.  It is larger than the prescaling applied during simulation, "
                                                  << LarProp->ScintPreScale()
                                                  << ".  Final QE must be equalt to or smaller than the QE applied at simulation time.";
        assert(false);
      }
    }
    
    art::ServiceHandle<geo::Geometry> geom;
    
    // sanity check: number of readout channels in geo service matches number of channels in parameters
    unsigned int nchannel_values = geom->NOpChannels() + kLogicNChannel;

    for(size_t i=0; i<kChConfigTypeMax; ++i)
      if( tmp_params.at(i).size() != nchannel_values )
	throw UBOpticalException(Form("ChConfigType_t %zu # values (%zu) != # channels (%d)!",
				      i,
				      tmp_params.at(i).size(),
				      nchannel_values));

    // finally fll fParams data member
    for ( size_t i=0; i<kChConfigTypeMax; ++i) {
      std::map< unsigned int, float > pardata;
      unsigned int ich = 0;
      for ( unsigned int ipmt=0; ipmt<geom->NOpDets(); ipmt++ ) {
	for ( unsigned int isplit=0; isplit<geom->NOpHardwareChannels( ipmt ); isplit++ ) {
	  unsigned int channel = geom->OpChannel( ipmt, isplit );
	  //fParams[.insert( std::pair< unsigned int, float >( channel, tmp_params.at( i ) ) );
	  pardata[ channel ] = tmp_params.at( i ).at( ich );
	  ich++;
	}
      }
      fParams[ i ] = pardata;
    }
    
  }


  DEFINE_ART_SERVICE(UBOpticalChConfig)
    
}
#endif
