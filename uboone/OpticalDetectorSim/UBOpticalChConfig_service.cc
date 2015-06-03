#ifndef UBOPTICALCHCONFIG_CXX
#define UBOPTICALCHCONFIG_CXX

#include "UBOpticalChConfig.h"
#include "Utilities/LArProperties.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Geometry/Geometry.h" // larcore
#include "Geometry/ExptGeoHelperInterface.h" // larcore
#include "uboone/Geometry/ChannelMapUBooNEAlg.h" // uboonecode

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
  //-----------------------------------------------------------
    reconfigure( _pset );
  }

  //-----------------------------------------------------------
  void UBOpticalChConfig::reconfigure(fhicl::ParameterSet const& pset)
  //-----------------------------------------------------------
  {
    art::ServiceHandle<geo::Geometry> geom;
    std::shared_ptr< const geo::ChannelMapUBooNEAlg > chanmap = std::dynamic_pointer_cast< const geo::ChannelMapUBooNEAlg >( geom->GetChannelMapAlg() );
    
    std::vector< std::vector< float    > > tmp_float_params;
    std::vector< std::vector< uint16_t > > tmp_int_params;
    tmp_float_params.resize ( kChConfigTypeMax );
    tmp_int_params.resize   ( kChConfigTypeMax );
    tmp_float_params.at ( kPedestalMean   ) = pset.get<std::vector< float   > >("PedestalMean");
    tmp_float_params.at ( kPedestalSpread ) = pset.get<std::vector< float   > >("PedestalSpread");
    tmp_float_params.at ( kQE             ) = pset.get<std::vector< float   > >("QE");
    tmp_float_params.at ( kPMTGain        ) = pset.get<std::vector< float   > >("PMTGain");
    tmp_float_params.at ( kSplitterGain   ) = pset.get<std::vector< float   > >("PMTGain");
    tmp_float_params.at ( kGainSpread     ) = pset.get<std::vector< float   > >("GainSpread");
    tmp_float_params.at ( kT0             ) = pset.get<std::vector< float   > >("T0");
    tmp_float_params.at ( kT0Spread       ) = pset.get<std::vector< float   > >("T0Spread");
    tmp_float_params.at ( kDarkRate       ) = pset.get<std::vector< float   > >("DarkRate");
    tmp_int_params.at   ( kDisc0Threshold ) = pset.get<std::vector< uint16_t> >("Disc0Threshold");
    tmp_int_params.at   ( kDisc1Threshold ) = pset.get<std::vector< uint16_t> >("Disc1Threshold");
    tmp_int_params.at   ( kDisc3Threshold ) = pset.get<std::vector< uint16_t> >("Disc3Threshold");

    std::vector< unsigned int >  channel_list = pset.get<std::vector<unsigned int> >("ChannelList");

    // ------------------------------------------------------------------------------------------------------
    // sanity check: number of readout channels in geo service matches number of channels in parameters
    unsigned int nchannel_values = geom->NOpChannels();

    for(size_t i=0; i<kChConfigTypeMax; ++i) {

      size_t nchannel_input = tmp_float_params.at(i).size();
      if(!nchannel_input)
	nchannel_input = tmp_int_params.at(i).size();

      if(nchannel_input != nchannel_values) 
	throw UBOpticalException(Form("ChConfigType_t enum=%zu # values (%zu) != # channels (%d)!",
				      i,
				      nchannel_input,
				      nchannel_values));
    }
    
    // ------------------------------------------------------------------------------------------------------
    
    // Correct QE by prescaling set in LArProperties
    art::ServiceHandle<util::LArProperties>   LarProp;
    auto tmp_QE = tmp_float_params.at( kQE );
    for (unsigned int i = 0; i < tmp_QE.size(); i++) {

      unsigned int chnum = channel_list.at(i);
      if ( chanmap->GetChannelGain( chnum )==opdet::LogicChannel )
	continue; // skip QE check for logic channels
      
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

    // Make sure FEM related settings are all 12bit
    std::vector<ChConfigType_t> fem_config_v;
    fem_config_v.push_back(kDisc0Threshold);
    fem_config_v.push_back(kDisc1Threshold);
    fem_config_v.push_back(kDisc3Threshold);
    for(auto const& config_type : fem_config_v) {
      for(auto const& value : tmp_int_params.at(config_type)) {
	
	if( (value>>12) ) {
	  mf::LogError("UBOpticalChConfig_service") << "FEM configuration type: "
						    << config_type
						    << " must be 12bit integer (found values >= 4096)!";
	  assert(false);
	}
      }
    }
    
    // finally fll fParams data member
    for ( size_t i=0; i<kChConfigTypeMax; ++i) {
      if(tmp_int_params[i].size()) {
	int nch = 0;
	for ( auto channel : channel_list ) {
	  fIntParams[ (ChConfigType_t)i ][ channel ] = tmp_int_params[ i ][ nch ];
	  ++nch;
	}
      }
      if(tmp_float_params[i].size()) {
	int nch = 0;
	for ( auto channel : channel_list ) {
	  fFloatParams[ (ChConfigType_t)i ][ channel ] = tmp_float_params[ i ][ nch ];
	  ++nch;
	}
      }
    }
    
  }
  
  
  DEFINE_ART_SERVICE(UBOpticalChConfig)
  
}
#endif
