#ifndef OVERLAY_DATAOVERLAY_OPDETWAVEFORMMIXER_CXX
#define OVERLAY_DATAOVERLAY_OPDETWAVEFORMMIXER_CXX

#include "OpDetWaveformMixer.h"
#include <limits>
#include <iostream>
#include <stdexcept>
#include <algorithm>


void mix::OpDetWaveformMixer::DeclareData(std::vector<raw::OpDetWaveform> const& dataVector,
					  std::vector<raw::OpDetWaveform> & outputVector){

  fChannelIndexMap.clear();
  outputVector.reserve(dataVector.size());

  for(auto const& od : dataVector){


    outputVector.emplace_back(od);

    if(od.size() < fMinSampleSize) continue;
    
    //we're going to keep the longest one ... JUST handling beam gate stuff for now
    auto my_channel = fChannelIndexMap.find(od.ChannelNumber());
    if( my_channel != fChannelIndexMap.end() && outputVector[my_channel->second].size() > od.size())
      continue;

    fChannelIndexMap[od.ChannelNumber()] = outputVector.size()-1;

  }
  
}

void mix::OpDetWaveformMixer::Mix(std::vector<raw::OpDetWaveform> const& mcVector,
				  std::unordered_map<raw::Channel_t,float> const& scale_map,
				  std::vector<raw::OpDetWaveform> & outputVector){


  for( auto const& od : mcVector){

    if(od.size() < fMinSampleSize) continue;

    auto it_ch = fChannelIndexMap.find(od.ChannelNumber());

    //if this channel is not in the data, skip this channel!
    if(it_ch==fChannelIndexMap.end())
      continue;

    size_t i_output = it_ch->second;

    fRDAdderAlg.SetPedestalInputs(2048,0.0); //HARDCODED PEDESTAL AT 2048!!!!!!!
    fRDAdderAlg.SetScaleInputs(scale_map.at(od.ChannelNumber()),1.0);
    
    //If the sizes are not the same...
    if(od.size() != outputVector[i_output].size()){

      if(_printWarnings)
	std::cout << "WARNING! Two collections don't have same number of samples:\t"
		  << outputVector[i_output].size() << " " << od.size() << std::endl;
      
      //if the samples is larger, we make a new vector of the right size, trimmed down appropriately
      if(od.size() > outputVector[i_output].size()){
	std::vector<short> const& mc_trimmed = std::vector<short>(od.begin(),
								  od.begin()+outputVector[i_output].size());
	fRDAdderAlg.AddRawDigits(mc_trimmed,outputVector[i_output]);
      }
      //if the samples is shorter, pad it out with the pedestal
      else if(od.size() < outputVector[i_output].size()){
	std::vector<short> mc_trimmed(outputVector[i_output].size(),0.0);
	std::copy(od.begin(),od.end(),mc_trimmed.begin());
	fRDAdderAlg.AddRawDigits(mc_trimmed,outputVector[i_output]);
      }
    }
    //Sizes are the same? Easy!
    else{
      fRDAdderAlg.AddRawDigits(od,outputVector[i_output]);
    }
  }

}

#endif
