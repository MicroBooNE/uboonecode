#ifndef OVERLAY_DATAOVERLAY_RAWDIGITMIXER_CXX
#define OVERLAY_DATAOVERLAY_RAWDIGITMIXER_CXX

#include "RawDigitMixer.h"
#include <limits>
#include <iostream>
#include <stdexcept>
#include <algorithm>


void mix::RawDigitMixer::DeclareData(std::vector<raw::RawDigit> const& dataVector){

  fChannelIndexMap.clear();
  fOutputWaveforms.clear();
  fOutputWaveforms.resize(dataVector.size());

  for(size_t i_rd=0; i_rd<dataVector.size(); i_rd++){

    if(dataVector[i_rd].Compression()!=raw::Compress_t::kNone)
      throw std::runtime_error("Error in RawDigitMixer::DeclareData : Compressed waveforms not yet supported!");
    
    fChannelIndexMap[dataVector[i_rd].Channel()] = i_rd;

    //initialize adc vector for output
    fOutputWaveforms[i_rd].waveform.resize(dataVector[i_rd].Samples(),0);
    fOutputWaveforms[i_rd].channel = dataVector[i_rd].Channel();
    fOutputWaveforms[i_rd].ped     = dataVector[i_rd].GetPedestal();
    fOutputWaveforms[i_rd].sigma   = dataVector[i_rd].GetSigma();
    
    
    //do the steps for filling that output vector
    //This is data, so set scales to one.
    fRDAdderAlg.SetScaleInputs(1.0,1.0);
    fRDAdderAlg.SetPedestalInputs(0.0,0.0);
    fRDAdderAlg.AddRawDigits(dataVector[i_rd].ADCs(),fOutputWaveforms[i_rd].waveform);
  }
  
}

void mix::RawDigitMixer::Mix(std::vector<raw::RawDigit> const& mcVector,
			     std::unordered_map<raw::ChannelID_t,float> const& scale_map){


  for( auto const& rd : mcVector){

    if(rd.Compression()!=raw::Compress_t::kNone)
      throw std::runtime_error("Error in RawDigitMixer::Mix : Compressed waveforms not yet supported!");
	
    auto it_ch = fChannelIndexMap.find(rd.Channel());

    //if this channel is not in the data, skip this channel!
    if(it_ch==fChannelIndexMap.end())
      continue;

    size_t i_output = it_ch->second;

    fRDAdderAlg.SetPedestalInputs(rd.GetPedestal(),0.0);
    fRDAdderAlg.SetScaleInputs(scale_map.at(rd.Channel()),1.0);
    
    //If the sizes are not the same...
    if(rd.Samples() != fOutputWaveforms[i_output].waveform.size()){
      if(_printWarnings)
	std::cout << "WARNING! Two collections don't have same number of samples:\t"
		  << fOutputWaveforms[i_output].waveform.size() << " " << rd.Samples() << std::endl;

      //if the samples is larger, we make a new vector of the right size, trimmed down appropriately
      if(rd.Samples() > fOutputWaveforms[i_output].waveform.size()){
	std::vector<short> const& mc_trimmed = std::vector<short>(rd.ADCs().begin(),
								  rd.ADCs().begin()+fOutputWaveforms[i_output].waveform.size());
	fRDAdderAlg.AddRawDigits(mc_trimmed,fOutputWaveforms[i_output].waveform);
      }
      //if the samples is shorter, pad it out with the pedestal
      else if(rd.Samples() < fOutputWaveforms[i_output].waveform.size()){
	std::vector<short> mc_trimmed(fOutputWaveforms[i_output].waveform.size(),rd.GetPedestal());
	std::copy(rd.ADCs().begin(),rd.ADCs().end(),mc_trimmed.begin());
	fRDAdderAlg.AddRawDigits(mc_trimmed,fOutputWaveforms[i_output].waveform);
      }
    }
    //Sizes are the same? Easy!
    else{
      fRDAdderAlg.AddRawDigits(rd.ADCs(),fOutputWaveforms[i_output].waveform);
    }

  }

}

void mix::RawDigitMixer::FillRawDigitOutput(std::vector<raw::RawDigit> & output){

  for(auto const& rd_i : fOutputWaveforms){
    
    //now emplace back onto output collection...
    output.emplace_back(rd_i.channel,
			rd_i.waveform.size(),
			rd_i.waveform);
    
    //set pedestal and rms to be same as data
    output.back().SetPedestal(rd_i.ped,rd_i.sigma);
  }
  
}


#endif
