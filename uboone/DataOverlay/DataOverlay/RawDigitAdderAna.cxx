#ifndef MIX_RAWDIGITADDERANA_CXX
#define MIX_RAWDIGITADDERANA_CXX

#include "RawDigitAdderAna.hh"

#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <sstream>

mix::RawDigitAdderAna::RawDigitAdderAna(size_t sample,
					std::vector<raw::ChannelID_t> special_channels,
					bool print_bad,
					std::string in1label,
					std::string in2label,
					std::string sumlabel)
  : fChannelSampleInterval(sample),
    fChannelsSpecial(special_channels),
    fPrintBadOverlays(print_bad),
    fInput1Label(in1label),
    fInput2Label(in2label),
    fSumLabel(sumlabel)
{
}

void mix::RawDigitAdderAna::ResetOutput(){
  fChannelsToPrint.clear();
  fChannelsBadOverlay.clear();
}

size_t mix::RawDigitAdderAna::CheckOverlay(std::vector<raw::RawDigit> const& in1,
					   std::vector<raw::RawDigit> const& in2,
					   std::vector<raw::RawDigit> const& sum)
{

  ResetOutput();
  
  //check size of all the rawdigit vectors
  if(in1.size()!=in2.size() ||
     in1.size()!=sum.size() )
    throw std::runtime_error("Error in RawDigitAdderAna::CheckOverlay : Input vector lists not equal size.");

  for(size_t i_ch=0; i_ch<sum.size(); i_ch++){

    //check to make sure channels are same
    if(in1[i_ch].Channel()!=in2[i_ch].Channel() ||
       in1[i_ch].Channel()!=sum[i_ch].Channel())
      throw std::runtime_error("Error in RawDigitAdderAna::CheckOverlay : Input vector lists out of order.");
    
    //check to make sure channels have same ADC vector
    if(in1[i_ch].ADCs().size()!=in2[i_ch].ADCs().size() ||
       in1[i_ch].ADCs().size()!=sum[i_ch].ADCs().size())
      throw std::runtime_error("Error in RawDigitAdderAna::CheckOverlay : Input channels have different sizes.");

    bool print= ( (in1[i_ch].Channel()%fChannelSampleInterval==0) ||
		  (std::find(fChannelsSpecial.begin(),fChannelsSpecial.end(),in1[i_ch].Channel())!=fChannelsSpecial.end()) );

    for(size_t i_t=0; i_t<sum[i_ch].ADCs().size(); i_t++){
      if( (in1[i_ch].ADCs()[i_t]+in2[i_ch].ADCs()[i_t])==sum[i_ch].ADCs()[i_t] )
	continue;

      fChannelsBadOverlay.emplace_back(in1[i_ch].Channel());
      if(fPrintBadOverlays) print=true;
      break;
    }

    if(print) fChannelsToPrint.emplace_back(in1[i_ch].Channel());
    
  }

  return fChannelsToPrint.size()*3;
  
}

void mix::RawDigitAdderAna::CreateHistogram(TH1S* histo,
					    raw::RawDigit const& waveform,
					    unsigned int run,
					    unsigned int event,
					    unsigned int channel,
					    std::string label)
{
    std::stringstream hname,htitle;
    hname << "h_" << label
	  << "_run" << run
	  << "_ev" << event
	  << "_ch" << channel;
    htitle << "Waveform, Channel " << channel
	   << ", Event " << event
	   << ", Run " << run
	   << ", Input " << label;
    histo->SetName(hname.str().c_str());
    histo->SetTitle(htitle.str().c_str());
    histo->SetBins(waveform.ADCs().size(),0,waveform.ADCs().size());

    for(size_t i_t=0; i_t<waveform.ADCs().size(); i_t++)
      histo->SetBinContent(i_t,waveform.ADCs()[i_t]);
    
}
					    

void mix::RawDigitAdderAna::CreateOutputHistograms(std::vector<TH1S*> const histoPtrVector,
						   std::vector<raw::RawDigit> const& in1,
						   std::vector<raw::RawDigit> const& in2,
						   std::vector<raw::RawDigit> const& sum,
						   unsigned int run, unsigned int event)
{

  if(histoPtrVector.size()!=fChannelsToPrint.size())
    throw std::runtime_error("Error in RawDigitAdderAna::CreateOutputHistograms : Histogram output vector not equal to channels to print.");
  

  //check size of all the rawdigit vectors
  if(in1.size()!=in2.size() ||
     in1.size()!=sum.size() )
    throw std::runtime_error("Error in RawDigitAdderAna::CheckOverlay : Input vector lists not equal size.");

  size_t histo_count=0;
  
  for(size_t i_ch=0; i_ch<sum.size(); i_ch++){
    
    //check to make sure channels are same
    if(in1[i_ch].Channel()!=in2[i_ch].Channel() ||
       in1[i_ch].Channel()!=sum[i_ch].Channel())
      throw std::runtime_error("Error in RawDigitAdderAna::CheckOverlay : Input vector lists out of order.");

    if(std::find(fChannelsToPrint.begin(),fChannelsToPrint.end(),in1[i_ch].Channel())==fChannelsToPrint.end())
      continue;
    
    //check to make sure channels have same ADC vector
    if(in1[i_ch].ADCs().size()!=in2[i_ch].ADCs().size() ||
       in1[i_ch].ADCs().size()!=sum[i_ch].ADCs().size())
      throw std::runtime_error("Error in RawDigitAdderAna::CheckOverlay : Input channels have different sizes.");

    CreateHistogram(histoPtrVector.at(histo_count),
		    in1[i_ch],
		    run,event,in1[i_ch].Channel(),
		    fInput1Label);
    histo_count++;
    CreateHistogram(histoPtrVector.at(histo_count),
		    in2[i_ch],
		    run,event,in2[i_ch].Channel(),
		    fInput2Label);
    histo_count++;
    CreateHistogram(histoPtrVector.at(histo_count),
		    sum[i_ch],
		    run,event,sum[i_ch].Channel(),
		    fSumLabel);
    histo_count++;
    
  }

}

#endif
