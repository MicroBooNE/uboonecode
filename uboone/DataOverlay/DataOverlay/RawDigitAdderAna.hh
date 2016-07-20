/**
 * \file RawDigitAdderAna.h
 *
 * 
 * \brief Little sample program for establishing a user analysis space.
 *
 * @author wketchum
*/

#ifndef MIX_RAWDIGITADDERANA_H
#define MIX_RAWDIGITADDERANA_H

#include<vector>
#include<string>

#include "lardataobj/RawData/RawDigit.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

#include "TH1S.h"


namespace mix{
  class RawDigitAdderAna;
}

class mix::RawDigitAdderAna{
  
public:
  
  /// Default constructor
  RawDigitAdderAna(size_t sample=100,
		   std::vector<raw::ChannelID_t> special_channels=std::vector<raw::ChannelID_t>(0),
		   bool print_bad=true,
		   std::string in1label="in1",
		   std::string in2label="in2",
		   std::string sumlabel="sum");

  /// Default destructor
  virtual ~RawDigitAdderAna(){};

  //returns number of histograms to put in output file
  size_t CheckOverlay(std::vector<raw::RawDigit> const& in1,
		      std::vector<raw::RawDigit> const& in2,
		      std::vector<raw::RawDigit> const& sum);
  void CreateOutputHistograms(std::vector<TH1S*> const histoPtrVector,
			      std::vector<raw::RawDigit> const& in1,
			      std::vector<raw::RawDigit> const& in2,
			      std::vector<raw::RawDigit> const& sum,
			      unsigned int run, unsigned int event);
  
 private:

  size_t                        fChannelSampleInterval;
  std::vector<raw::ChannelID_t> fChannelsSpecial;  
  bool                          fPrintBadOverlays;
  std::string                   fInput1Label;
  std::string                   fInput2Label;
  std::string                   fSumLabel;
  std::vector<raw::ChannelID_t> fChannelsToPrint;
  std::vector<raw::ChannelID_t> fChannelsBadOverlay;

  void ResetOutput();
  void CreateHistogram(TH1S* histo,
		       raw::RawDigit const& waveform,
		       unsigned int run,
		       unsigned int event,
		       unsigned int channel,
		       std::string label);
};

#endif
