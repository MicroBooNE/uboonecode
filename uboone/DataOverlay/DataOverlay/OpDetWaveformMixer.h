/**
 * \file OpDetWaveformMixer.h
 *
 * \ingroup DataOverlay
 * 
 * \brief Mixer function for putting together two raw digit collections
 *
 * @author wketchum
 */

/** \addtogroup DataOverlay

    @{*/
#ifndef OVERLAY_DATAOVERLAY_OPDETWAVEFORMMIXER_H
#define OVERLAY_DATAOVERLAY_OPDETWAVEFORMMIXER_H

#include <vector>
#include <string>
#include <unordered_map>

#include "lardataobj/RawData/OpDetWaveform.h"
#include "RawDigitAdder_HardSaturate.h"


/**
   \class OpDetWaveformMixer
   Add two raw digit collections together.
   
*/

namespace mix {
  class OpDetWaveformMixer;
}

class mix::OpDetWaveformMixer{

public:

  /// Default constructor
  OpDetWaveformMixer(bool p=false):
  _printWarnings(p){};

  void DeclareData(std::vector<raw::OpDetWaveform> const& dataVector,
		   std::vector<raw::OpDetWaveform> & outputVector);
  void Mix(std::vector<raw::OpDetWaveform> const& mcVector,
	   std::unordered_map<raw::Channel_t,float> const& map,
	   std::vector<raw::OpDetWaveform> & outputVector);
  
  void SetSaturationPoint(short x)
  { fRDAdderAlg.SetSaturationPoint(x); }
  
  void SetMinSampleSize(size_t x)
  { fMinSampleSize = x; }

  /// Default destructor
  virtual ~OpDetWaveformMixer(){};
  
  
 private:
  
  bool _printWarnings;

  std::unordered_map<raw::Channel_t,size_t> fChannelIndexMap;
  
  RawDigitAdder_HardSaturate fRDAdderAlg;

  size_t fMinSampleSize;
  
};

#endif
/** @} */ // end of doxygen group 

