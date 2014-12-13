#ifndef WAVEFORMPROPERTIESALG_H
#define WAVEFORMPROPERTIESALG_H

/*!
 * Title:   WaveformPropertiesAlg
 * Author:  wketchum@lanl.gov
 * Inputs:  Waveform (from raw::RawDigit)
 * Outputs: various waveform properties
 *
 * Description:
 * This algorithm is intended to provide tools for investigating a waveform.
 */

#include <vector>
#include <string>
#include <exception>

#include "fhiclcpp/ParameterSet.h"
#include "RawData/RawDigit.h"
#include "ROIAlg.h"

namespace util{

  class WaveformPropertiesAlg{

    typedef short                    Digit;
    typedef std::vector<Digit>       Waveform;
    typedef Waveform::const_iterator Tick;    
    typedef std::pair<Tick,Tick>     Region;
    
  public:
    
    WaveformPropertiesAlg(fhicl::ParameterSet const& p)
      { fROIAlgPtr.swap(ROIAlg::MakeROIAlg(p.get<fhicl::ParameterSet>("ROIAlgParams"))); }
    
    float GetPedestal(raw::RawDigit const&,
		      size_t min_tick=0,
		      size_t max_tick=Waveform::max_size());
    float GetNoise(raw::RawDigit const&,
		   size_t min_tick=0,
		   size_t max_tick=Waveform::max_size());
    Digit GetMax(raw::RawDigit const&,
		 size_t min_tick=0,
		 size_t max_tick=Waveform::max_size());
    Digit GetMin(raw::RawDigit const&,
		 size_t min_tick=0,
		 size_t max_tick=Waveform::max_size());
    TickIter GetMaxLocation(raw::RawDigit const&,
			    size_t min_tick=0,
			    size_t max_tick=Waveform::max_size());
    TickIter GetMinLocation(raw::RawDigit const&,
			    size_t min_tick=0,
			    size_t max_tick=Waveform::max_size());

    
    std::vector<Region> GetSignalRegions(raw::RawDigit const&,
					 size_t min_tick=0,
					 size_t max_tick=Waveform::max_size());
    std::vector<Region> GetBaselineRegions(raw::RawDigit const&,
					   size_t min_tick=0,
					   size_t max_tick=Waveform::max_size());
    
  private:
    
    std::unique_ptr<ROIAlg> fROIAlgPtr;

    typedef struct {
      unsigned int channel;
      size_t min_tick;
      size_t max_tick;
      float pedestal;
      float noise;
      Digit max;
      Digit min;
      TickIter max_iterator;
      TickIter min_iterator;
      std::vector<Region> signal_regions;
      std::vector<Region> baseline_regions;
    } ChannelProperties_t;
    
    std::vector<ChannelProperties_t> fChannelMemory;
    
    void ClearChannelMemory();
    
  };

}//end namespace util


#endif
