#ifndef WAVEFORMPROPERTIESALG_H
#define WAVEFORMPROPERTIESALG_H

/*!
 * Title:   WaveformPropertiesAlg
 * Author:  wketchum@lanl.gov
 * Inputs:  Waveform (vector of digits, like short or float)
 * Outputs: various waveform properties
 *
 * Description:
 * This algorithm is intended to provide tools for investigating a waveform.
 */

#include <vector>
#include <string>
#include <exception>

#include "fhiclcpp/ParameterSet.h"
#include "ROIAlg.h"

namespace util{

  template <class Digit>
  class WaveformPropertiesAlg{

  public:
    
    typedef std::vector<Digit> Waveform;
    typedef typename Waveform::const_iterator Tick;
    typedef Range<Tick> Region;
    
    WaveformPropertiesAlg(fhicl::ParameterSet const& p):
      fCurrentProcessedRegion(Tick(),Tick())
      { 
	fhicl::ParameterSet pset_ROIAlg;
	if( p.get_if_present<fhicl::ParameterSet>("ROIAlgParams",pset_ROIAlg)){
	  std::unique_ptr< ROIAlg<Digit> > new_ptr(ROIAlg<Digit>::MakeROIAlg(pset_ROIAlg));
	  fROIAlgPtr.swap(new_ptr);
	}
      }
    
    //these methods work regardless of whether waveform has been processed or not
    Digit GetSum(Region const&);
    Digit GetSum(Waveform const& waveform) { return GetSum(Region(waveform.cbegin(),waveform.cend())); }

    float GetAverage(Waveform const& waveform) { return GetAverage(Region(waveform.cbegin(),waveform.cend())); }
    float GetAverage(Region const&);

    float GetRMS(Waveform const& waveform) { return GetRMS(Region(waveform.cbegin(),waveform.cend())); }
    float GetRMS(Region const&);

    Digit GetMax(Waveform const& waveform) { return GetMax(Region(waveform.cbegin(),waveform.cend())); }
    Digit GetMax(Region const& r) { return *(GetMaxLocation(r)); }
    Digit GetMin(Waveform const& waveform) { return GetMin(Region(waveform.cbegin(),waveform.cend())); }
    Digit GetMin(Region const& r) { return *(GetMinLocation(r)); } 
    
    Tick GetMaxLocation(Waveform const& waveform) { return GetMaxLocation(Region(waveform.cbegin(),waveform.cend())); }
    Tick GetMaxLocation(Region const&);

    Tick GetMinLocation(Waveform const& waveform) { return GetMinLocation(Region(waveform.cbegin(),waveform.cend())); }
    Tick GetMinLocation(Region const&);

    void ProcessWaveform(Waveform const& w) { ProcessWaveform(Region(w.cbegin(),w.cend())); }
    void ProcessWaveform(Region const&);

    //these return pedestal and noise from all baseline regions
    //requires that waveform be processed
    float GetWaveformPedestal(Waveform const& w) { return GetWaveformPedestal(Region(w.cbegin(),w.cend())); }
    float GetWaveformPedestal(Region const&);
    float GetWaveformNoise(Waveform const& w) { return GetWaveformNoise(Region(w.cbegin(),w.cend())); }
    float GetWaveformNoise(Region const&);

    //these return pedestal from region based on iterator
    //if iterator in baseline region, then return noise/pedestal for that region
    //if iterator in signal region, return average noise/pedestal in surrounding baseline regions
    float GetLocalPedestal(Waveform const& w, Tick const& t) { return GetLocalPedestal(t,Region(w.cbegin(),w.cend())); }
    float GetLocalPedestal(Tick const&, Region const&);
    float GetLocalNoise(Waveform const& w, Tick const& t) { return GetLocalNoise(t,Region(w.cbegin(),w.cend())); }
    float GetLocalNoise(Tick const&, Region const&);

    //true if Tick in signal region, false if not
    bool IsSignalRegion(Waveform const& w, Tick const& t) { return IsSignalRegion(t,Region(w.cbegin(),w.cend())); }
    bool IsSignalRegion(Tick const&, Region const&);

    //position of range in signal unique range set, -1 if not there
    int GetSignalRegionNumber(Waveform const& w, Tick const& t) { return GetSignalRegionNumber(t,Region(w.cbegin(),w.cend())); }
    int GetSignalRegionNumber(Tick const&, Region const&);

    //return range for given tick
    Region GetRegion(Waveform const& w, Tick const& t) { return GetRegion(t,Region(w.cbegin(),w.cend())); }
    Region GetRegion(Tick const&, Region const&);
    typename ROIAlg<Digit>::SignalBaselineTrio GetSignalBaselineTrio(Waveform const& w, Tick const& t) { return GetRegion(t,Region(w.cbegin(),w.cend())); }
    typename ROIAlg<Digit>::SignalBaselineTrio GetSignalBaselineTrio(Tick const&, Region const&);

    UniqueRangeSet<Tick> const& GetSignalRegions(Waveform const& w) { return GetSignalRegions(Region(w.cbegin(),w.cend())); }
    UniqueRangeSet<Tick> const& GetSignalRegions(Region const&);
    UniqueRangeSet<Tick> const& GetBaselineRegions(Waveform const& w) { return GetBaselineRegions(Region(w.cbegin(),w.cend())); }
    UniqueRangeSet<Tick> const& GetBaselineRegions(Region const&);

    const size_t GetNSignalRegions(Waveform const& w) { return GetNSignalRegions(Region(w.cbegin(),w.cend())); }
    const size_t GetNSignalRegions(Region const&);
    const size_t GetNBaselineRegions(Waveform const& w) { return GetNSignalRegions(Region(w.cbegin(),w.cend())); }
    const size_t GetNBaselineRegions(Region const&);

  private:
    
    std::unique_ptr< ROIAlg<Digit> > fROIAlgPtr;

    Region fCurrentProcessedRegion;
    bool   RegionIsCurrent(Region const&);

  };

}//end namespace util


#endif
