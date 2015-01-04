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
#include <algorithm>

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
    double GetSum(Waveform const& waveform) { return GetSum(Region(waveform.cbegin(),waveform.cend())); }
    double GetSum(Region const& r){
      double sum=0;
      for(Tick t=r.Start(); t!=r.End(); t++)
	sum += (double)*t;
      
      return sum;
    }

    double GetAverage(Waveform const& waveform) { return GetAverage(Region(waveform.cbegin(),waveform.cend())); }
    double GetAverage(Region const& r){
      if(r.Start()==r.End()) return 0;
      double sum=GetSum(r);
      return sum/std::distance(r.Start(),r.End());
    }
    
    double GetRMS(Waveform const& waveform) { return GetRMS(Region(waveform.cbegin(),waveform.cend())); }
    double GetRMS(Region const& r){
      if(r.Start()==r.End()) return 0;
      double sum2=0;
      const double average = GetAverage(r);
      
      for(Tick t=r.Start(); t!=r.End(); t++)
	sum2 += ((double)*t - average)*((double)*t - average);
      
      return std::sqrt(sum2/std::distance(r.Start(),r.End()));
    }


    Digit GetMax(Waveform const& waveform) { return GetMax(Region(waveform.cbegin(),waveform.cend())); }
    Digit GetMax(Region const& r) { return *(GetMaxLocation(r)); }
    Digit GetMin(Waveform const& waveform) { return GetMin(Region(waveform.cbegin(),waveform.cend())); }
    Digit GetMin(Region const& r) { return *(GetMinLocation(r)); } 
    
    Tick GetMaxLocation(Waveform const& waveform) { return GetMaxLocation(Region(waveform.cbegin(),waveform.cend())); }
    Tick GetMaxLocation(Region const& r) { return std::max_element(r.Start(),r.End()); }

    Tick GetMinLocation(Waveform const& waveform) { return GetMinLocation(Region(waveform.cbegin(),waveform.cend())); }
    Tick GetMinLocation(Region const& r) { return std::min_element(r.Start(),r.End()); }

    void ProcessWaveform(Waveform const& w) { ProcessWaveform(Region(w.cbegin(),w.cend())); }
    void ProcessWaveform(Region const& r){

      if(!fROIAlgPtr)
	throw std::runtime_error("ERROR in WaveformPropertiesAlg: ROIAlg called but no params given.");
      
      fCurrentProcessedRegion = r;
      fROIAlgPtr->ProcessWaveform(r);
    }

    //these return pedestal and noise from all baseline regions
    //requires that waveform be processed
    double GetWaveformPedestal(Waveform const& w) { return GetWaveformPedestal(Region(w.cbegin(),w.cend())); }
    double GetWaveformPedestal(Region const& r){
      if(!RegionIsCurrent(r))
	ProcessWaveform(r);
      
      double  sum = 0;
      size_t total_entries=0;
      for( auto const& range : fROIAlgPtr->GetBaselineRegions()){
	sum += GetSum(range);
	total_entries += std::distance(range.Start(),range.End());
      }
      
      if(total_entries!=0) sum = sum/total_entries;
      
      return sum;
    }
    
    double GetWaveformNoise(Waveform const& w) { return GetWaveformNoise(Region(w.cbegin(),w.cend())); }
    double GetWaveformNoise(Region const& r){
      if(!RegionIsCurrent(r))
	ProcessWaveform(r);
      
      double  sum2 = 0;
      const double pedestal = GetWaveformPedestal(r);
      size_t total_entries=0;
      for( auto const& range : fROIAlgPtr->GetBaselineRegions()){
	
	total_entries += std::distance(range.Start(),range.End());
	for(Tick t=range.Start(); t!=range.End(); t++)
	  sum2 += ((double)*t - pedestal)*((double)*t - pedestal);
	
      }
      
      if(total_entries!=0) sum2 = std::sqrt(sum2/total_entries);
      
      return sum2;
      
    }


    //these return pedestal from region based on iterator
    //if iterator in baseline region, then return noise/pedestal for that region
    //if iterator in signal region, return average noise/pedestal in surrounding baseline regions
    double GetLocalPedestal(Waveform const& w, Tick const& t) { return GetLocalPedestal(t,Region(w.cbegin(),w.cend())); }
    double GetLocalPedestal(Tick const& t, Region const& r){

      if(t<r.Start() || t>r.End())
	throw std::runtime_error("Error in WaveformPropertiesAlg: Tick outside given range.");
      
      if(!RegionIsCurrent(r))
	ProcessWaveform(r);
      
      //if tick is in baseline region
      if(!IsSignalRegion(t,r)){
	return GetAverage( GetRegion(t,r) );
      }
      //else it's a signal region!
      else{
	//std::cout << "Is signal region!" << std::endl;
	typename ROIAlg<Digit>::SignalBaselineTrio sbt = GetSignalBaselineTrio(t,r);
	double sum1 = GetSum(sbt.BaselineRegion_Pre);
	size_t n1 = std::distance(sbt.BaselineRegion_Pre.Start(),sbt.BaselineRegion_Pre.End());
	double sum2 = GetSum(sbt.BaselineRegion_Post);
	size_t n2 = std::distance(sbt.BaselineRegion_Post.Start(),sbt.BaselineRegion_Post.End());
	double average = (double)(sum1+sum2)/(n1+n2);
	//std::cout << sum1 << "/" << n1 << " " << sum2 << "/" << n2 << " " << average << std::endl;;
	return average;
      }
      
    }

    double GetLocalNoise(Waveform const& w, Tick const& t) { return GetLocalNoise(t,Region(w.cbegin(),w.cend())); }
    double GetLocalNoise(Tick const& t, Region const& r){

      if(t<r.Start() || t>r.End())
	throw std::runtime_error("Error in WaveformPropertiesAlg: Tick outside given range.");
      
      if(!RegionIsCurrent(r))
	ProcessWaveform(r);
      
      //if tick is in baseline region
      if(!IsSignalRegion(t,r)){
	return GetRMS( GetRegion(t,r) );
      }
      //else it's a signal region!
      else{
	typename ROIAlg<Digit>::SignalBaselineTrio sbt = GetSignalBaselineTrio(t,r);
	double sum1 = GetSum(sbt.BaselineRegion_Pre);
	size_t n1 = std::distance(sbt.BaselineRegion_Pre.Start(),sbt.BaselineRegion_Pre.End());
	double sum2 = GetSum(sbt.BaselineRegion_Post);
	size_t n2 = std::distance(sbt.BaselineRegion_Post.Start(),sbt.BaselineRegion_Post.End());
	double pedestal = (double)(sum1+sum2)/(n1+n2);
	
	double sum_sq=0;
	for(Tick t=sbt.BaselineRegion_Pre.Start(); t!=sbt.BaselineRegion_Pre.End(); t++)
	  sum_sq += ((double)*t - pedestal)*((double)*t - pedestal);
	for(Tick t=sbt.BaselineRegion_Post.Start(); t!=sbt.BaselineRegion_Post.End(); t++)
	  sum_sq += ((double)*t - pedestal)*((double)*t - pedestal);
	
	return std::sqrt(sum_sq/(n1+n2));
	
      }
      
    }


    //true if Tick in signal region, false if not
    bool IsSignalRegion(Waveform const& w, Tick const& t) { return IsSignalRegion(t,Region(w.cbegin(),w.cend())); }
    bool IsSignalRegion(Tick const& t, Region const& r){
      if(t<r.Start() || t>r.End())
	throw std::runtime_error("Error in WaveformPropertiesAlg: Tick outside given range.");
      
      if(!RegionIsCurrent(r))
	ProcessWaveform(r);
      
      typename UniqueRangeSet<Tick>::iterator iter;
      iter=fROIAlgPtr->GetSignalRegions().find(Region(t,t+1));
      if (iter!=(fROIAlgPtr->GetSignalRegions()).end()) return true;
      
      iter=fROIAlgPtr->GetBaselineRegions().find(Region(t,t+1));
      if (iter!=(fROIAlgPtr->GetBaselineRegions()).end()) return false;
      
      throw std::runtime_error("Error in WaveformPropertiesAlg: Tick not found in RangeSets.");
      
    }

    //position of range in signal unique range set, -1 if not there
    int GetSignalRegionNumber(Waveform const& w, Tick const& t) { return GetSignalRegionNumber(t,Region(w.cbegin(),w.cend())); }
    int GetSignalRegionNumber(Tick const& t, Region const& r){
      if(t<r.Start() || t>r.End())
	throw std::runtime_error("Error in WaveformPropertiesAlg: Tick outside given range.");
      
      if(!RegionIsCurrent(r))
	ProcessWaveform(r);
      
      typename UniqueRangeSet<Tick>::iterator iter;
      iter=fROIAlgPtr->GetSignalRegions().find(Region(t,t+1));
      if (iter!=(fROIAlgPtr->GetSignalRegions()).end()) 
	return std::distance(fROIAlgPtr->GetSignalRegions().begin(),iter);
      
      //else, not in signal region
      return -1;
    }

    //return range for given tick
    Region GetRegion(Waveform const& w, Tick const& t) { return GetRegion(t,Region(w.cbegin(),w.cend())); }
    Region GetRegion(Tick const& t, Region const& r){
      
      if(t<r.Start() || t>r.End())
	throw std::runtime_error("Error in WaveformPropertiesAlg: Tick outside given range.");
      
      if(!RegionIsCurrent(r))
	ProcessWaveform(r);
      
      typename UniqueRangeSet<Tick>::iterator iter;
      iter=fROIAlgPtr->GetSignalRegions().find(Region(t,t+1));
      if (iter!=(fROIAlgPtr->GetSignalRegions()).end()) return *iter;;
      
      iter=fROIAlgPtr->GetBaselineRegions().find(Region(t,t+1));
      if (iter!=(fROIAlgPtr->GetBaselineRegions()).end()) return *iter;
      
      throw std::runtime_error("Error in WaveformPropertiesAlg: Tick not found in RangeSets.");
      
    }

    typename ROIAlg<Digit>::SignalBaselineTrio GetSignalBaselineTrio(Waveform const& w, Tick const& t) { return GetSignalBaselineTrio(t,Region(w.cbegin(),w.cend())); }
    typename ROIAlg<Digit>::SignalBaselineTrio GetSignalBaselineTrio(Tick const& t, Region const& r){
      
      if(t<r.Start() || t>r.End())
	throw std::runtime_error("Error in WaveformPropertiesAlg: Tick outside given range.");
      
      int i_roi = GetSignalRegionNumber(t,r);
      if(i_roi<0)
	throw std::runtime_error("Error in WaveformPropertiesAlg: Tick not in signal region.");
      
      typename util::ROIAlg<Digit>::SignalBaselineTrio sbt;
      fROIAlgPtr->GetSignalAndBaselineRegions(i_roi,sbt);
      return sbt;
    }

    UniqueRangeSet<Tick> const& GetSignalRegions(Waveform const& w) { return GetSignalRegions(Region(w.cbegin(),w.cend())); }
    UniqueRangeSet<Tick> const& GetSignalRegions(Region const& r){
      if(!RegionIsCurrent(r))
	ProcessWaveform(r);
      
      return fROIAlgPtr->GetSignalRegions();
    }

    UniqueRangeSet<Tick> const& GetBaselineRegions(Waveform const& w) { return GetBaselineRegions(Region(w.cbegin(),w.cend())); }
    UniqueRangeSet<Tick> const& GetBaselineRegions(Region const& r){
      if(!RegionIsCurrent(r))
	ProcessWaveform(r);
      
      return fROIAlgPtr->GetBaselineRegions();
    }

    const size_t GetNSignalRegions(Waveform const& w) { return GetNSignalRegions(Region(w.cbegin(),w.cend())); }
    const size_t GetNSignalRegions(Region const& r){
      if(!RegionIsCurrent(r))
	ProcessWaveform(r);
      
      return fROIAlgPtr->GetNSignalRegions();
    }
    
    const size_t GetNBaselineRegions(Waveform const& w) { return GetNSignalRegions(Region(w.cbegin(),w.cend())); }
    const size_t GetNBaselineRegions(Region const& r){
      if(!RegionIsCurrent(r))
	ProcessWaveform(r);
      
      return fROIAlgPtr->GetNBaselineRegions();
    }

  private:
    
    std::unique_ptr< ROIAlg<Digit> > fROIAlgPtr;

    Region fCurrentProcessedRegion;
    bool   RegionIsCurrent(Region const& region){
      return ( region.Start()==fCurrentProcessedRegion.Start() &&
	       region.End()==fCurrentProcessedRegion.End());
    }

  };

}//end namespace util


#endif
