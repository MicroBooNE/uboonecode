/*!
 * Title:   WaveformPropertiesAlg
 * Author:  wketchum@lanl.gov
 * Inputs:  Waveform (vector of digits, like short or float)
 * Outputs: various waveform properties
 *
 * Description:
 * This algorithm is intended to provide tools for investigating a waveform.
 */

#include "WaveformPropertiesAlg.h"
#include <algorithm>

template <class Digit>
bool util::WaveformPropertiesAlg<Digit>::RegionIsCurrent(Region const& region){
  return ( region.Start()==fCurrentProcessedRegion.Start() &&
	   region.End()==fCurrentProcessedRegion.End());
}

template <class Digit>
void util::WaveformPropertiesAlg<Digit>::ProcessWaveform(Region const& r){
  fCurrentProcessedRegion = r;
  fROIAlgPtr->ProcessWaveform(r);
}

template <class Digit>
Digit util::WaveformPropertiesAlg<Digit>::GetSum(Region const& r){
  Digit sum=0;
  for(Tick const& t=r.Start(); t!=r.End(); t++)
    sum += (float)*t;

  return sum;
}

template <class Digit>
float util::WaveformPropertiesAlg<Digit>::GetAverage(Region const& r){
  if(r.Start()==r.End()) return 0;
  float sum=(float)GetSum(r);
  return sum/std::distance(r.Start(),r.End());
}

template <class Digit>
float util::WaveformPropertiesAlg<Digit>::GetRMS(Region const& r){
  if(r.Start()==r.End()) return 0;
  float sum2=0;
  const float average = GetAverage(r);

  for(Tick const& t=r.Start(); t!=r.End(); t++)
    sum2 += ((float)*t - average)*((float)*t - average);

  return std::sqrt(sum2/std::distance(r.Start(),r.End()));
}

template <class Digit>
typename util::WaveformPropertiesAlg<Digit>::Tick util::WaveformPropertiesAlg<Digit>::GetMaxLocation(Region const& r){
  return std::max_element(r.Start(),r.End());
}

template <class Digit>
typename util::WaveformPropertiesAlg<Digit>::Tick util::WaveformPropertiesAlg<Digit>::GetMinLocation(Region const& r){
  return std::min_element(r.Start(),r.End());
}

template <class Digit>
float util::WaveformPropertiesAlg<Digit>::GetWaveformPedestal(Region const& r){
  if(!RegionIsCurrent(r))
    ProcessWaveform(r);

  float  sum = 0;
  size_t total_entries=0;
  for( auto const& range : fROIAlgPtr->GetBaselineRegions()){
    sum += GetSum(range);
    total_entries += std::distance(range.Start(),range.End());
  }

  if(total_entries!=0) sum = sum/total_entries;

  return sum;
}

template <class Digit>
float util::WaveformPropertiesAlg<Digit>::GetWaveformNoise(Region const& r){
  if(!RegionIsCurrent(r))
    ProcessWaveform(r);

  float  sum2 = 0;
  const float pedestal = GetWaveformPedestal(r);
  size_t total_entries=0;
  for( auto const& range : fROIAlgPtr->GetBaselineRegions()){

    total_entries += std::distance(range.Start(),range.End());
    for(Tick const& t=range.Start(); t!=range.End(); t++)
      sum2 += ((float)*t - pedestal)*((float)*t - pedestal);
    
  }

  if(total_entries!=0) sum2 = std::sqrt(sum2/total_entries);

  return sum2;

}

template <class Digit>
bool util::WaveformPropertiesAlg<Digit>::IsSignalRegion(Tick const& t, Region const& r){
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

template <class Digit>
int util::WaveformPropertiesAlg<Digit>::GetSignalRegionNumber(Tick const& t, Region const& r){
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

template <class Digit>
typename util::WaveformPropertiesAlg<Digit>::Region util::WaveformPropertiesAlg<Digit>::GetRegion(Tick const& t, Region const& r){

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

template <class Digit>
typename util::ROIAlg<Digit>::SignalBaselineTrio util::WaveformPropertiesAlg<Digit>::GetSignalBaselineTrio(Tick const& t, Region const& r){

  if(t<r.Start() || t>r.End())
    throw std::runtime_error("Error in WaveformPropertiesAlg: Tick outside given range.");

  int i_roi = GetSignalRegionNumber(t,r);
  if(i_roi<0)
    throw std::runtime_error("Error in WaveformPropertiesAlg: Tick not in signal region.");

  typename util::ROIAlg<Digit>::SignalBaselineTrio sbt;
  return fROIAlgPtr->GetSignalAndBaselineRegions(i_roi,sbt);
}

template <class Digit>
float util::WaveformPropertiesAlg<Digit>::GetLocalPedestal(Tick const& t, Region const& r){

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
    typename ROIAlg<Digit>::SignalBaselineTrio sbt = GetSignalBaselineTrio(t,r);
    float sum1 = GetSum(sbt.BaselineRegion_Pre);
    size_t n1 = std::distance(sbt.BaselineRegion_Pre.Start(),sbt.BaselineRegion_Pre.End());
    float sum2 = GetSum(sbt.BaselineRegion_Post);
    size_t n2 = std::distance(sbt.BaselineRegion_Post.Start(),sbt.BaselineRegion_Post.End());
    return (sum1+sum2)/(n1+n2);
  }

}

template <class Digit>
float util::WaveformPropertiesAlg<Digit>::GetLocalNoise(Tick const& t, Region const& r){

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
    float sum1 = GetSum(sbt.BaselineRegion_Pre);
    size_t n1 = std::distance(sbt.BaselineRegion_Pre.Start(),sbt.BaselineRegion_Pre.End());
    float sum2 = GetSum(sbt.BaselineRegion_Post);
    size_t n2 = std::distance(sbt.BaselineRegion_Post.Start(),sbt.BaselineRegion_Post.End());
    float pedestal = (sum1+sum2)/(n1+n2);

    float sum_sq=0;
    for(Tick const& t=sbt.BaselineRegion_Pre.Start(); t!=sbt.BaselineRegion_Pre.End(); t++)
      sum_sq += ((float)*t - pedestal)*((float)*t - pedestal);
    for(Tick const& t=sbt.BaselineRegion_Post.Start(); t!=sbt.BaselineRegion_Post.End(); t++)
      sum_sq += ((float)*t - pedestal)*((float)*t - pedestal);
    
    return std::sqrt(sum_sq/(n1+n2));
    
  }

}

template <class Digit>
typename util::UniqueRangeSet< typename util::WaveformPropertiesAlg<Digit>::Tick > const&
util::WaveformPropertiesAlg<Digit>::GetSignalRegions(Region const& r){
  if(!RegionIsCurrent(r))
    ProcessWaveform(r);

  return fROIAlgPtr->GetSignalRegions();
}

template <class Digit>
typename util::UniqueRangeSet< typename util::WaveformPropertiesAlg<Digit>::Tick > const&
util::WaveformPropertiesAlg<Digit>::GetBaselineRegions(Region const& r){
  if(!RegionIsCurrent(r))
    ProcessWaveform(r);

  return fROIAlgPtr->GetBaselineRegions();
}
