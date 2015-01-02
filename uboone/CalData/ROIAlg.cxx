/*!
 * Title:   ROIAlg base class
 * Author:  wketchum@lanl.gov
 * Inputs:  
 * Outputs: 
 *
 * Description:
 * This is a base class for ROIAlgs.
 */

#include "ROIAlg.h"
#include "ROIAlg_DigitAboveThreshold.h"

template <class Digit>
std::unique_ptr< util::ROIAlg<Digit> > util::ROIAlg<Digit>::MakeROIAlg(fhicl::ParameterSet const& p){

  fAlgName = p.get<std::string>("AlgName");

  std::unique_ptr< ROIAlg<Digit> > ptr;
  if(fAlgName.compare("DigitAboveThreshold")){
    std::unique_ptr< ROIAlg<Digit> > new_ptr(new ROIAlg_DigitAboveThreshold<Digit>(p));
    ptr.swap(new_ptr);
  }
  else
    throw std::runtime_error("ERROR in ROIAlg: No registered ROIAlg with that name.");

  ptr->ClearRangeSets();

  return std::move(ptr);
}

template <class Digit>
void util::ROIAlg<Digit>::ClearRangeSets(){
  fSignalRangeSet.clear();
  fBaselineRangeSet.clear();
}

template <class Digit>
void util::ROIAlg<Digit>::ProcessWaveform(Region const& r){
  ClearRangeSets();
  fWaveformStart = r.Start();
  fWaveformEnd   = r.End();
  AnalyzeWaveform(r);
  CreateBaselineRangeSet();

  if(fBaselineRangeSet.size()!=(fSignalRangeSet.size()+1))
    throw std::runtime_error("ERROR in ROIAlg: BaselineRangeSet size is not exactly one greater than SignalRangeSet size.");
}


template <class Digit>
void util::ROIAlg<Digit>::GetSignalAndBaselineRegions( size_t const i_roi,
						       SignalBaselineTrio& trio){
  ThrowIfNoBaselineRegions();

  if( (i_roi+1) > fSignalRangeSet.size() )
    throw std::runtime_error("ERROR in ROIAlg: asked for roi iter that's larger than SignalRangeSet size");

  typename UniqueRangeSet<Tick>::iterator iter_signal = fSignalRangeSet.begin();
  typename UniqueRangeSet<Tick>::iterator iter_baseline_prev = fBaselineRangeSet.begin();
  typename UniqueRangeSet<Tick>::iterator iter_baseline_next = fBaselineRangeSet.begin();

  std::advance(iter_signal,i_roi);
  std::advance(iter_baseline_prev,i_roi);
  std::advance(iter_baseline_next,i_roi+1);

  trio.SignalRegion = *iter_signal;
  trio.BaselineRegion_Pre = *iter_baseline_prev;
  trio.BaselineRegion_Post = *iter_baseline_next;
  
}

//first tuple element is signal, then b_prev, then b_next
template <class Digit>
void util::ROIAlg<Digit>::GetAllSignalAndBaselineRegions( std::vector< SignalBaselineTrio >& trio_vector){

  ThrowIfNoBaselineRegions();

  trio_vector.clear(); trio_vector.reserve( fSignalRangeSet.size() );
  typename UniqueRangeSet<Tick>::iterator iter_signal = fSignalRangeSet.begin();
  typename UniqueRangeSet<Tick>::iterator iter_baseline_prev = fBaselineRangeSet.begin();
  typename UniqueRangeSet<Tick>::iterator iter_baseline_next = fBaselineRangeSet.begin();

  for(size_t i_signal=0; i_signal < fSignalRangeSet.size(); i_signal++){
    iter_baseline_next++;
    trio_vector[i_signal].SignalRegion = *iter_signal;
    trio_vector[i_signal].BaselineRegion_Pre = *iter_baseline_prev;
    trio_vector[i_signal].BaselineRegion_Post = *iter_baseline_next;
    iter_signal++;
    iter_baseline_prev++;
  }
  
  
}

template <class Digit>
const util::UniqueRangeSet<typename util::ROIAlg<Digit>::Tick> util::ROIAlg<Digit>::GetSignalRegions(){ 
  ThrowIfNoBaselineRegions();
  return fSignalRangeSet; 
}

template <class Digit>
const size_t util::ROIAlg<Digit>::GetNSignalRegions(){ 
  ThrowIfNoBaselineRegions();
  return fSignalRangeSet.size(); 
}

template <class Digit>
const util::UniqueRangeSet<typename util::ROIAlg<Digit>::Tick> util::ROIAlg<Digit>::GetBaselineRegions(){ 
  ThrowIfNoBaselineRegions();
  return fBaselineRangeSet;
}

template <class Digit>
const size_t util::ROIAlg<Digit>::GetNBaselineRegions(){ 
  ThrowIfNoBaselineRegions();
  return fBaselineRangeSet.size(); 
}

template <class Digit>
void util::ROIAlg<Digit>::ThrowIfNoBaselineRegions(){
  if(fBaselineRangeSet.size()==0)
    throw std::runtime_error("ERROR in ROIAlg: ProcessWaveform not yet called, but trying to access results!");
}

template <class Digit>
void util::ROIAlg<Digit>::CreateBaselineRangeSet(){
  fBaselineRangeSet.clear();
  Tick start,end;

  start = fWaveformStart; //initialize start
  //go through signal ranges
  //signal begin = baseline end, and baseline end = signal begin
  for(auto const& range : fSignalRangeSet){
    end = range.Start();
    fBaselineRangeSet.emplace(start,end);
    start = range.End();
  }
  //now need to do last baseline region: last signal end to waveform end
  end = fWaveformEnd;
  fBaselineRangeSet.emplace(start,end);
  
}
