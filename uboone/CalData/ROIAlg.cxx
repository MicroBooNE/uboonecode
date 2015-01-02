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

std::unique_ptr<util::ROIAlg> util::ROIAlg::MakeROIAlg(fhicl::ParameterSet const& p){

  fAlgName = p.get<std::string>("AlgName");

  std::unique_ptr<ROIAlg> ptr;
  if(fAlgName.compare("DigitAboveThreshold")){
    std::unique_ptr<ROIAlg> new_ptr(new ROIAlg_DigitAboveThreshold(p));
    ptr.swap(new_ptr);
  }
  else
    throw std::runtime_error("ERROR in ROIAlg: No registered ROIAlg with that name.");

  return std::move(ptr);
}

void util::ROIAlg::ProcessWaveform(Waveform const& waveform){
  fSignalRangeSet.clear();
  fWaveformStart = waveform.cbegin();
  fWaveformEnd   = waveform.cend();
  AnalyzeWaveform(waveform);
}


void util::ROIAlg::GetSignalAndBaselineRegions( size_t const i_roi,
						SignalBaselineTrio& trio){
  if( (i_roi+1) > fSignalRangeSet.size() )
    throw std::runtime_error("ERROR in ROIAlg: asked for roi iter that's larger than SignalRangeSet size");

  UniqueRangeSet<Tick>::iterator iter = fSignalRangeSet.begin();

  //set prev baseline region, also iter is incremented to signal region
  if(i_roi==0){
    trio.BaselineRegion_Pre = Region(fWaveformStart,iter->Start());
  }
  else{
    std::advance(iter,i_roi-1);
    Tick prev_start = iter->End();
    std::advance(iter,1);
    trio.BaselineRegion_Pre = Region(prev_start,iter->Start());
  }

  //set signal region
  trio.SignalRegion = *iter;

  //now, set next baseline region
  if( i_roi==(fSignalRangeSet.size()-1) ){
    trio.BaselineRegion_Post = Region(iter->End(),fWaveformEnd);
  }
  else{
    Tick next_start = iter->End();
    std::advance(iter,1);
    trio.BaselineRegion_Post = Region(next_start,iter->Start());
  }
  
}

//first tuple element is signal, then b_prev, then b_next
void util::ROIAlg::GetAllSignalAndBaselineRegions( std::vector< SignalBaselineTrio >& trio_vector){

  trio_vector.clear(); trio_vector.reserve( fSignalRangeSet.size() );
  UniqueRangeSet<Tick> baselineRangeSet = CreateBaselineRangeSet();

  UniqueRangeSet<Tick>::iterator iter_signal = fSignalRangeSet.begin();
  UniqueRangeSet<Tick>::iterator iter_baseline_prev = baselineRangeSet.begin();
  UniqueRangeSet<Tick>::iterator iter_baseline_next = baselineRangeSet.begin();

  for(size_t i_signal=0; i_signal < fSignalRangeSet.size(); i_signal++){
    iter_baseline_next++;
    trio_vector[i_signal].SignalRegion = *iter_signal;
    trio_vector[i_signal].BaselineRegion_Pre = *iter_baseline_prev;
    trio_vector[i_signal].BaselineRegion_Post = *iter_baseline_next;
    iter_signal++;
    iter_baseline_prev++;
  }
  
  
}

const util::UniqueRangeSet<util::Tick> util::ROIAlg::GetSignalRegions() { return fSignalRangeSet; }

const size_t util::ROIAlg::GetNSignalRegions() { return fSignalRangeSet.size(); }

const util::UniqueRangeSet<util::Tick> util::ROIAlg::GetBaselineRegions() { return CreateBaselineRangeSet(); }

const size_t util::ROIAlg::GetNBaselineRegions() { return CreateBaselineRangeSet().size(); }

util::UniqueRangeSet<util::Tick> util::ROIAlg::CreateBaselineRangeSet() const{
  UniqueRangeSet<Tick> baselineRangeSet;
  Tick start,end;

  start = fWaveformStart; //initialize start
  //go through signal ranges
  //signal begin = baseline end, and baseline end = signal begin
  for(auto const& range : fSignalRangeSet){
    end = range.Start();
    baselineRangeSet.emplace(start,end);
    start = range.End();
  }
  //now need to do last baseline region: last signal end to waveform end
  end = fWaveformEnd;
  baselineRangeSet.emplace(start,end);
  
  return baselineRangeSet;
}
