/*!
 * Title:   ROIAlg_DigitAboveThreshold
 * Author:  wketchum@lanl.gov
 * Inputs:  
 * Outputs: 
 *
 * Description:
 * Very simple digit above threshold ROI finder
 */

#include "ROIAlg_DigitAboveThreshold.h"

template <class Digit>
void util::ROIAlg_DigitAboveThreshold<Digit>::AnalyzeWaveform(Region const& region){
  
  unsigned int above_threshold_counter=0;
  Tick start_tick,end_tick;
  
  for( Tick tick=region.Start(); tick!=region.End(); tick++){
    
    //first, if we were above threshold, and now we're below
    if( above_threshold_counter>=fMinWidth && !(PassesThreshold(*tick)) ){
      end_tick = tick;
      InsertSignalRegion(start_tick,end_tick);
      above_threshold_counter=0;
      continue;
    }
    
    if( PassesThreshold(*tick) ){
      if(above_threshold_counter==0) start_tick = tick;
      above_threshold_counter++;
    }
    
    
  }//end loop over waveform
  
  //handle case where we end in signal region
  if(above_threshold_counter>fMinWidth)
    InsertSignalRegion(start_tick,region.End());
  
}
