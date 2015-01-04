#ifndef ROIALGDIGITABOVETHRESHOLD_H
#define ROIALGDIGITABOVETHRESHOLD_H

/*!
 * Title:   ROIAlg_DigitAboveThreshold
 * Author:  wketchum@lanl.gov
 * Inputs:  
 * Outputs: 
 *
 * Description:
 * Very simple digit above threshold ROI finder
 */

#include <vector>
#include <string>
#include <exception>

namespace util{

  template <class Digit>
  class ROIAlg_DigitAboveThreshold : public ROIAlg<Digit> {
    
  public:

    //this should be copied from ROIAlg.h
    typedef std::vector<Digit> Waveform;
    typedef typename Waveform::const_iterator  Tick;
    typedef Range<Tick> Region;

    ROIAlg_DigitAboveThreshold(fhicl::ParameterSet const& p){
      fThresholdVal  = p.get<Digit>("ThresholdVal");
      fNegativePulse = p.get<bool>("NegativePulse",false);
      fMinWidth      = p.get<unsigned int>("MinWidth");
      
      //can't have minwidth of zero: minimum is 1
      if(fMinWidth==0) fMinWidth=1;
    }
    ~ROIAlg_DigitAboveThreshold(){}
    
  protected: 
    void AnalyzeWaveform(Region const& region){
  
      unsigned int above_threshold_counter=0;
      Tick start_tick,end_tick;
      
      for( Tick tick=region.Start(); tick!=region.End(); tick++){

	//first, if we were above threshold, and now we're below
	if( above_threshold_counter>=fMinWidth && !(PassesThreshold(*tick)) ){
	  end_tick = tick;
	  this->InsertSignalRegion(start_tick,end_tick);
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
	this->InsertSignalRegion(start_tick,region.End());
      
    }


  private:
    Digit         fThresholdVal;
    bool          fNegativePulse;
    unsigned int  fMinWidth;

    bool PassesThreshold(Digit const& val){
      if(fNegativePulse)
	return (val<=fThresholdVal);
      else
	return (val>=fThresholdVal);
    }
    
  };

} //end namespace util

#endif
