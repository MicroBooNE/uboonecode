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

#include "ROIAlg.h"

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
    void AnalyzeWaveform(Region const&);

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
