#ifndef ROIALGDIGITABOVEBASELINE_H
#define ROIALGDIGITABOVEBASELINE_H

/*!
 * Title:   ROIAlg_DigitAboveBaseline
 * Author:  wketchum@lanl.gov
 * Inputs:  
 * Outputs: 
 *
 * Description:
 * Very simple digit above baseline ROI finder
 */

#include <vector>
#include <string>
#include <exception>

namespace util{

  template <class Digit>
  class ROIAlg_DigitAboveBaseline : public ROIAlg<Digit> {
    
  public:

    //this should be copied from ROIAlg.h
    typedef std::vector<Digit> Waveform;
    typedef typename Waveform::const_iterator  Tick;
    typedef Range<Tick> Region;

    ROIAlg_DigitAboveBaseline(fhicl::ParameterSet const& p){
      fPaddingSize = p.get<unsigned int>("PaddingSize");
      fThresholdInRMS = p.get<double>("ThresholdInRMS");
    }
    ~ROIAlg_DigitAboveBaseline(){}
    
  protected: 
    void AnalyzeWaveform(Region const& region){

      //average of the region
      if(region.Start()==region.End()) return;
      const Digit baseline_threshold = DetermineBaselineThreshold(region);
      
      Tick start_tick,end_tick;
      
      for( Tick tick=region.Start(); tick!=region.End(); tick++){

	//first, if we were above threshold, and now we're below
	if( *tick < baseline_threshold ){
	  if( std::distance(tick,region.End()) < fPaddingSize)
	    end_tick = region.End();
	  else
	    end_tick = tick + fPaddingSize;

	  this->InsertSignalRegion(start_tick,end_tick);
	  continue;
	}
	
	if( *tick >= baseline_threshold ){
	  if( std::distance(region.Start(),tick) < fPaddingSize)
	    start_tick = region.Start();
	  else
	    start_tick = tick - fPaddingSize;
	}
	
      }//end loop over waveform
      
      //handle case where we end in signal region
      if(start_tick > end_tick)
	this->InsertSignalRegion(start_tick,region.End());
      
    }


  private:
    unsigned int  fPaddingSize;
    double        fThresholdInRMS;

    Digit DetermineBaselineThreshold(Region const& r)
    {
      double sum=0;
      for(Tick t=r.Start(); t!=r.End(); t++)
	sum += (double)*t;
      double average = sum/std::distance(r.Start(),r.End());

      double sum2=0;
      for(Tick t=r.Start(); t!=r.End(); t++)
	sum2 += ((double)*t - average)*((double)*t - average);
      
      double rms =  std::sqrt(sum2/std::distance(r.Start(),r.End()));

      return average + fThresholdInRMS*rms;
      
    }

  };

} //end namespace util

#endif
