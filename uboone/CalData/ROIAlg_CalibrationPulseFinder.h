#ifndef ROIALGCALIBRATIONPULSEFINDER_H
#define ROIALGCALIBRATIONPULSEFINDER_H

/*!
 * Title:   ROIAlg_CalibrationPulseFinder
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
  class ROIAlg_CalibrationPulseFinder : public ROIAlg<Digit> {
    
  public:

    //this should be copied from ROIAlg.h
    typedef std::vector<Digit> Waveform;
    typedef typename Waveform::const_iterator  Tick;
    typedef Range<Tick> Region;

    ROIAlg_CalibrationPulseFinder(fhicl::ParameterSet const& p){
      fPrePulseBins  = p.get<unsigned int>("PrePulseBins");
      fThreshold     = p.get<Digit>("Threshold");
      fMinWidth      = p.get<unsigned int>("MinWidth");
      
      //can't have minwidth of zero: minimum is 1
      if(fMinWidth==0) fMinWidth=1;
    }
    ~ROIAlg_CalibrationPulseFinder(){}
    
  protected: 
    void AnalyzeWaveform(Region const& region){
      
      // If size of region is less than fPrePulseBins
      // error! send out a warning
      if (std::distance(region.Start(),region.End()) < fPrePulseBins)
	throw std::runtime_error("ERROR in ROIAlg_CalibrationPulseFinder: size of Baseline region requested larger than size of region being analyzed.");

      unsigned int above_threshold_counter=0;
      Tick start_tick,end_tick;

      unsigned int tickCounter = 0;

      fBaseline = 0;
 
      // Set iterator to beginning of region
      Tick tick = region.Start();
      // First find the baseline
      while (tickCounter < fPrePulseBins){
	fBaseline += *tick;
	tickCounter += 1;
	tick++;
      }
      fBaseline /= fPrePulseBins;
      // now continue with the rest of the waveform
      while (tick != region.End() ){
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
	tick++;
      }// end loop over waveform
      
      //handle case where we end in signal region
      if(above_threshold_counter>fMinWidth)
	this->InsertSignalRegion(start_tick,region.End());
      
    }

    
  private:
    unsigned int  fPrePulseBins;
    Digit         fThreshold;
    unsigned int  fMinWidth;
    double        fBaseline;

    bool PassesThreshold(Digit const& val){
      if ( ((val-fBaseline) > fThreshold) || ((fBaseline-val) > fThreshold) )
	return true;
      return false;
    }
    
  };

} //end namespace util

#endif
