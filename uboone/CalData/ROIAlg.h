#ifndef ROIALG_H
#define ROIALG_H

/*!
 * Title:   ROIAlg base class
 * Author:  wketchum@lanl.gov
 * Inputs:  
 * Outputs: 
 *
 * Description:
 * This is a base class for ROIAlgs.
 */

#include <vector>
#include <string>
#include <exception>

#include "fhiclcpp/ParameterSet.h"
#include "Utilities/UniqueRangeSet.h"

namespace util{

  typedef short                    Digit;
  typedef std::vector<Digit>       Waveform;
  typedef Waveform::const_iterator Tick;    
  typedef Range<Tick>              Region;    

  struct SignalBaselineTrio{
    Region BaselineRegion_Pre;
    Region SignalRegion;
    Region BaselineRegion_Post;
  };
  
  class ROIAlg{
    
  public: 
    
    std::unique_ptr<ROIAlg> MakeROIAlg(fhicl::ParameterSet const&);
    virtual ~ROIAlg(){}

    std::string         GetName() { return fAlgName; }

    void                ProcessWaveform(Waveform const&);

    const UniqueRangeSet<Tick> GetSignalRegions();
    const size_t               GetNSignalRegions();
    const UniqueRangeSet<Tick> GetBaselineRegions();
    const size_t               GetNBaselineRegions();

    void GetSignalAndBaselineRegions( size_t const i_roi,
				      SignalBaselineTrio&);
    void GetAllSignalAndBaselineRegions( std::vector< SignalBaselineTrio >& );
    
  protected: 
    virtual void AnalyzeWaveform(Waveform const&) = 0;
    void         InsertSignalRegion(Tick const& start, Tick const& end)
    { fSignalRangeSet.emplace(start,end); }
    
  private:

    
    std::string           fAlgName;
    Tick                  fWaveformStart;
    Tick                  fWaveformEnd;
    UniqueRangeSet<Tick>  fSignalRangeSet;

    UniqueRangeSet<Tick>  CreateBaselineRangeSet() const;
  };

} //end namespace util

#endif
