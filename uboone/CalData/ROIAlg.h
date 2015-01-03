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

  template <class Digit>
  class ROIAlg{
        
  public: 

    typedef std::vector<Digit> Waveform;
    typedef typename Waveform::const_iterator Tick;
    typedef Range<Tick> Region;
    
    struct SignalBaselineTrio{
      Region BaselineRegion_Pre;
      Region SignalRegion;
      Region BaselineRegion_Post;
    };

    static std::unique_ptr< ROIAlg<Digit> > MakeROIAlg(fhicl::ParameterSet const&);
    ROIAlg(){}
    virtual ~ROIAlg(){}

    std::string         GetName() { return fAlgName; }

    void ProcessWaveform(Waveform const& w){ ProcessWaveform(Region(w.cbegin(),w.cend())); }
    void ProcessWaveform(Region const&);

    UniqueRangeSet<Tick> const& GetSignalRegions();
    const size_t                GetNSignalRegions();
    UniqueRangeSet<Tick> const& GetBaselineRegions();
    const size_t                GetNBaselineRegions();

    void GetSignalAndBaselineRegions( size_t const i_roi,
				      SignalBaselineTrio&);
    void GetAllSignalAndBaselineRegions( std::vector< SignalBaselineTrio >& );
    
  protected: 
    virtual void AnalyzeWaveform(Region const&) = 0;
    void         InsertSignalRegion(Tick const& start, Tick const& end)
    { fSignalRangeSet.emplace(start,end); }
    
  private:
    std::string           fAlgName;
    Tick                  fWaveformStart;
    Tick                  fWaveformEnd;
    UniqueRangeSet<Tick>  fSignalRangeSet;
    UniqueRangeSet<Tick>  fBaselineRangeSet;

    void CreateBaselineRangeSet();
    void ThrowIfNoBaselineRegions();
    void ClearRangeSets();
  };

} //end namespace util

#endif
