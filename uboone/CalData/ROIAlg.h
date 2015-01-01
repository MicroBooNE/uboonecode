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
#include "Range.h"

namespace util{

  class ROIALg{
    
  public: 
    typedef short                    Digit;
    typedef std::vector<Digit>       Waveform;
    typedef Waveform::const_iterator Tick;    
    typedef Range<Tick>              Region;
    

    static std::unique_ptr<ROIAlg> MakeROIAlg(fhicl::ParameterSet const&);

    std::string         GetName() { return fAlgName; }
    std::vector<Region> GetSignalRegions();
    std::vector<Region> GetBaselineRegions();
    size_t              GetNSignalRegions();
    size_t              GetNBaselineRegions();

  protected: 
    virtual void AnalyzeWaveform(Waveform const&) = 0;
    void         InsertRegion(Region region, Region_t type);
    void         InsertSignalRegion(Region region);
    void         InsertBaselineRegion(Region region);
    
  private:
    std::string           fAlgName;
    UniqueRangeSet<Tick>  fRangeSet;
    
  };

} //end namespace util

#endif
