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

namespace util{

  class ROIALg{

    struct setcomp{
      bool operator() ( Region const& lhs, Region const& rhs) const
	return (lhs.first.first < rhs.first.first);
    };
    
  public: 
    typedef short                    Digit;
    typedef std::vector<Digit>       Waveform;
    typedef Waveform::const_iterator Tick;    
    typedef std::pair<Tick,Tick>     Region;

    typedef enum{
      kBackground,
      kSignal
    } Region_t;

    static std::unique_ptr<ROIAlg> MakeROIAlg(fhicl::ParameterSet const&);

    std::string         GetName() { return fAlgName; }
    std::vector<Region> GetSignalRegions();
    std::vector<Region> GetBackgroundRegions();
    size_t              GetNSignalRegions();
    size_t              GetNBackgroundRegions();

  protected: 
    virtual void AnalyzeWaveform(Waveform const&) = 0;
    void         InsertRegion( Region region, Region_t type, Tick hint=fRegions.cend());
    
  private:
    std::string                                     fAlgName;
    std::set< std::pair<Region,Region_t>, setcomp > fRegions;

    
    
  };

} //end namespace util

#endif
