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

/*!
 * When you add a new ROIAlg, you should do the following:
 * (1) Forward declare that ROIAlg in the block right below this comment.
 * (2) Define your ROIAlg in its own header file. Keep everything in that one 
 *     header file (don't split into a cxx file).
 * (3) INCLUDE THAT HEADER FILE AT THE BOTTOM OF THIS FILE!
 * (4) Add the fcl parameters for your alg to the roialg.fcl file.
 *     Make sure there is a parameter "AlgName" defined.
 * (5) Add a declaration of your alg in MakeROIAlg below, using your alg
 *     name to pick your particular alg. Needless to say, your alg name 
 *     must be unique then.
 * 
 * Luv ya,
 * Wes
 */
namespace util{
template <class Digit>
  class ROIAlg_DigitAboveThreshold;
}
namespace util{
template <class Digit>
  class ROIAlg_CalibrationPulseFinder;
}

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
     SignalBaselineTrio():
      BaselineRegion_Pre(Tick(),Tick()),
      SignalRegion(Tick(),Tick()),
      BaselineRegion_Post(Tick(),Tick()) {}
    };

    static std::unique_ptr< ROIAlg<Digit> > MakeROIAlg(fhicl::ParameterSet const& p){

      std::string algName = p.get<std::string>("AlgName");
      
      std::unique_ptr< ROIAlg<Digit> > ptr;
      if(algName.compare("DigitAboveThreshold")==0){
	std::unique_ptr< ROIAlg<Digit> > new_ptr(new ROIAlg_DigitAboveThreshold<Digit>(p));
	ptr.swap(new_ptr);
      }
      if(algName.compare("CalibrationPulseFinder")==0){
	std::unique_ptr< ROIAlg<Digit> > new_ptr(new ROIAlg_CalibrationPulseFinder<Digit>(p));
	ptr.swap(new_ptr);
      }
      else{
	std::cout << "Algname is ... " << algName << std::endl;
	throw std::runtime_error("ERROR in ROIAlg: No registered ROIAlg with that name.");
      }
      ptr->ClearRangeSets();
      
      return std::move(ptr);
    }

    ROIAlg(){}
    virtual ~ROIAlg(){}

    std::string         GetName() { return fAlgName; }

    void ProcessWaveform(Waveform const& w){ ProcessWaveform(Region(w.cbegin(),w.cend())); }
    void ProcessWaveform(Region const& r){
      ClearRangeSets();
      fWaveformStart = r.Start();
      fWaveformEnd   = r.End();
      AnalyzeWaveform(r);
      CreateBaselineRangeSet();
      
      if(fBaselineRangeSet.size()!=(fSignalRangeSet.size()+1))
	throw std::runtime_error("ERROR in ROIAlg: BaselineRangeSet size is not exactly one greater than SignalRangeSet size.");
    }


    UniqueRangeSet<Tick> const& GetSignalRegions()    { ThrowIfNoBaselineRegions(); return fSignalRangeSet; }
    const size_t                GetNSignalRegions()   { ThrowIfNoBaselineRegions(); return fSignalRangeSet.size(); }
    UniqueRangeSet<Tick> const& GetBaselineRegions()  { ThrowIfNoBaselineRegions(); return fBaselineRangeSet; }
    const size_t                GetNBaselineRegions() { ThrowIfNoBaselineRegions(); return fBaselineRangeSet.size(); }

    void GetSignalAndBaselineRegions( size_t const i_roi,
				      SignalBaselineTrio& trio){
      ThrowIfNoBaselineRegions();
      
      if( (i_roi+1) > fSignalRangeSet.size() )
	throw std::runtime_error("ERROR in ROIAlg: asked for roi iter that's larger than SignalRangeSet size");
      
      typename UniqueRangeSet<Tick>::iterator iter_signal = fSignalRangeSet.begin();
      typename UniqueRangeSet<Tick>::iterator iter_baseline_prev = fBaselineRangeSet.begin();
      typename UniqueRangeSet<Tick>::iterator iter_baseline_next = fBaselineRangeSet.begin();
      
      std::advance(iter_signal,i_roi);
      std::advance(iter_baseline_prev,i_roi);
      std::advance(iter_baseline_next,i_roi+1);
      
      trio.SignalRegion = *iter_signal;
      trio.BaselineRegion_Pre = *iter_baseline_prev;
      trio.BaselineRegion_Post = *iter_baseline_next;
      
    }

    void GetAllSignalAndBaselineRegions( std::vector< SignalBaselineTrio >& trio_vector){

      ThrowIfNoBaselineRegions();
      
      trio_vector.clear(); trio_vector.reserve( fSignalRangeSet.size() );
      typename UniqueRangeSet<Tick>::iterator iter_signal = fSignalRangeSet.begin();
      typename UniqueRangeSet<Tick>::iterator iter_baseline_prev = fBaselineRangeSet.begin();
      typename UniqueRangeSet<Tick>::iterator iter_baseline_next = fBaselineRangeSet.begin();
      
      for(size_t i_signal=0; i_signal < fSignalRangeSet.size(); i_signal++){
	iter_baseline_next++;
	trio_vector[i_signal].SignalRegion = *iter_signal;
	trio_vector[i_signal].BaselineRegion_Pre = *iter_baseline_prev;
	trio_vector[i_signal].BaselineRegion_Post = *iter_baseline_next;
	iter_signal++;
	iter_baseline_prev++;
      }
    }

    
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

    void CreateBaselineRangeSet(){
      fBaselineRangeSet.clear();
      Tick start,end;
      
      start = fWaveformStart; //initialize start
      //go through signal ranges
      //signal begin = baseline end, and baseline end = signal begin
      for(auto const& range : fSignalRangeSet){
	end = range.Start();
	fBaselineRangeSet.emplace(start,end);
	start = range.End();
      }
      //now need to do last baseline region: last signal end to waveform end
      end = fWaveformEnd;
      fBaselineRangeSet.emplace(start,end);
    }
    
    void ThrowIfNoBaselineRegions(){
      if(fBaselineRangeSet.size()==0)
	throw std::runtime_error("ERROR in ROIAlg: ProcessWaveform not yet called, but trying to access results!");
    }
    void ClearRangeSets(){
      fSignalRangeSet.clear();
      fBaselineRangeSet.clear();
    }
    
  };
  
} //end namespace util

#include "ROIAlg_DigitAboveThreshold.h"
#include "ROIAlg_CalibrationPulseFinder.h"

#endif
