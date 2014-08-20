////////////////////////////////////////////////////////////////////////
// Class:       WireROIInspector
// Module Type: analyzer
// File:        WireROIInspector_module.cc
//
// Generated at Tue Aug 12 07:05:12 2014 by Kazuhiro Terao using artmod
// from cetpkgsupport v1_05_04.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include <iostream>

#include "RecoBase/Wire.h"

class WireROIInspector;

class WireROIInspector : public art::EDAnalyzer {
public:
  explicit WireROIInspector(fhicl::ParameterSet const & p);
  virtual ~WireROIInspector();

  void analyze(art::Event const & e) override;


private:

  // Declare member data here.

};


WireROIInspector::WireROIInspector(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

WireROIInspector::~WireROIInspector()
{
  // Clean up dynamic memory and other resources here.
}

void WireROIInspector::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  art::Handle<std::vector<recob::Wire> > wireHandle;
  e.getByLabel("caldata",wireHandle);

  if(wireHandle.isValid()) {

    for(auto const& w : *wireHandle) {

      auto const& roi = w.SignalROI();

      auto const& ranges = roi.get_ranges();
      
      std::cout<<roi.size()<<" ... "<<ranges.size()<<std::endl;
      for(auto const& r : ranges){
	
	std::cout<<r.begin_index()<<" => "<<r.end_index()<<std::endl;
	for(size_t i=r.begin_index(); i<=r.end_index(); ++i)

	  std::cout<< roi[i] << " " << std::flush;

	std::cout<<std::endl;

      }
	//std::cout<<r.first<<" : "<<r.second<<std::endl;

    }
  }

}

DEFINE_ART_MODULE(WireROIInspector)
