#include "uboone/CRT/CRTSubRunFilter.hh"

namespace crt{
  CRTSubRunFilter::CRTSubRunFilter(fhicl::ParameterSet const & pset):
    //fRawDigitModuleLabel(p.get<std::string>("RawDigitModuleLabel")),
    art::EDFilter()
  {
        
  }
  CRTSubRunFilter::~CRTSubRunFilter()
  {
    // Check to see if the text file is still open
    // If it is, then <check with either Ketchum or Kirby here>
  }
  bool CRTSubRunFilter::beginSubRun(art::SubRun& subrun)
  {
    //open text file        
    return false;
  }
  bool CRTSubRunFilter::filter(art::Event & evt)
  {
    // extract CRT Hits
    // Check to see if in subrun text file
    // pop the ones that fail the test
    return false;
  }
  bool CRTSubRunFilter::endSubRun(art::SubRun& subrun){
    // Close the text file       
    return false;
  }
}
