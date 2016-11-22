#include "uboone/"

namespace crt{
  CRTSubRunFilter::CRTSubRunFilter(fhicl::ParameterSet const & pset):
    fRawDigitModuleLabel(p.get<std::string>("RawDigitModuleLabel")),
    art::EDFilter(pset)
  {
        
  }
  CRTSubRunFilter::~CRTSubRunFilter()
  {
    // Check to see if the text file is still open
    // If it is, then <check with either Ketchum or Kirby here>
  }
  void CRTSubRunFilter::beginSubRun(art::SubRun& subrun)
  {
    //open text file        
  }
  bool CRTSubRunFilter::filter(art::Event & e) override
  {
    // extract CRT Hits
    // Check to see if in subrun text file
    // pop the ones that fail the test
  }
  void CRTSubRunFilter::endSubRun(art::SubRun& subrun){
    // Close the text file        
  }
}