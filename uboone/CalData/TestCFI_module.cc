#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geometry/Geometry.h"

#include "CalibrationDBI/Interface/IChannelFilterService.h"
#include "CalibrationDBI/Interface/IChannelFilterProvider.h"

#include <fstream>

class TestCFI : public art::EDAnalyzer {

  public:
    explicit TestCFI(fhicl::ParameterSet const& p) :
      EDAnalyzer(p) {};
    
    void analyze(art::Event const& e) override;
    
};

void TestCFI::analyze(art::Event const& evt) {
  
  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<lariov::IChannelFilterService> cf; 
  const lariov::IChannelFilterProvider& cp = cf->GetFilter();
  
  std::ofstream f("channelStatus.txt",std::ofstream::out);
  
  for (unsigned int i=0; i < geo->Nchannels(); ++i) {
    
    unsigned short status = cp.Status(i);
    f <<i<<" "<<status<<std::endl;
  }
  
  f.close();
}

DEFINE_ART_MODULE(TestCFI)
    
  
  

