#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"

#include "larevt/CalibrationDBI/Interface/IChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/IChannelStatusProvider.h"
#include "larevt/CalibrationDBI/IOVData/ChannelStatus.h"

#include <fstream>

class TestCFI : public art::EDAnalyzer {

  public:
    explicit TestCFI(fhicl::ParameterSet const& p) :
      EDAnalyzer(p) {};
    
    void analyze(art::Event const& e) override;
    
};

void TestCFI::analyze(art::Event const& evt) {
  
  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<lariov::IChannelStatusService> cf; 
  const lariov::IChannelStatusProvider& cp = cf->GetProvider();
  
  std::ofstream f("channelStatus.txt",std::ofstream::out);
  
  for (unsigned int i=0; i < geo->Nchannels(); ++i) {
    
    unsigned short status = cp.Status(i);
    if (status == lariov::kNOISY) f <<i<<" "<<status<<std::endl;
  }
  
  f.close();
}

DEFINE_ART_MODULE(TestCFI)
    
  
  

