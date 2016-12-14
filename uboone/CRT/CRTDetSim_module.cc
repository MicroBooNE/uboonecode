#include "uboone/CRT/CRTData.hh"
#include "uboone/CRT/CRTDetSim.hh"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/AuxDetGeometry.h"

#include <vector>
#include <memory>
#include <sstream>
#include <string>
#include <math.h>

namespace crt{

  CRTDetSim::CRTDetSim(const fhicl::ParameterSet& pSet): 
    fThreshold(pSet.get<uint32_t>("Threshold", 1)),
    fConversionFactor(pSet.get<float>("Calibration", 1.e4)),
    fT1Precision(pSet.get<float>("T1Precision", 1.e4)),
    fProducerName(pSet.get<std::string>("ProducerName", "largeant"))
  {
    produces< std::vector<CRTData> >();
  }

  CRTDetSim::~CRTDetSim()
  {

  }

  void CRTDetSim::produce(art::Event& evt)
  {
    art::ServiceHandle<geo::AuxDetGeometry> geo;
    art::Handle< std::vector<sim::AuxDetSimChannel> > channels;
    evt.getByLabel(this->fProducerName,channels);
    if(!channels.isValid()){
      mf::LogWarning(__log_name__)<<"Cannot get the AuxDetChannels";
      return;
    }

    mf::LogInfo(__log_name__)<<" Number of Channels Hit: "<<channels->size();
    std::unique_ptr< std::vector<CRTData> > hits(new std::vector<CRTData>);
    for(auto it = channels->begin(); it!= channels->end(); ++it){
      uint32_t id = it->AuxDetID();
      uint32_t sens_id = it->AuxDetSensitiveID();
      mf::LogInfo("CRTDetSim")<<"Found AuxDetData: "<<id<<" , "<<sens_id;
      std::vector< sim::AuxDetIDE > ides = it->AuxDetIDEs();
      mf::LogInfo("CRTDetSim")<<"Number of IDEs in this event: "<<ides.size();

      for(auto ideIt = ides.begin(); ideIt!= ides.end(); ++ideIt){
        float adc = (this->fConversionFactor * ideIt->energyDeposited); 
        if(adc<fThreshold) continue;
        /// t0 is currently computed as the average time betwen entry and exit.
        float t0 = (ideIt->entryT+ideIt->exitT)/2.0;
        /// t1 is computed as the same with lower precision
        float t1 = trunc(t0/fT1Precision)*fT1Precision;
        CRTData dat(id, t0, t1, adc);
        hits->push_back(dat);
      }
    }
    evt.put(std::move(hits));
  }

  DEFINE_ART_MODULE(CRTDetSim)
}
