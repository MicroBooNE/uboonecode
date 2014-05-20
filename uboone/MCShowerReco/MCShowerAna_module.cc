
#ifndef MCShowerAna_H
#define MCShowerAna_H

#include "MCShowerRecoAlg.h"

// LArSoft includes
//#include "Geometry/Geometry.h"

// ART includes.
#include "art/Framework/Core/EDAnalyzer.h"

namespace larreco {
 
  class MCShowerAna : public art::EDAnalyzer{
  public:
 
    MCShowerAna(const fhicl::ParameterSet&);
    virtual ~MCShowerAna();

    void beginJob();

    void analyze (const art::Event&); 

  private:

    MCShowerRecoAlg fShowerRecoAlg;

  };

} 

#endif//  MCShowerAna_H

// MCShowerAna.cc

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

namespace larreco {
  DEFINE_ART_MODULE(MCShowerAna)
}


// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Simulation/SimListUtils.h"

namespace larreco {

  //-----------------------------------------------------------------------
  // Constructor
  MCShowerAna::MCShowerAna(fhicl::ParameterSet const& pset) 
    : EDAnalyzer(pset),
      fShowerRecoAlg(pset.get< fhicl::ParameterSet >("MCShowerRecoAlg"))
  {
    
  }

  //-----------------------------------------------------------------------
  // Destructor
  MCShowerAna::~MCShowerAna(){}
   
  //-----------------------------------------------------------------------
  void MCShowerAna::beginJob(){}
   

  //-----------------------------------------------------------------------
  void MCShowerAna::analyze(const art::Event& evt) 
  {
    fShowerRecoAlg.RecoMCShower(evt);
  }
} // namespace opdet


