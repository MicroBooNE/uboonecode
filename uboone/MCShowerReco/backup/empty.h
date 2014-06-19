#ifndef MCSHOWERRECOALG_H
#define MCSHOWERRECOALG_H

// ART includes
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"

// LArSoft
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/SpacePoint.h"
#include "RecoAlg/SpacePointAlg.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"
#include "Geometry/Geometry.h"
// STL
#include <set>
#include <vector>
#include <sstream>

// ROOT
#include <TString.h>
#include <TTree.h>

namespace larreco
{

  class MCShowerRecoAlg {

  public:

    /// Default constructor with fhicl parameters
    MCShowerRecoAlg(fhicl::ParameterSet const& pset);

    /// Default destructor
    virtual ~MCShowerRecoAlg(){};

  protected:


  }; // class MCShowerHitRecoAlg
  
} //namespace cluster
#endif
