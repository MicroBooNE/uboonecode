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
#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"
#include "Geometry/Geometry.h"

#include "MCShowerRecoPart.h"
#include "MCShowerRecoEdep.h"
#include "MCBase/MCShower.h"

// STL
#include <set>
#include <vector>
#include <sstream>

// ROOT
#include <TString.h>
#include <TTree.h>

namespace sim
{

  class MCShowerRecoAlg {

  public:

    /// Default constructor with fhicl parameters
    MCShowerRecoAlg(fhicl::ParameterSet const& pset);

    /// Default destructor
    virtual ~MCShowerRecoAlg(){};

    void RecoMCShower(const art::Event &evt);

    size_t NumShowers() const
    {return fMCShower.size();}

    void ValidateIndex(const size_t shower_index) const
    {
      if(shower_index >= fMCShower.size()) 
	throw cet::exception(__FUNCTION__) << Form("Invalid shower index %zu",shower_index);
    }

    const ::sim::MCShower& ShowerProfile(const size_t shower_index) const
    {
      ValidateIndex(shower_index);
      return fMCShower.at(shower_index);
    }

    const std::vector<unsigned int> &Daughters(const size_t shower_index) const
    {
      ValidateIndex(shower_index);
      return fMCShower.at(shower_index).DaughterTrackID();
    }

    int ShowerIndex(const unsigned int track_id) const
    {
      unsigned int part_index = fPartAlg.TrackToParticleIndex(track_id);
      if(part_index == MCShowerRecoPart::kINVALID_UINT) 
	throw cet::exception(__FUNCTION__) << Form("Track ID %d not found in the list...",track_id);

      return fPartAlg.ShowerIndex(part_index);
    }
    
    void Compute3DAngle(const double x, const double y, const double z,
			double &theta, double&phi) const;

    void Compute2DAngle(const geo::View_t view,
			const double x, const double y, const double z,
			double angle);

  protected:

    bool             fDebugMode;
    std::string      fG4ModName;
    MCShowerRecoEdep fEdepAlg;
    MCShowerRecoPart fPartAlg;
    
    std::vector<sim::MCShower> fMCShower;

    double fMinShowerEnergy;
    unsigned int fMinNumDaughters;

  }; // class MCShowerHitRecoAlg
  
} //namespace cluster
#endif
