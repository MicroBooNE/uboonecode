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

#include "MCShowerRecoPart.h"
#include "MCShowerRecoEdep.h"

// STL
#include <set>
#include <vector>
#include <sstream>

// ROOT
#include <TString.h>
#include <TTree.h>

namespace larreco
{

  class MCShowerProfile {

  public:

    static const int kINVALID_INT;
    static const unsigned int kINVALID_UINT;

  public:

    MCShowerProfile(){ Clear(); }
    ~MCShowerProfile(){}

    //---- Mother info ----//
    int momPdgCode;                 ///< mother PDG code
    unsigned int momTrackId;        ///< mother G4 Track ID
    std::vector<double> vtxMother;  ///< mother position 4-vector @ generation
    std::vector<double> momMother;  ///< mother momentum 4-vector @ generation
    /// mother 3D angle phi (along shower angle definition, not ordinary coord. system)
    double phiMother;
    /// mother 3D angle theta (along shower angle definition, not ordinary coord. system)
    double thetaMother;
    /// mother 2D angle on U-plane
    double uAngleMother;
    /// mother 2D angle on V-plane
    double vAngleMother;
    /// mother 2D angle on W-plane
    double wAngleMother;

    //---- Daughter info ----//
    std::vector<unsigned int> daughterTrackId; ///< Daughters' track ID
    std::vector<float> momDaughter;            ///< Daughters' deposit sum momentum 4-vector
    /// daughter 3D angle phi (along shower angle definition, not ordinary coord. system)
    float phiDaughter;
    /// daughter 3D angle theta (along shower angle definition, not ordinary coord. system)
    float thetaDaughter;
    /// daughter 2D angle on U-plane
    float uAngleDaughter;
    /// daughter 2D angle on V-plane
    float vAngleDaughter;
    /// daughter 2D angle on W-plane
    float wAngleDaughter;

    //---- Charge per plane ----//
    float qU; ///< Charge deposit on U plane
    float qV; ///< Charge deposit on V plane
    float qW; ///< Charge deposit on W plane

    /// Energy deposit point (x,y,z,dE) by daughters
    std::vector<std::vector<float> > vtxEdep;

    void Clear() {

      momPdgCode = kINVALID_INT;
      momTrackId = kINVALID_UINT;      
      vtxMother.clear();
      vtxMother.resize(4,0);
      momMother.clear();
      momMother.resize(4,0);
      thetaMother = phiMother = TMath::Pi()*4;
      uAngleMother = vAngleMother = wAngleMother = TMath::Pi()*4;

      daughterTrackId.clear();
      momDaughter.clear();
      momDaughter.resize(4,0);
      phiDaughter = thetaDaughter = TMath::Pi()*4;
      uAngleDaughter = vAngleDaughter = wAngleDaughter = TMath::Pi()*4;

      qU = qV = qW = 0;
      vtxEdep.clear();

    }

  };

  class MCShowerRecoAlg {

  public:

    /// Default constructor with fhicl parameters
    MCShowerRecoAlg(fhicl::ParameterSet const& pset);

    /// Default destructor
    virtual ~MCShowerRecoAlg(){};

    void RecoMCShower(const art::Event &evt);

    size_t NumShowers() const
    {return fShowerProf.size();}

    void ValidateIndex(const size_t shower_index) const
    {
      if(shower_index >= fShowerProf.size()) 
	throw cet::exception(__FUNCTION__) << Form("Invalid shower index %zu",shower_index);
    }

    const larreco::MCShowerProfile& ShowerProfile(const size_t shower_index) const
    {
      ValidateIndex(shower_index);
      return fShowerProf.at(shower_index);
    }

    const std::vector<unsigned int> &Daughters(const size_t shower_index) const
    {
      ValidateIndex(shower_index);
      return fShowerProf.at(shower_index).daughterTrackId;
    }

    int ShowerIndex(const unsigned int track_id) const
    {
      unsigned int part_index = fPartAlg.TrackToParticleIndex(track_id);
      if(part_index == MCShowerRecoPart::kINVALID_UINT) 
	throw cet::exception(__FUNCTION__) << Form("Track ID %d not found in the list...",track_id);

      return fPartAlg.ShowerIndex(part_index);
    }
    
    const std::vector<float> HitPurity(recob::Hit hit) const;

    const std::vector<float> ClusterPurity(recob::Hit hit) const;

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
    
    std::vector<larreco::MCShowerProfile> fShowerProf;

    double fMinShowerEnergy;
    unsigned int fMinNumDaughters;

  }; // class MCShowerHitRecoAlg
  
} //namespace cluster
#endif
