#ifndef MCSHOWERRECOPART_H
#define MCSHOWERRECOPART_H

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
// STL
#include <set>
#include <vector>
#include <sstream>

// ROOT
#include <TString.h>
#include <TTree.h>

namespace sim
{

  class MCShowerRecoPart {

  public:
    static const unsigned int kINVALID_UINT;
    static const int kINVALID_INT;

  public:

    /// Default constructor with fhicl parameters
    MCShowerRecoPart(fhicl::ParameterSet const& pset);
    //ClusterMergeAlg(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

    /// Default destructor
    virtual ~MCShowerRecoPart(){};

    /// Main function to read-in data and fill variables in this algorithm to reconstruct MC shower
    void ConstructShower(const art::Handle<std::vector<simb::MCParticle> > mcpArray);

    /**
       Returns a list ot daughter particle index numbers for the specified shower 
       with the shower index number as an input
     */
    const std::vector<unsigned int>& ShowerDaughters(const unsigned int shower_id) const
    {
      if(shower_id >= _shower_daughters.size()) throw cet::exception(__FUNCTION__) << "Invalid shower index!";
      return _shower_daughters.at(shower_id);
    }

    /**
       Returns a list of shower-mother's particle index. Order respects shower index number.
    */
    const std::vector<unsigned int> ShowerMothers() const
    {
      std::vector<unsigned int> mothers(_shower_index.size(),0);
      for(auto mother_iter = _shower_index.begin(); mother_iter!=_shower_index.end(); ++mother_iter)
	mothers.at((*mother_iter).second) = (*mother_iter).first;
      return mothers;
    }

    //--------------- Particle Information Getters -----------------//

    /**
       Take particle index number and returns shower index number to which this particle belongs.
       Returns -1 if a particle does not belong to any MC shower. Returns kINVALID_INT if input is invalid.
    */
    int ShowerIndex(const unsigned int part_index) const
    {
      if(_shower_id.size() <= part_index) return kINVALID_INT;
      return _shower_id.at(part_index);
    }

    /*
      Take TrackID and returns the corresponding particle unique index number (MCParticle array index)
      Returns kINVALID_UINT if nothing found. 
    */
    unsigned int TrackToParticleIndex(const unsigned int track_id) const
    { 
      auto const iter (_track_index.find(track_id));
      if(iter==_track_index.end()) return kINVALID_UINT;
      return (*iter).second;
    }

    /// Returns particle's PDG code for a specified particle index. Returns kINVALID_INT for invalid part_index.
    int PdgCode(const unsigned int part_index) const
    { if(_pdgcode.size() > part_index) return _pdgcode.at(part_index);
      return kINVALID_INT;
    }

    /// Returns particle's track id for a specified particle index. Returns kINVALID_INT for invalid part_index.
    unsigned int TrackId(const unsigned int part_index) const
    { if(_track_id.size() > part_index) return _track_id.at(part_index);
      return kINVALID_UINT;
    }
    
    /// Returns mother's track ID for a specified particle index. Returns kINVALID_UINT for invalid part_index.
    unsigned int MotherTrackId(const unsigned int part_index) const
    { if(_mother.size() > part_index) return _mother.at(part_index);
      return kINVALID_UINT;
    }

    /**
       Fills the provided vector with particle's position std::vector<double>, (x,y,z,t)
     */
    void Position(const unsigned int part_index, std::vector<double> &xyzt) const
    {
      if(_start_vtx.size() > part_index) xyzt = _start_vtx.at(part_index);
    }
    
    /**
       Fills the provided vector with particle's momentum std::vector<double>, (px,py,pz,E)
    */
    void Momentum(const unsigned int part_index, std::vector<double> &mom) const
    {
      if(_start_mom.size() > part_index) mom = _start_mom.at(part_index);
    }
    
  protected:

    /// lots of stdout stream
    bool _debug_mode;

    //
    // particle-indexed-variables
    //
    /// Track ID => Index Map
    std::map<unsigned int, unsigned int> _track_index;

    /// Track ID
    std::vector<unsigned int> _track_id;

    /// Mother track ID
    std::vector<unsigned int> _mother;

    /// Ancestor track ID
    std::vector<unsigned int> _ancestor;

    /// PDGID
    std::vector<int> _pdgcode;

    /// Stard XYZT
    std::vector<std::vector<double> > _start_vtx;

    /// Start momentum + E
    std::vector<std::vector<double> > _start_mom;

    /// End XYZ
    std::vector<std::vector<double> > _end_vtx;

    /// Set of daughters' track IDs
    std::vector<std::set<unsigned int> > _daughters;

    /// Track index to shower index map
    std::vector<int> _shower_id;


    //
    // shower-indexed-variables
    //
    /// Shower Primary Index ID => Shower Index Map
    std::map<unsigned int, unsigned int> _shower_index;

    /// Shower time-ordered daughters
    std::vector<std::vector<unsigned int> > _shower_daughters;


  protected:

    void GetTrackStartInfo(const unsigned int &index,
			   double &start_x,
			   double &start_y,
			   double &start_z,
			   double &start_time);

    void ClearAndReservePartArray(size_t n);

    void AddParticle(unsigned int track_id,
		     unsigned int mother_track_id,
		     int pdgcode,
		     const TLorentzVector &start_vtx,
		     const TLorentzVector &end_vtx,
		     const std::set<unsigned int> &daughters);

    void AddParticles(const art::Handle<std::vector<simb::MCParticle> > mcpArray);

    void ConstructGranularShower();
    
  }; // class MCShowerRecoPart
  
} //namespace cluster
#endif
