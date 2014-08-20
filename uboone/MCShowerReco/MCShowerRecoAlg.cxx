////////////////////////////////////////////////////////////////////////
//
//  MCShowerRecoAlg source
//
////////////////////////////////////////////////////////////////////////

#ifndef MCSHOWERRECOALG_CC
#define MCSHOWERRECOALG_CC

#include "MCShowerRecoAlg.h"

namespace sim {

  //##################################################################
  MCShowerRecoAlg::MCShowerRecoAlg(fhicl::ParameterSet const& pset) 
    : fEdepAlg(pset.get< fhicl::ParameterSet >("MCShowerRecoEdep")),
      fPartAlg(pset.get< fhicl::ParameterSet >("MCShowerRecoPart"))
  //##################################################################
  {
    fDebugMode = pset.get<bool>("DebugMode");
    fG4ModName = pset.get<std::string>("G4ModName");
    fMinShowerEnergy = pset.get<double>("MinShowerEnergy");
    fMinNumDaughters = pset.get<unsigned int>("MinNumDaughters");
  }

  void MCShowerRecoAlg::RecoMCShower(const art::Event &evt)
  {
    fMCShower.clear();

    // Reconstruct shower's particles
    art::Handle<std::vector<simb::MCParticle> > mcpArray;
    evt.getByLabel(fG4ModName, mcpArray);
    fPartAlg.ConstructShower(mcpArray);    

    // Reconstruct shower's energy deposition
    art::Handle<std::vector<sim::SimChannel> > schArray;
    evt.getByLabel(fG4ModName,schArray);
    fEdepAlg.MakeMCShowerEdep(schArray);
    
    // Get shower info from grouped particles
    const std::vector<unsigned int> mother_index = fPartAlg.ShowerMothers();
    fMCShower.reserve(mother_index.size());
    for(size_t shower_index = 0; shower_index < mother_index.size(); ++shower_index) {

      unsigned int mother_candidate = mother_index.at(shower_index);

      ::sim::MCShower shower_prof;
      fPartAlg.Position(mother_candidate, shower_prof.vtxMother);
      fPartAlg.Momentum(mother_candidate, shower_prof.momMother);
      for(auto& v : shower_prof.momMother) v *= 1.e3;
      
      shower_prof.momPdgCode = fPartAlg.PdgCode(mother_candidate);
      shower_prof.momTrackId = fPartAlg.TrackId(mother_candidate);

      shower_prof.daughterTrackId = fPartAlg.ShowerDaughters(shower_index);

      fMCShower.push_back(shower_prof);

    }
    
    // Next, loop over deposited energy and group them into showers
    std::map<unsigned int,size_t> edep_index_map(fEdepAlg.TrackIndexMap());
    for(auto track_iter = edep_index_map.begin();
	track_iter != edep_index_map.end();
	++track_iter) {

      unsigned int part_index = fPartAlg.TrackToParticleIndex((*track_iter).first);
      if(part_index == fPartAlg.kINVALID_UINT) continue;
      unsigned int shower_index = fPartAlg.ShowerIndex(part_index);
      if(shower_index < 0 || shower_index == fPartAlg.kINVALID_UINT) continue;

      for(auto const edep : fEdepAlg.GetEdepArrayAt((*track_iter).second)) {
	
	std::vector<float> vtx(4,0);
	vtx.at(0) = edep.x;
	vtx.at(1) = edep.y;
	vtx.at(2) = edep.z;
	vtx.at(3) = edep.energy;
	fMCShower.at(shower_index).vtxEdep.push_back(vtx);

	// Compute unit vector to this energy deposition
	std::vector<float> mom(4,0);
	mom.at(0) = vtx.at(0) - fMCShower.at(shower_index).vtxMother.at(0);
	mom.at(1) = vtx.at(1) - fMCShower.at(shower_index).vtxMother.at(1);
	mom.at(2) = vtx.at(2) - fMCShower.at(shower_index).vtxMother.at(2);

	float magnitude = sqrt(pow(mom.at(0),2) + pow(mom.at(1),2) + pow(mom.at(2),2));
	mom.at(0) *= edep.energy/magnitude;
	mom.at(1) *= edep.energy/magnitude;
	mom.at(2) *= edep.energy/magnitude;
	mom.at(3) = edep.energy;

	fMCShower.at(shower_index).momDaughter.at(0) += mom.at(0);
	fMCShower.at(shower_index).momDaughter.at(1) += mom.at(1);
	fMCShower.at(shower_index).momDaughter.at(2) += mom.at(2);
	fMCShower.at(shower_index).momDaughter.at(3) += mom.at(3);
	
	fMCShower.at(shower_index).qU   += edep.qU;
	fMCShower.at(shower_index).qV   += edep.qV;
	fMCShower.at(shower_index).qW   += edep.qW;

      }
    }

    if(fDebugMode) {

      for(auto const prof : fMCShower) {
	
	std::cout
	  
	  << Form("  Mother PDG=%d Start @ (%g,%g,%g,%g) with Momentum (%g,%g,%g,%g)",
		  prof.momPdgCode, 
		  prof.vtxMother.at(0), prof.vtxMother.at(1), prof.vtxMother.at(2), prof.vtxMother.at(3),
		  prof.momMother.at(0), prof.momMother.at(1), prof.momMother.at(2), prof.momMother.at(3))
	  << std::endl
	  << Form("    ... with %zu daughters with sum momentum (%g,%g,%g,%g)",
		  prof.daughterTrackId.size(),
		  prof.momDaughter.at(0),prof.momDaughter.at(1),prof.momDaughter.at(2),prof.momDaughter.at(3))
	  << std::endl
	  << Form("    Charge per plane: %g ... %g ... %g", prof.qU, prof.qV, prof.qW)
	  << std::endl
	  << std::endl;
      }
    }
  }

  void MCShowerRecoAlg::Compute3DAngle(const double x, const double y, const double z,
				       double &theta, double&phi) const
  {
    double magnitude = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
    phi = TMath::ATan(z/x);
    theta = TMath::ACos( y / magnitude);
  }
  
  void MCShowerRecoAlg::Compute2DAngle(const geo::View_t view,
				       const double x, const double y, const double z,
				       double angle)
  {
    
  }

}


#endif
