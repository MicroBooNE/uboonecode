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
    
    art::ServiceHandle<geo::Geometry> geo;

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
    std::map<size_t,int> stored_mcs_index;
    for(size_t shower_index = 0; shower_index < mother_index.size(); ++shower_index) {

      unsigned int mother_candidate = mother_index.at(shower_index);

      auto mcs_index_iter = stored_mcs_index.insert(std::make_pair(shower_index,-1));

      ::sim::MCShower shower_prof;
      std::vector<double> mother_vtx;
      std::vector<double> mother_mom;
      fPartAlg.Position(mother_candidate, mother_vtx);
      fPartAlg.Momentum(mother_candidate, mother_mom);
      for(auto& v : mother_mom) v *= 1.e3;

      if(fDebugMode)

	std::cout << "Found MCShower with mother energy: " << mother_mom.at(3) << " MeV";

      // Skip if mother energy is less than the enery threshold
      if(mother_mom.at(3) < fMinShowerEnergy) {
	if(fDebugMode)
	  std::cout << " ... below energy threshold: skipping!"<<std::endl;
	continue;
      }else if(fPartAlg.ShowerDaughters(shower_index).size() < fMinNumDaughters) {
	if(fDebugMode)
	  std::cout << " ... below # daughter particle count threshold: skipping!"<<std::endl;
	continue;
      }

      (*mcs_index_iter.first).second = fMCShower.size();
      if(fDebugMode)
	
	std::cout << " index " << fMCShower.size() << " => " << shower_index
		  << std::endl;

      shower_prof.MotherPosition(mother_vtx);
      shower_prof.MotherMomentum(mother_mom);
      
      shower_prof.MotherPDGID   ( fPartAlg.PdgCode(mother_candidate) );
      shower_prof.MotherTrackID ( fPartAlg.TrackId(mother_candidate) );
      shower_prof.MotherCreationProcess( fPartAlg.Process(mother_candidate) );

      std::vector<unsigned int> daughter_track_id;
      daughter_track_id.reserve( fPartAlg.ShowerDaughters(shower_index).size() );
      std::vector<double>       daughter_vtx;

      bool first = true;
      for(auto const& index : fPartAlg.ShowerDaughters(shower_index)) {

	daughter_track_id.push_back(fPartAlg.TrackId(index));
	if(first) {
	  fPartAlg.Position(index, daughter_vtx);
	  first=false;
	}

      }

      shower_prof.DaughterTrackID(daughter_track_id);
      shower_prof.DaughterPosition(daughter_vtx);

      fMCShower.push_back(shower_prof);
    }
    
    // Next, loop over deposited energy and group them into showers
    std::map<unsigned int,size_t> edep_index_map(fEdepAlg.TrackIndexMap());

    std::vector<std::vector<double> > mcs_daughter_mom_v(fMCShower.size(),std::vector<double>(4,0));
    std::vector<std::vector<double> > plane_charge_v(fMCShower.size(),std::vector<double>(geo->Nplanes(),0));

    for(auto track_iter = edep_index_map.begin();
	track_iter != edep_index_map.end();
	++track_iter) {

      unsigned int part_index = fPartAlg.TrackToParticleIndex((*track_iter).first);
      if(part_index == fPartAlg.kINVALID_UINT) continue;
      unsigned int shower_index = fPartAlg.ShowerIndex(part_index);
      if(shower_index < 0 || shower_index == fPartAlg.kINVALID_UINT) continue;

      auto mcs_index_iter = stored_mcs_index.find(shower_index);
      if(mcs_index_iter == stored_mcs_index.end()) {
	std::cerr<<"Logic error: invalid index of stored MCShower!"<<std::endl;
	throw std::exception();
      }

      if((*mcs_index_iter).second < 0) continue;

      int edep_index = fEdepAlg.TrackToEdepIndex((*track_iter).first);
      if(edep_index < 0) continue;

      auto& daughter_mom = mcs_daughter_mom_v [(*mcs_index_iter).second];
      auto& plane_charge = plane_charge_v     [(*mcs_index_iter).second];

      for(auto const& edep : fEdepAlg.GetEdepArrayAt(edep_index)) {
	
	std::vector<double> vtx(4,0);
	vtx.at(0) = edep.x;
	vtx.at(1) = edep.y;
	vtx.at(2) = edep.z;
	vtx.at(3) = edep.energy;

	//fMCShower.at(shower_index).vtxEdep.push_back(vtx);

	// Compute unit vector to this energy deposition
	std::vector<double> mom(4,0);
	mom.at(0) = vtx.at(0) - fMCShower.at((*mcs_index_iter).second).DaughterPosition().at(0);
	mom.at(1) = vtx.at(1) - fMCShower.at((*mcs_index_iter).second).DaughterPosition().at(1);
	mom.at(2) = vtx.at(2) - fMCShower.at((*mcs_index_iter).second).DaughterPosition().at(2);

	double magnitude = sqrt(pow(mom.at(0),2) + pow(mom.at(1),2) + pow(mom.at(2),2));
	mom.at(0) *= (edep.energy/magnitude);
	mom.at(1) *= (edep.energy/magnitude);
	mom.at(2) *= (edep.energy/magnitude);
	mom.at(3) = edep.energy;

	daughter_mom[0] += mom.at(0);
	daughter_mom[1] += mom.at(1);
	daughter_mom[2] += mom.at(2);
	daughter_mom[3] += mom.at(3);

	plane_charge.at(geo::kU) += edep.qU;
	plane_charge.at(geo::kV) += edep.qV;
	plane_charge.at(geo::kW) += edep.qW;
	
      }
    }

    // Store plane charge & daughter momentum
    for(size_t mcs_index=0; mcs_index<fMCShower.size(); ++mcs_index) {

      // Correct for energy deposition normalization
      auto& daughter_mom = mcs_daughter_mom_v[mcs_index];

      double magnitude = sqrt(pow(daughter_mom[0],2)+pow(daughter_mom[1],2)+pow(daughter_mom[2],2));

      daughter_mom[0] *= (daughter_mom[3]/magnitude);
      daughter_mom[1] *= (daughter_mom[3]/magnitude);
      daughter_mom[2] *= (daughter_mom[3]/magnitude);

      fMCShower.at(mcs_index).DaughterMomentum(daughter_mom);
      fMCShower.at(mcs_index).Charge(plane_charge_v[mcs_index]);
    }

    if(fDebugMode) {

      for(auto const& prof : fMCShower) {
	
	std::cout
	  
	  << Form("  Mother PDG=%d Start @ (%g,%g,%g,%g) with Momentum (%g,%g,%g,%g)",
		  prof.MotherPDGID(),
		  prof.MotherPosition().at(0), prof.MotherPosition().at(1), prof.MotherPosition().at(2), prof.MotherPosition().at(3),
		  prof.MotherMomentum().at(0), prof.MotherMomentum().at(1), prof.MotherMomentum().at(2), prof.MotherMomentum().at(3))
	  << std::endl
	  << Form("    ... with %zu daughters with sum momentum (%g,%g,%g,%g)",
		  prof.DaughterTrackID().size(),
		  prof.DaughterMomentum().at(0),prof.DaughterMomentum().at(1),prof.DaughterMomentum().at(2),prof.DaughterMomentum().at(3))
	  << std::endl
	  << "    Charge per plane: ";

	for(size_t i=0; i<prof.Charge().size(); ++i) {

	  std::cout << " | Plane " << i << std::flush;
	  std::cout << " ... Q = " << prof.Charge(i) << std::flush;

	}

	std::cout<<std::endl<<std::endl;
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
