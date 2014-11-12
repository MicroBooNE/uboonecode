////////////////////////////////////////////////////////////////////////
//
//  MCTrackRecoAlg source
//
////////////////////////////////////////////////////////////////////////

#ifndef MCTRACKRECOALG_CXX
#define MCTRACKRECOALG_CXX

#include "MCTrackRecoAlg.h"

namespace sim {

  //##################################################################
  MCTrackRecoAlg::MCTrackRecoAlg(fhicl::ParameterSet const& pset) 
  //##################################################################
  {
    fDebugMode = pset.get<bool>("DebugMode");
  }

  void MCTrackRecoAlg::Reconstruct(const MCRecoPart& part_v)
  {

    fMCTrack.clear();
    for(size_t i=0; i<part_v.size(); ++i) {

      auto const& mini_part = part_v[i];

      if( part_v._pdg_list.find(mini_part._pdgcode) == part_v._pdg_list.end() ) continue;

      ::sim::MCTrack mini_track;
      
      mini_track.PdgCode   ( mini_part._pdgcode  );
      mini_track.G4TrackID ( mini_part._track_id );
      mini_track.Process   ( mini_part._process  );
      mini_track.G4Start   ( MCStep( mini_part._start_vtx, mini_part._start_mom ) );
      mini_track.G4End     ( MCStep( mini_part._end_vtx,   mini_part._end_mom   ) );

      unsigned int mother_index   = i;
      MCMiniPart   mother_part    = mini_part;
      unsigned int ancestor_index = i;
      MCMiniPart   ancestor_part  = mini_part;

      if(mini_part._mother) {
	mother_index   = part_v.TrackToParticleIndex(mini_part._mother);
	mother_part    = part_v.at(mother_index);
	ancestor_index = part_v.TrackToParticleIndex(part_v.AncestorTrackID(i));
	ancestor_part  = part_v.at(ancestor_index);
      }

      mini_track.MotherPdgCode   ( mother_part._pdgcode  );
      mini_track.MotherG4TrackID ( mother_part._track_id );
      mini_track.MotherProcess   ( mother_part._process  );
      mini_track.MotherG4Start   ( MCStep( mother_part._start_vtx, mother_part._start_mom ) );
      mini_track.MotherG4End     ( MCStep( mother_part._end_vtx,   mother_part._end_mom   ) );

      mini_track.AncestorPdgCode   ( ancestor_part._pdgcode  );
      mini_track.AncestorG4TrackID ( ancestor_part._track_id );
      mini_track.AncestorProcess   ( ancestor_part._process  );
      mini_track.AncestorG4Start   ( MCStep( ancestor_part._start_vtx, ancestor_part._start_mom ) );
      mini_track.AncestorG4End     ( MCStep( ancestor_part._end_vtx,   ancestor_part._end_mom   ) );

      for(auto const& vtx_mom : mini_part._det_path)

	mini_track.push_back(MCStep(vtx_mom.first,vtx_mom.second));

      fMCTrack.push_back(mini_track);
    }
    
    if(fDebugMode) {

      for(auto const& prof : fMCTrack) {
	
	std::cout
	  
	  << Form("  Track particle:     PDG=%d Start @ (%g,%g,%g,%g) with Momentum (%g,%g,%g,%g)",
		  prof.PdgCode(),
		  prof.G4Start().X(),prof.G4Start().Y(),prof.G4Start().Z(),prof.G4Start().T(),
		  prof.G4Start().Px(),prof.G4Start().Py(),prof.G4Start().Pz(),prof.G4Start().E())
	  << std::endl
	  << Form("    Mother particle:   PDG=%d Start @ (%g,%g,%g,%g) with Momentum (%g,%g,%g,%g)",
		  prof.MotherPdgCode(),
		  prof.MotherG4Start().X(),prof.MotherG4Start().Y(),prof.MotherG4Start().Z(),prof.MotherG4Start().T(),
		  prof.MotherG4Start().Px(),prof.MotherG4Start().Py(),prof.MotherG4Start().Pz(),prof.MotherG4Start().E())
	  << std::endl
	  << Form("    Ancestor particle: PDG=%d Start @ (%g,%g,%g,%g) with Momentum (%g,%g,%g,%g)",
		  prof.AncestorPdgCode(),
		  prof.AncestorG4Start().X(),prof.AncestorG4Start().Y(),prof.AncestorG4Start().Z(),prof.AncestorG4Start().T(),
		  prof.AncestorG4Start().Px(),prof.AncestorG4Start().Py(),prof.AncestorG4Start().Pz(),prof.AncestorG4Start().E())
	  << std::endl
	  << Form("    ... with %zu trajectory points!",prof.size())
	  << std::endl
	  << Form("        Start @ (%g,%g,%g,%g) with Momentum (%g,%g,%g,%g)",
		  prof[0].X(), prof[0].Y(), prof[0].Z(), prof[0].T(),
		  prof[0].Px(), prof[0].Py(), prof[0].Pz(), prof[0].E())
	  << std::endl
	  << Form("        End @ (%g,%g,%g,%g) with Momentum (%g,%g,%g,%g)",
		  (*prof.rbegin()).X(), (*prof.rbegin()).Y(), (*prof.rbegin()).Z(), (*prof.rbegin()).T(),
		  (*prof.rbegin()).Px(), (*prof.rbegin()).Py(), (*prof.rbegin()).Pz(), (*prof.rbegin()).E())
	  << std::endl;
      }

      std::cout<<std::endl<<std::endl;
    }
  }
}


#endif
