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

      mini_track.Origin  ( mini_part._origin   );
      mini_track.PdgCode ( mini_part._pdgcode  );
      mini_track.TrackID ( mini_part._track_id );
      mini_track.Process ( mini_part._process  );
      mini_track.Start   ( MCStep( mini_part._start_vtx, mini_part._start_mom ) );
      mini_track.End     ( MCStep( mini_part._end_vtx,   mini_part._end_mom   ) );

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

      mini_track.MotherPdgCode ( mother_part._pdgcode  );
      mini_track.MotherTrackID ( mother_part._track_id );
      mini_track.MotherProcess ( mother_part._process  );
      mini_track.MotherStart   ( MCStep( mother_part._start_vtx, mother_part._start_mom ) );
      mini_track.MotherEnd     ( MCStep( mother_part._end_vtx,   mother_part._end_mom   ) );

      mini_track.AncestorPdgCode ( ancestor_part._pdgcode  );
      mini_track.AncestorTrackID ( ancestor_part._track_id );
      mini_track.AncestorProcess ( ancestor_part._process  );
      mini_track.AncestorStart   ( MCStep( ancestor_part._start_vtx, ancestor_part._start_mom ) );
      mini_track.AncestorEnd     ( MCStep( ancestor_part._end_vtx,   ancestor_part._end_mom   ) );

      for(auto const& vtx_mom : mini_part._det_path)

	mini_track.push_back(MCStep(vtx_mom.first,vtx_mom.second));

      fMCTrack.push_back(mini_track);
    }
    
    if(fDebugMode) {

      for(auto const& prof : fMCTrack) {
	
	std::cout
	  
	  << Form("  Track particle:     PDG=%d Start @ (%g,%g,%g,%g) with Momentum (%g,%g,%g,%g)",
		  prof.PdgCode(),
		  prof.Start().X(),prof.Start().Y(),prof.Start().Z(),prof.Start().T(),
		  prof.Start().Px(),prof.Start().Py(),prof.Start().Pz(),prof.Start().E())
	  << std::endl
	  << Form("    Mother particle:   PDG=%d Start @ (%g,%g,%g,%g) with Momentum (%g,%g,%g,%g)",
		  prof.MotherPdgCode(),
		  prof.MotherStart().X(),prof.MotherStart().Y(),prof.MotherStart().Z(),prof.MotherStart().T(),
		  prof.MotherStart().Px(),prof.MotherStart().Py(),prof.MotherStart().Pz(),prof.MotherStart().E())
	  << std::endl
	  << Form("    Ancestor particle: PDG=%d Start @ (%g,%g,%g,%g) with Momentum (%g,%g,%g,%g)",
		  prof.AncestorPdgCode(),
		  prof.AncestorStart().X(),prof.AncestorStart().Y(),prof.AncestorStart().Z(),prof.AncestorStart().T(),
		  prof.AncestorStart().Px(),prof.AncestorStart().Py(),prof.AncestorStart().Pz(),prof.AncestorStart().E())
	  << std::endl
	  << Form("    ... with %zu trajectory points!",prof.size())
	  << std::endl;

	if(prof.size()) {
	  std::cout 
	    << Form("        Start @ (%g,%g,%g,%g) with Momentum (%g,%g,%g,%g)",
		    prof[0].X(), prof[0].Y(), prof[0].Z(), prof[0].T(),
		    prof[0].Px(), prof[0].Py(), prof[0].Pz(), prof[0].E())
	    << std::endl
	    << Form("        End @ (%g,%g,%g,%g) with Momentum (%g,%g,%g,%g)",
		    (*prof.rbegin()).X(), (*prof.rbegin()).Y(), (*prof.rbegin()).Z(), (*prof.rbegin()).T(),
		    (*prof.rbegin()).Px(), (*prof.rbegin()).Py(), (*prof.rbegin()).Pz(), (*prof.rbegin()).E())
	    << std::endl;
	}
      }

      std::cout<<std::endl<<std::endl;
    }
  }
}


#endif
