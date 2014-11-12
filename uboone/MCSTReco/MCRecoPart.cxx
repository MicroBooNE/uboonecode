////////////////////////////////////////////////////////////////////////
//
//  MCRecoPart source
//
////////////////////////////////////////////////////////////////////////

#ifndef MCRECOPART_CXX
#define MCRECOPART_CXX

#include "MCRecoPart.h"

namespace sim {

  //--------------------------------------------------------------------------------------------
  MCRecoPart::MCRecoPart(fhicl::ParameterSet const& pset)
  //--------------------------------------------------------------------------------------------
  { 
    this->clear();
    _track_index.clear();
    _pdg_list.clear();
    for(auto const& id : pset.get<std::vector<int> >("SavePathPDGList"))

      _pdg_list.insert(id);
  }

  //--------------------------------------------------------------------------------------------
  unsigned int MCRecoPart::AncestorTrackID(const unsigned int part_index) const
  //--------------------------------------------------------------------------------------------
  {
    unsigned int result = kINVALID_UINT;

    if(part_index >= this->size()) return result;

    result = this->at(part_index)._mother;

    auto mother_index = TrackToParticleIndex(result);

    while(mother_index != kINVALID_UINT){

      auto const& new_result = this->at(mother_index)._mother;

      if(!new_result) break;

      result = new_result;

      mother_index = TrackToParticleIndex(result);

    }

    return result;
  }

  //--------------------------------------------------------------------------------------------
  void MCRecoPart::AddParticles(const std::vector<simb::MCParticle>& mcp_v,
				const std::vector<simb::Origin_t>&   orig_v)
  //--------------------------------------------------------------------------------------------
  {
    if(orig_v.size() != mcp_v.size()) throw cet::exception(__FUNCTION__) << "MCParticle and Origin_t vector size not same!";

    art::ServiceHandle<geo::Geometry> geo;
    double y_max = geo->DetHalfHeight();
    double y_min = (-1.) * y_max;
    double z_min = 0;
    double z_max = geo->DetLength();
    double x_min = 0;
    double x_max = 2.*(geo->DetHalfWidth());

    this->clear();
    _track_index.clear();

    for(size_t i=0; i < mcp_v.size(); ++i) {
      
      auto const& mcp = mcp_v[i];

      _track_index.insert(std::make_pair((size_t)(mcp.TrackId()),(size_t)(this->size())));

      this->push_back(MCMiniPart());

      auto& mini_mcp = (*this->rbegin());

      for(size_t i=0; i<(size_t)(mcp.NumberDaughters()); ++i)
	mini_mcp._daughters.insert(mcp.Daughter(i));

      mini_mcp._track_id  = mcp.TrackId();
      mini_mcp._pdgcode   = mcp.PdgCode();
      mini_mcp._mother    = mcp.Mother();
      mini_mcp._process   = mcp.Process();
      mini_mcp._start_vtx = mcp.Position();
      mini_mcp._start_mom = mcp.Momentum();
      mini_mcp._end_vtx   = mcp.EndPosition();
      mini_mcp._end_mom   = mcp.EndMomentum();
      mini_mcp._origin    = orig_v[i];

      // Change units to LArSoft (MeV, cm, us)
      for(size_t i=0; i<4; ++i) {
	mini_mcp._start_mom[i] *= 1.e3;
	mini_mcp._end_mom[i]   *= 1.e3;
      }
      /*
      for(size_t i=0; i<3; ++i) {
	mini_mcp._start_vtx[i] /= 10.;
	mini_mcp._end_vtx[i] /= 10.;
      }
      mini_mcp.start_vtx[3] /= 1.e-3;
      mini_mcp.end_vtx[3]   /= 1.e-3;
      */

      if(_pdg_list.find(mcp.PdgCode()) != _pdg_list.end()) {

	std::set<size_t> det_path_index;

	for(size_t i=0; i<mcp.NumberTrajectoryPoints(); ++i) {
	  
	  if( mcp.Vx(i) > x_max || mcp.Vx(i) < x_min || 
	      mcp.Vz(i) > z_max || mcp.Vz(i) < z_min || 
	      mcp.Vy(i) > y_max || mcp.Vy(i) < y_min ) continue;

	  det_path_index.insert(i);
	}

	if(det_path_index.size()) {
	  if( (*det_path_index.begin()) ) 
	    det_path_index.insert( (*det_path_index.begin())-1 );
	  if( det_path_index.size()>1 ) {
	    if( ((*det_path_index.rbegin())+1) < mcp.NumberTrajectoryPoints() ) 
	      det_path_index.insert( (*det_path_index.rbegin())+1 );
	  }
	  mini_mcp._det_path.reserve(det_path_index.size());
	  for(auto const& index : det_path_index) {

	    TLorentzVector vec(mcp.Momentum(index));
	    for(size_t i=0; i<4; ++i) vec[i] *= 1.e3;

	    mini_mcp._det_path.push_back(std::make_pair(mcp.Position(index),vec));

	  }
	}
      }
    }
  }
}
#endif
