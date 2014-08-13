////////////////////////////////////////////////////////////////////////
//
//  MCShowerRecoPart source
//
////////////////////////////////////////////////////////////////////////

#ifndef MCSHOWERRECOPART_CC
#define MCSHOWERRECOPART_CC

#include "MCShowerRecoPart.h"

namespace sim {

  const unsigned int MCShowerRecoPart::kINVALID_UINT = std::numeric_limits<unsigned int>::max();
  const int MCShowerRecoPart::kINVALID_INT = std::numeric_limits<int>::max();

  //##################################################################
  MCShowerRecoPart::MCShowerRecoPart(fhicl::ParameterSet const& pset)
  //##################################################################
  {
    _debug_mode = pset.get<bool>("DebugMode");
  }

  //----------------------------------------------------------------------------------------------- 
  void MCShowerRecoPart::ConstructShower(const art::Handle<std::vector<simb::MCParticle> > mcpArray)
  //-----------------------------------------------------------------------------------------------
  {
    if(!mcpArray.isValid()) 

      throw cet::exception(__FUNCTION__) << "Invalid particle array provided. Nothing done!";
    
    AddParticles(mcpArray);

    ConstructGranularShower();

  }

  //----------------------------------------------------------------
  void MCShowerRecoPart::ClearAndReservePartArray(size_t n)
  //----------------------------------------------------------------
  {
    _track_index.clear();
    _track_id.clear();
    _mother.clear();
    _ancestor.clear();
    _pdgcode.clear();
    _daughters.clear();
    _shower_id.clear();
    _start_vtx.clear();
    _start_mom.clear();

    _track_id.reserve(n);
    _mother.reserve(n);
    _pdgcode.reserve(n);
    _start_vtx.reserve(n);
    _start_mom.reserve(n);
    _daughters.reserve(n);
    _shower_id.reserve(n);
  }

  //----------------------------------------------------------------
  void MCShowerRecoPart::AddParticle(unsigned int track_id,
				    unsigned int mother_track_id,
				    int pdgcode,
				    const TLorentzVector &start_vtx,
				    const TLorentzVector &start_mom,
				    const std::set<unsigned int> &daughters)
  //----------------------------------------------------------------
  {
    unsigned int index = _track_index.size();

    _track_index.insert(std::pair<unsigned int, unsigned int>(track_id,index));
    _track_id.push_back(track_id);
    _mother.push_back(mother_track_id);
    _pdgcode.push_back(pdgcode);
    _start_vtx.push_back(std::vector<double>(4,0.));
    _start_mom.push_back(std::vector<double>(4,0.));
    _daughters.push_back(std::set<unsigned int>(daughters));

    _start_vtx.at(index).at(0)=start_vtx.X();
    _start_vtx.at(index).at(1)=start_vtx.Y();
    _start_vtx.at(index).at(2)=start_vtx.Z();
    _start_vtx.at(index).at(3)=start_vtx.T();

    _start_mom.at(index).at(0)=start_mom.X();
    _start_mom.at(index).at(1)=start_mom.Y();
    _start_mom.at(index).at(2)=start_mom.Z();
    _start_mom.at(index).at(3)=start_mom.T();
    
  }
  
  //--------------------------------------------------------------------------------------------
  void MCShowerRecoPart::AddParticles(const art::Handle<std::vector<simb::MCParticle> > mcpArray)
  //--------------------------------------------------------------------------------------------
  {

    //art::Handle<std::vector<simb::MCParticle> > mcpArray;
    //evt.getByLabel(fG4ModName,mcpArray);

    // Read in all particles' information
    ClearAndReservePartArray(mcpArray->size());
    for(size_t i=0; i < mcpArray->size(); ++i) {

      const art::Ptr<simb::MCParticle> mcp_ptr(mcpArray,i);


      std::set<unsigned int> daughters;
      for(size_t i=0; i<(size_t)(mcp_ptr->NumberDaughters()); ++i)
	daughters.insert(mcp_ptr->Daughter(i));

      this->AddParticle((unsigned int)(mcp_ptr->TrackId()),
			(unsigned int)(mcp_ptr->Mother()),
			mcp_ptr->PdgCode(),
			mcp_ptr->Position(),
			mcp_ptr->Momentum(),
			daughters);
    }
  }

  void MCShowerRecoPart::GetTrackStartInfo(const unsigned int &index,
					  double &start_x,
					  double &start_y,
					  double &start_z,
					  double &start_time)
  {
    if(index > _track_id.size())

      throw cet::exception(__FUNCTION__) << Form("Particle index %d not found!",index);

    start_x = _start_vtx.at(index).at(0);
    start_y = _start_vtx.at(index).at(1);
    start_z = _start_vtx.at(index).at(2);
    start_time = _start_vtx.at(index).at(3);

  }

  //----------------------------------------------------------------
  void MCShowerRecoPart::ConstructGranularShower()
  //----------------------------------------------------------------
  {

    if(!_mother.size()) return;

    _shower_id.clear();
    _shower_id.resize(_start_vtx.size(),-1);
    _shower_index.clear();
    _shower_daughters.clear();

    // Construct MCShower
    std::vector<std::multimap<double,unsigned int> > daughter_map;
    for(size_t i=0; i<_track_id.size(); ++i) {

      int candidate_mom_index=-1;
      if( _pdgcode.at(i) == 22 ||
          _pdgcode.at(i) == 11 ||
          _pdgcode.at(i) == -11 )
	candidate_mom_index = i;

      unsigned int mom_track = _mother.at(i);
      auto mom_iter = _track_index.find(mom_track);
      while(mom_iter != _track_index.end()) {

        unsigned int mom_index = (*mom_iter).second;

        if( _pdgcode.at(mom_index) == 22 || _pdgcode.at(mom_index) == 11 || _pdgcode.at(mom_index) == -11 )

          candidate_mom_index = mom_index;
	
        mom_iter = _track_index.find(_mother.at(mom_index));

      }

      if(candidate_mom_index >= 0) {
	
	auto candidate_mom_iter = _shower_index.find(candidate_mom_index);
	if(candidate_mom_iter == _shower_index.end()) {
	  _shower_index.insert(std::pair<unsigned int,unsigned int>(candidate_mom_index,_shower_index.size()));
	  daughter_map.push_back(std::multimap<double,unsigned int>());
	}
	unsigned int shower_index = (*_shower_index.find(candidate_mom_index)).second;
	daughter_map.at(shower_index).insert(std::pair<double,unsigned int>(_start_vtx.at(i).at(3),i));
	_shower_id.at(i) = shower_index;
	
      } else if(_debug_mode) {

	std::cout
	  << "Found a particle that does not belong to a shower!" << std::endl
	  << Form(" PDGID: %d ... Track %d @ (%g,%g,%g,%g) with (%g,%g,%g,%g)",
		  _pdgcode.at(i),
		  _track_id.at(i),
		  _start_vtx.at(i).at(0),
		  _start_vtx.at(i).at(1),
		  _start_vtx.at(i).at(2),
		  _start_vtx.at(i).at(3),
		  _start_mom.at(i).at(0),
		  _start_mom.at(i).at(1),
		  _start_mom.at(i).at(2),
		  _start_mom.at(i).at(3))
	  << std::endl << std::endl;

      }
    }


    if(_debug_mode)
      std::cout
	<< Form("Found %zu MCShowers....",_shower_index.size()) << std::endl;

    _shower_daughters.resize(_shower_index.size(),std::vector<unsigned int>());
    for(const auto &mom : _shower_index) {

      _shower_daughters.at(mom.second).reserve(daughter_map.at(mom.second).size());
      for(auto const &part_index : daughter_map.at(mom.second))

	_shower_daughters.at(mom.second).push_back(part_index.second);

      if(_debug_mode) 
	std::cout 
	  << Form("PDGID: %d ... Track %d @ (%g,%g,%g,%g) with (%g,%g,%g,%g) ... %zu daughters!",
		  _pdgcode.at(mom.first),
		  _track_id.at(mom.first),
		  _start_vtx.at(mom.first).at(0),
		  _start_vtx.at(mom.first).at(1),
		  _start_vtx.at(mom.first).at(2),
		  _start_vtx.at(mom.first).at(3),
		  _start_mom.at(mom.first).at(0),
		  _start_mom.at(mom.first).at(1),
		  _start_mom.at(mom.first).at(2),
		  _start_mom.at(mom.first).at(3),
		  _shower_daughters.at(mom.second).size())
	  << std::endl;
    }

    if(_debug_mode)
      std::cout<<std::endl;
  }
  
}
#endif
