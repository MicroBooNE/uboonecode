////////////////////////////////////////////////////////////////////////
//
//  MCShowerRecoEdep source
//
////////////////////////////////////////////////////////////////////////

#ifndef MCSHOWERRECOEDEP_CC
#define MCSHOWERRECOEDEP_CC

#include "MCShowerRecoEdep.h"

namespace sim {

  //const unsigned short MCShowerHit::kINVALID_USHORT = std::numeric_limits<unsigned short>::max();

  //const short MCShowerEdep::kINVALID_SHORT = std::numeric_limits<short>::max();

  //##################################################################
  MCShowerRecoEdep::MCShowerRecoEdep(fhicl::ParameterSet const& pset)
  //##################################################################
  {
    _debug_mode = pset.get<bool>("DebugMode");
    _save_mchit = pset.get<bool>("SaveMCHit");
  }

  void MCShowerRecoEdep::FillTrackList(std::vector<unsigned int> &track_list) const
  {
    track_list.resize(_track_index.size(),0);
    for(auto track_iter = _track_index.begin(); track_iter != _track_index.end(); ++track_iter)
      track_list.at((*track_iter).second)=(*track_iter).first;
  }

  const std::vector<sim::MCShowerEdep>& MCShowerRecoEdep::GetEdepArrayAt(size_t edep_index)
  {
    if(edep_index >= _mc_edeps.size()) 
      
      throw cet::exception(__FUNCTION__) << Form("Track ID %zu not found!",edep_index);

    return _mc_edeps.at(edep_index);
  }

  std::vector<sim::MCShowerEdep>& MCShowerRecoEdep::__GetEdepArray__(unsigned int track_id)
  { 
    if(ExistTrack(track_id)) return _mc_edeps.at((*_track_index.find(track_id)).second);

    _track_index.insert(std::pair<unsigned int,size_t>(track_id,_mc_edeps.size()));
    
    _mc_edeps.push_back(std::vector<sim::MCShowerEdep>());
    
    return (*(_mc_edeps.rbegin()));
  }


  void MCShowerRecoEdep::MakeMCShowerEdep(const art::Handle<std::vector<sim::SimChannel> > schArray)
  {
    if(!schArray.isValid()) throw cet::exception(__FUNCTION__) << "Invalid handle given...";

    _mc_edeps.clear();
    _track_index.clear();

    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::DetectorProperties> detp;

    // Key map to identify a unique particle energy deposition point
    std::map<unsigned short, std::map<double, std::map<double, std::map<double, int> > > > hit_index_m;

    if(_debug_mode) std::cout<<"Processing "<<schArray->size()<<" channels..."<<std::endl;
    // Loop over channels
    for(size_t i=0; i<schArray->size(); ++i) {

      // Get data to loop over
      const art::Ptr<sim::SimChannel> sch_ptr(schArray,i);
      const std::map<unsigned short,std::vector<sim::IDE> > sch_map(sch_ptr->TDCIDEMap());
      // Channel
      UInt_t ch = sch_ptr->Channel();

      // Loop over ticks
      for(auto tdc_iter = sch_map.begin(); tdc_iter!=sch_map.end(); ++tdc_iter) {
	
	unsigned short hit_time = (*tdc_iter).first;

	// Loop over IDEs
	for(auto const ide : (*tdc_iter).second) {

	  int track_id = ide.trackID;
	  if(track_id < 0) track_id = track_id * (-1);
	  unsigned int real_track_id = track_id;

	  //short x_mm = (short)(ide.x+0.5);
	  //short y_mm = (short)(ide.y+0.5);
	  //short z_mm = (short)(ide.z+0.5);
	  double x_mm = ide.x;
	  double y_mm = ide.y;
	  double z_mm = ide.z;

	  int hit_index = -1;
	  auto hit_index_track_iter = hit_index_m.find(real_track_id);
	  if(hit_index_track_iter == hit_index_m.end()) {
	    // create new entry here
	    //_track_index.insert(std::pair<unsigned int,size_t>(real_track_id,_mc_edeps.size()));
	    //_mc_edeps.push_back(std::vector<sim::MCShowerEdep>());
	    int new_hit_index = this->__GetEdepArray__(real_track_id).size();
	    std::map<double,int> tmp_z_index;
	    std::map<double,std::map<double,int> > tmp_y_z_index;
	    std::map<double,std::map<double,std::map<double,int> > > tmp_x_y_z_index;
            tmp_z_index[z_mm] = new_hit_index;
            tmp_y_z_index[y_mm] = tmp_z_index;
            tmp_x_y_z_index[x_mm] = tmp_y_z_index;
            hit_index_m[real_track_id]=tmp_x_y_z_index;

	  }else{

	    auto hit_index_x_iter = (*hit_index_track_iter).second.find(x_mm);

	    if(hit_index_x_iter == (*hit_index_track_iter).second.end()) {
	      // Track exists but this "x" doesn't exist. Create.
	      std::map<double,int> tmp_z_index;
	      std::map<double,std::map<double,int> > tmp_y_z_index;
	      std::map<double,std::map<double,std::map<double,int> > > tmp_x_y_z_index;
	      //int new_hit_index = _mc_edeps.at((*(_track_index.find(real_track_id))).second).size();
	      int new_hit_index = this->__GetEdepArray__(real_track_id).size();
	      tmp_z_index[z_mm] = new_hit_index;
	      tmp_y_z_index[y_mm] = tmp_z_index;
	      (*hit_index_track_iter).second.insert(std::pair<double,std::map<double,std::map<double,int> > >(x_mm,tmp_y_z_index));

	    }else{

	      auto hit_index_y_iter = (*hit_index_x_iter).second.find(y_mm);

	      if(hit_index_y_iter == (*hit_index_x_iter).second.end()) {
		// "x" exists but this "y" doesn't exist. Create.
		std::map<double,int> tmp_z_index;
		//int new_hit_index = _mc_edeps.at((*(_track_index.find(real_track_id))).second).size();
		int new_hit_index = this->__GetEdepArray__(real_track_id).size();
		tmp_z_index[z_mm] = new_hit_index;
		(*hit_index_x_iter).second.insert(std::pair<double,std::map<double,int> >(y_mm,tmp_z_index));

	      }else{

		auto hit_index_z_iter = (*hit_index_y_iter).second.find(z_mm);

		if(hit_index_z_iter == (*hit_index_y_iter).second.end()) {
		  // "y" exists but this "z" does not exist. Create.
		  //int new_hit_index = _mc_edeps.at((*(_track_index.find(real_track_id))).second).size();
		  int new_hit_index = this->__GetEdepArray__(real_track_id).size();
		  (*hit_index_y_iter).second.insert(std::pair<double,int>(z_mm,new_hit_index));
		}

		else

		  hit_index = (*hit_index_z_iter).second;

	      }
	    }
	  }
	  if(hit_index < 0) {
	    // This particle energy deposition is never recorded so far. Create a new Edep
	    MCShowerEdep edep;
	    // Fill Edep
	    edep.x = x_mm;
	    edep.y = y_mm;
	    edep.z = z_mm;
	    //float charge = ide.numElectrons * detp->ElectronsToADC();
	    double charge = ide.numElectrons;

	    geo::View_t view = geom->View(ch);
	    if(view == geo::kU){
	      edep.energy = ide.energy;
	      edep.qU = charge;
	    }
	    else if(view == geo::kV) edep.qV = charge;
	    else if(view == geo::kW || view == geo::kZ) edep.qW = charge;
	    else
	      
	      throw cet::exception(__FUNCTION__) << "Unknown view type: " << view;

	    // If configured to save MC hits, do so
	    if(_save_mchit) {
	      MCShowerHit mchit;
	      mchit.timeStart = hit_time;
	      mchit.timeMax = hit_time;
	      mchit.timeEnd = hit_time;
	      mchit.qSum = charge;
	      mchit.qMax = charge;
	      edep.mchits.insert(std::pair<unsigned short,sim::MCShowerHit>(ch,mchit));
	    }

	    this->__GetEdepArray__(real_track_id).push_back(edep);

	  }else{
	    
	    // Append charge to the relevant edep (@ hit_index)
	    //float charge = ide.numElectrons * detp->ElectronsToADC();
	    double charge = ide.numElectrons;
	    MCShowerEdep &edep = this->__GetEdepArray__(real_track_id).at(hit_index);

	    geo::View_t view = geom->View(ch);
	    if(view == geo::kU) {
	      edep.energy += ide.energy;
	      edep.qU += charge;
	    }
	    else if(view == geo::kV) edep.qV += charge;
	    else if(view == geo::kW || view == geo::kZ) edep.qW += charge;
	    else
	      
	      throw cet::exception(__FUNCTION__) << "Unknown view type: " << view;

	    // If configured to store hit, store
	    if(_save_mchit) {
	      
	      // Check if channel is already stored. If so, append. Else store new.
	      auto hit_iter = edep.mchits.find(ch);
	      if( hit_iter == edep.mchits.end() ) {
		MCShowerHit mchit;
		mchit.timeStart = hit_time;
		mchit.timeMax = hit_time;
		mchit.timeEnd = hit_time;
		mchit.qSum = charge;
		mchit.qMax = charge;
		edep.mchits.insert(std::pair<unsigned short,sim::MCShowerHit>(ch,mchit));

	      }else{

		if( (*hit_iter).second.timeStart > hit_time )
		  throw cet::exception(__FUNCTION__) << "Found hit before the already-registered hit's start time!";
		(*hit_iter).second.qSum += charge;
		(*hit_iter).second.timeEnd = hit_time;
		if( (*hit_iter).second.qMax < charge ) {
		  (*hit_iter).second.qMax = charge;
		  (*hit_iter).second.timeMax = hit_time;
		}
		
	      }
	    }
	  }
	  
	} // end looping over ides in this tick
	
      } // end looping over ticks in this channel
    
    }// end looping over channels

    if(_debug_mode) {

      std::cout<< Form("  Collected %zu particles' energy depositions...",_mc_edeps.size()) << std::endl;
      for(auto const& track_id_index : _track_index ) {

	auto track_id   = track_id_index.first;
	auto edep_index = track_id_index.second;

	std::cout<< Form("    Track ID: %d ... %zu Edep!", track_id, edep_index) << std::endl;

      }
      std::cout<<std::endl;
    }
    


  }



}

#endif
