
#ifndef UBMCSHOWERRECO_H
#define UBMCSHOWERRECO_H

// Framework
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft
#include "Simulation/SimChannel.h"
#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Shower.h"
#include "AnalysisBase/Calorimetry.h"
#include "Geometry/Geometry.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
// ROOT
#include <TStopwatch.h>
#include <TH1D.h>
#include <TTree.h>

//#include "DataFormat-TypeDef.hh"
//#include "Base-TypeDef.hh"

// STD
#include <iostream>

namespace ana {

  struct TmpHitContainer_t {

    double x,y,z;
    unsigned int start_time;
    unsigned int peak_time;
    unsigned int end_time;
    unsigned int wire_number;
    double sum_charge;
    double max_charge;

    TmpHitContainer_t()
    {
      start_time = peak_time = end_time = 0;
      sum_charge = max_charge = -1;
      x = y = z = 0;
      wire_number = 1e9;
    }

  };
 
  class UBMCShowerReco : public art::EDProducer{
  public:
 
    UBMCShowerReco(const fhicl::ParameterSet&);
    virtual ~UBMCShowerReco();

    void beginJob();

    void endJob();

    void produce (art::Event&); 

  private:

    void GetTrackInfo(const unsigned int &index,
		      double &start_x,
		      double &start_y,
		      double &start_z,
		      double &start_time,
		      double &end_x,
		      double &end_y,
		      double &end_z,
		      double &end_time);

    void GetTrackStartInfo(const unsigned int &index,
			   double &start_x,
			   double &start_y,
			   double &start_z,
			   double &start_time);

    void GetTrackEndInfo(const unsigned int &index,
			   double &end_x,
			   double &end_y,
			   double &end_z,
			   double &end_time);

    void GetDistanceAndShowerAngle(const std::vector<double> &vtx1, const std::vector<double> &vtx2,
				   double &dist, double &theta, double &phi) const;

    void GetDistanceAndSphericalAngle(const std::vector<double> &vtx1, const std::vector<double> &vtx2,
				     double &dist, double &theta, double &phi) const;
    
    void ClearShowerStorage();

    void ConstructGranularShower(const art::Event&);
    
    void CombineGranularShower();

    void MakeMCShowerHits(const art::Event &evt, 
			  std::vector<recob::Shower> &showers_v,
			  std::vector<std::vector<recob::Cluster> >&clusters_v,
			  std::vector<std::vector<std::vector<recob::Hit> > >&hits_v,
			  std::vector<double> &shower_energy_v);

    /// lots of stdout stream
    bool _debug_mode;

    /// time separation cut
    double _comb_time_cut;

    /// spatial separation cut
    double _comb_dist2_cut;

    std::string fG4ModName;

    //
    // particle-indexed-variables
    //
    /// Track ID => Index Map
    std::map<unsigned int, unsigned int> _track_index;

    /// Track ID
    std::vector<unsigned int> _track_id;

    /// Mother track ID
    std::vector<unsigned int> _mother;

    /// PDGID
    std::vector<int> _pdgcode;

    /// Stard XYZ
    std::vector<std::vector<double> > _start_vtx;

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

  };

} 

#endif//  UBMCShowerReco_H

// UBMCShowerReco.cc

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

namespace ana {
  DEFINE_ART_MODULE(UBMCShowerReco)
}

namespace ana {

  //-----------------------------------------------------------------------
  // Constructor
  UBMCShowerReco::UBMCShowerReco(fhicl::ParameterSet const& pset) 
  {
    _debug_mode = pset.get<bool>("DebugMode");
    fG4ModName = pset.get<std::string>("G4ModName");
    _comb_time_cut  = pset.get<double>("TimeCut");
    _comb_dist2_cut = pset.get<double>("Dist2Cut");

    produces< std::vector<recob::Hit> >();
    produces< std::vector<recob::Shower> >();
    produces< std::vector<anab::Calorimetry> >();    
    //produces< std::vector<simb::MCTruth> >();
    produces< std::vector<recob::Cluster> >();

    produces< art::Assns<recob::Shower, recob::Cluster>  >();
    //produces< art::Assns<recob::Shower, anab::Calorimetry> >();
    produces< art::Assns<recob::Cluster, recob::Hit>  >();

  }

  //-----------------------------------------------------------------------
  // Destructor
  UBMCShowerReco::~UBMCShowerReco(){}
   
  //-----------------------------------------------------------------------
  void UBMCShowerReco::beginJob()
  {
    art::ServiceHandle<art::TFileService>  fileService;
  }

  //-----------------------------------------------------------------------
  void UBMCShowerReco::endJob(){

  }
   
  //--------------------------------------
  void UBMCShowerReco::ClearShowerStorage()
  //--------------------------------------
  {

    // Parcile IDs
    _track_index.clear();
    _track_id.clear();
    _mother.clear();
    _pdgcode.clear();
    _daughters.clear();
    _start_vtx.clear();
    _end_vtx.clear();

    // Shower IDs
    _shower_index.clear();
    _shower_daughters.clear();
    _shower_id.clear();
    // Charge
  }


  //-----------------------------------------------------------------------
  void UBMCShowerReco::produce(art::Event& evt) 
  {

    std::unique_ptr<std::vector<anab::Calorimetry> > out_calorimetries_v(new std::vector<anab::Calorimetry>);
    std::unique_ptr<std::vector<recob::Shower> > out_showers_v(new std::vector<recob::Shower>);
    std::unique_ptr<std::vector<recob::Cluster> > out_clusters_v(new std::vector<recob::Cluster>);
    std::unique_ptr<std::vector<recob::Hit> > out_hits_v(new std::vector<recob::Hit>);

    std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> > assn_cluster_hit(new art::Assns<recob::Cluster, recob::Hit>);
    std::unique_ptr< art::Assns<recob::Shower, recob::Cluster> > assn_shower_cluster(new art::Assns<recob::Shower, recob::Cluster>);
    //std::unique_ptr< art::Assns<recob::Shower, anab::Calorimetry> > assn_shower_calorimetry(new art::Assns<recob::Shower, anab::Calorimetry>);

    std::vector<recob::Shower>     showers_v;
    std::vector<simb::MCTruth>     mctruths_v;
    std::vector<anab::Calorimetry> calorimetries_v;
    std::vector<std::vector<recob::Cluster> > clusters_v;
    std::vector<std::vector<std::vector<recob::Hit> > > hits_v;

    art::Handle<std::vector<simb::MCParticle> > mcpArray;
    evt.getByLabel(fG4ModName,mcpArray);

    if(!mcpArray.isValid()) return;

    // Clear event-wise data product
    ClearShowerStorage();

    // Construct granular showers from track daughter-mother PDGID
    ConstructGranularShower(evt);

    // Combine showers 
    CombineGranularShower();

    // Create MC Shower Hits
    std::vector<double> shower_energy_v; /// this holds energy sum from all particles for this shower
    MakeMCShowerHits(evt, showers_v, clusters_v, hits_v, shower_energy_v);

    // Let's store now
    std::vector<std::pair<size_t,size_t> > hit_cluster_assn;
    std::vector<std::pair<size_t,size_t> > cluster_shower_assn;
    for(size_t shower_index=0; shower_index<showers_v.size(); ++shower_index) {

      if(!(clusters_v.at(shower_index).size())) continue;

      size_t cluster_shower_assn_start = out_clusters_v->size();
      size_t cluster_shower_assn_end = out_clusters_v->size();
      for(size_t cluster_index=0; cluster_index<clusters_v.at(shower_index).size(); ++cluster_index) {
	
	out_clusters_v->push_back(clusters_v.at(shower_index).at(cluster_index));

	cluster_shower_assn_end = out_clusters_v->size() - 1;

	size_t hit_cluster_assn_start = out_hits_v->size();
	size_t hit_cluster_assn_end = out_hits_v->size();
	for(size_t hit_index=0; hit_index<hits_v.at(shower_index).at(cluster_index).size(); ++hit_index) {

	  out_hits_v->push_back(hits_v.at(shower_index).at(cluster_index).at(hit_index));
	  hit_cluster_assn_end = out_hits_v->size() - 1;

	}
	hit_cluster_assn.push_back(std::pair<size_t,size_t>(hit_cluster_assn_start,hit_cluster_assn_end));
	
      }
      cluster_shower_assn.push_back(std::pair<size_t,size_t>(cluster_shower_assn_start, cluster_shower_assn_end));

      out_showers_v->push_back(showers_v.at(shower_index));

      out_calorimetries_v->push_back(anab::Calorimetry(shower_energy_v.at(shower_index),
						       std::vector<double>(),
						       std::vector<double>(),
						       std::vector<double>(),
						       std::vector<double>(),
						       0,0));
    }

    if(hit_cluster_assn.size()!=out_clusters_v->size() ||
       cluster_shower_assn.size()!=out_showers_v->size() )

      throw cet::exception(__FUNCTION__) << "Logic error for storing associations!";

    for(size_t i=0; i<out_showers_v->size(); ++i) {

      util::CreateAssn(*this, evt,
		       //*(out_showers_v.get()),
		       //*(out_clusters_v.get()),
		       *out_showers_v,
		       *out_clusters_v,
		       *(assn_shower_cluster.get()),
		       cluster_shower_assn.at(i).first,
		       cluster_shower_assn.at(i).second,
		       i);

    }

    for(size_t i=0; i<out_clusters_v->size(); ++i) {

      util::CreateAssn(*this, evt,
		       //*(out_clusters_v.get()),
		       //*(out_hits_v.get()),
		       *out_clusters_v,
		       *out_hits_v,
		       *(assn_cluster_hit.get()),
		       hit_cluster_assn.at(i).first,
		       hit_cluster_assn.at(i).second,
		       i);
    }

    // Store it
    evt.put(std::move(out_hits_v));
    evt.put(std::move(out_clusters_v));
    evt.put(std::move(out_showers_v));
    evt.put(std::move(out_calorimetries_v));

    evt.put(std::move(assn_shower_cluster));
    evt.put(std::move(assn_cluster_hit));

    // Create MC Shower Clusters
    //MakeMCShowerClusters(hits_v, clusters_v);

    // Create 3D Shower Truth
    //MakeMCShowers(evt,showers_v);

  }

  void UBMCShowerReco::GetDistanceAndShowerAngle(const std::vector<double> &vtx1, const std::vector<double> &vtx2,
						 double &dist, double &theta, double &phi) const
  {
    // Note these theta/phi are for showers
    dist  = sqrt( pow(vtx2.at(0) - vtx1.at(0),2) +
		  pow(vtx2.at(1) - vtx1.at(1),2) +
		  pow(vtx2.at(2) - vtx1.at(2),2) );
    phi   = TMath::ATan( (vtx2.at(2) - vtx1.at(2)) / (vtx2.at(0) - vtx1.at(0)));
    theta = TMath::ACos( (vtx2.at(1) - vtx1.at(1)) / dist );
  }

  void UBMCShowerReco::GetDistanceAndSphericalAngle(const std::vector<double> &vtx1, const std::vector<double> &vtx2,
						    double &dist, double &theta, double &phi) const
  {
    // Note these theta/phi are for spherical coordinate system
    dist  = sqrt( pow(vtx2.at(0) - vtx1.at(0),2) +
		  pow(vtx2.at(1) - vtx1.at(1),2) +
		  pow(vtx2.at(2) - vtx1.at(2),2) );
    theta = TMath::ACos( (vtx2.at(2) - vtx1.at(2)) / dist );
    phi   = TMath::ACos( (vtx2.at(1) - vtx1.at(1)) / (vtx2.at(0) - vtx1.at(0)) );
  }


  void UBMCShowerReco::MakeMCShowerHits(const art::Event &evt, 
					std::vector<recob::Shower> &showers_v,
					std::vector<std::vector<recob::Cluster> >&clusters_v,
					std::vector<std::vector<std::vector<recob::Hit> > >&hits_v,
					std::vector<double> &shower_energy_v)
  {
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::DetectorProperties> detp;
    art::Handle<std::vector<sim::SimChannel> > schArray;
    evt.getByLabel(fG4ModName,schArray);
    if(!schArray.isValid()) return;

    // Shower-wise variables initialized here
    shower_energy_v.clear();
    shower_energy_v.resize(_shower_index.size(),0);
    std::vector<std::vector<double> > shower_vtx(_shower_index.size(),std::vector<double>());
    std::vector<std::vector<std::vector<TmpHitContainer_t> > > tmp_hits_v(_shower_index.size(),
									  std::vector<std::vector<TmpHitContainer_t> >(geom->Nplanes(),
														       std::vector<TmpHitContainer_t>())
									  );
    hits_v.clear();
    hits_v.resize(_shower_index.size(),std::vector<std::vector<recob::Hit> >());
    clusters_v.clear();
    clusters_v.resize(_shower_index.size(),std::vector<recob::Cluster>());
    showers_v.clear();
    showers_v.reserve(_shower_index.size());
    for(auto shower_index_iter = _shower_index.begin();
	shower_index_iter     != _shower_index.end();
	++shower_index_iter) {

      clusters_v[(*shower_index_iter).second].reserve(geom->Nplanes());
      hits_v[(*shower_index_iter).second].resize(geom->Nplanes(),std::vector<recob::Hit>());
      shower_vtx[(*shower_index_iter).second] = _start_vtx.at((*shower_index_iter).first); 

    }
    if(_debug_mode) std::cout<<"Processing "<<schArray->size()<<" channels..."<<std::endl;
    // Loop over channels
    for(size_t i=0; i<schArray->size(); ++i) {

      // Get data to loop over
      const art::Ptr<sim::SimChannel> sch_ptr(schArray,i);
      const std::map<unsigned short,std::vector<sim::IDE> > sch_map(sch_ptr->TDCIDEMap());
      // Channel
      UInt_t ch = sch_ptr->Channel();
      // Wire ID
      std::vector<geo::WireID> wids = geom->ChannelToWire(ch);
      geo::WireID wire_id = wids[0];
      // SignalType
      geo::SigType_t signal_type = geom->SignalType(ch);
      // View
      geo::View_t view = geom->View(ch);
      // Plane
      size_t plane=wire_id.Plane;
      // Channel



      // Initialize variables to construct hits
      std::map<unsigned int, std::map<double, std::map<double, std::map<double, int> > > > hit_index_m;
      std::vector<unsigned int> hit_shower_index_v;
      std::vector<std::vector<double> > hit_vtx_v;
      std::vector<std::vector<unsigned int> > hit_time_v;
      std::vector<std::vector<double> > hit_charge_v;

      // Loop over ticks
      for(auto tdc_iter = sch_map.begin(); tdc_iter!=sch_map.end(); ++tdc_iter) {

	unsigned int hit_time = (*tdc_iter).first;

	// Loop over IDEs
	for(auto const ide : (*tdc_iter).second) {

	  int track_id = ide.trackID;
	  if(track_id < 0) track_id = track_id * (-1);
	  unsigned int real_track_id = track_id;

	  auto const track_index_iter = _track_index.find(real_track_id);
	  if(track_index_iter == _track_index.end()) {

	    std::cerr << " Unknown track ID in IDE: "<< ide.trackID << std::endl;
	  
	    continue;
	  }
	  
	  const int shower_index = _shower_id.at((*track_index_iter).second);
	  
	  if(shower_index < 0) {
	    //std::cout << Form("Ignoring sim::SimChannel for non-shower track: %d",(*track_index_iter).first)<<std::endl;
	    continue;
	  }else if((size_t)shower_index >= _shower_index.size()){
	    std::cerr << Form("INVALID SHOWER INDEX: %d ... track = %d, PDGCODE=%d",
			      shower_index, real_track_id, _pdgcode.at((*track_index_iter).second))<<std::endl;

	    for(auto tmpitr = _shower_index.begin(); tmpitr!=_shower_index.end(); ++tmpitr)
	      std::cerr<<(*tmpitr).first<< " => "<<(*tmpitr).second<<std::endl;
	    throw cet::exception("Logic Error");
	    continue;
	  }
	  shower_energy_v[shower_index] += ide.energy;

	  int hit_index = -1;
	  auto const hit_index_track_iter = hit_index_m.find(real_track_id);
	  if(hit_index_track_iter != hit_index_m.end()) {
	    
	    auto const hit_index_x_iter = (*hit_index_track_iter).second.find(ide.x);

	    if(hit_index_x_iter != (*hit_index_track_iter).second.end()) {

	      auto const hit_index_y_iter = (*hit_index_x_iter).second.find(ide.y);

	      if(hit_index_y_iter != (*hit_index_x_iter).second.end()) {
		
		auto const hit_index_z_iter = (*hit_index_y_iter).second.find(ide.z);

		if(hit_index_z_iter != (*hit_index_y_iter).second.end()) {

		  hit_index = (*hit_index_z_iter).second;
		}
	      }
	    }
	  }

	  if(hit_index < 0) {
	    // Record this IDE as the 1st sample
	    std::map<double,int> tmp_z_index;
	    std::map<double,std::map<double,int> > tmp_y_z_index;
	    std::map<double,std::map<double,std::map<double,int> > > tmp_x_y_z_index;
	    tmp_z_index[ide.z] = (int)(hit_shower_index_v.size());
	    tmp_y_z_index[ide.y] = tmp_z_index;
	    tmp_x_y_z_index[ide.x] = tmp_y_z_index;
	    hit_index_m[real_track_id]=tmp_x_y_z_index;
	    
	    hit_shower_index_v.push_back(shower_index);
	    hit_time_v.push_back(std::vector<unsigned int>(1,hit_time));
	    hit_charge_v.push_back(std::vector<double>(1,ide.numElectrons * detp->ElectronsToADC()));
	    // Record this particle-step-wise hit 
	    std::vector<double> vtx(3,0);
	    vtx[0] = ide.x;
	    vtx[1] = ide.y;
	    vtx[2] = ide.z;
	    hit_vtx_v.push_back(vtx);
	  }else{
	    // append this hit time and charge
	    hit_time_v[hit_index].push_back(hit_time);
	    hit_charge_v[hit_index].push_back(ide.numElectrons * detp->ElectronsToADC());
	  }
	  
	} // end looping over ides in this tick

      } // end looping over ticks in this channel

      // Now let's loop over found set of energy depositions and construct hits
      // If we find more than one hit at the same tdc and same shower id, we combine those hits.
      std::map<unsigned int,std::map<unsigned int,TmpHitContainer_t> > tmp_hit_m;
      for(size_t hit_index=0; hit_index<hit_shower_index_v.size(); ++hit_index) {

	TmpHitContainer_t tmp_hit;
	tmp_hit.start_time=(*(hit_time_v.at(hit_index).begin()));
	tmp_hit.end_time=(*(hit_time_v.at(hit_index).rbegin()));
	tmp_hit.peak_time=tmp_hit.start_time;
	tmp_hit.max_charge=0;
	tmp_hit.sum_charge=0;
	tmp_hit.x = hit_vtx_v.at(hit_index).at(0);
	tmp_hit.y = hit_vtx_v.at(hit_index).at(1);
	tmp_hit.z = hit_vtx_v.at(hit_index).at(2);
	tmp_hit.wire_number = wids[0].Wire;
	for(size_t tick=0; tick<hit_time_v.at(hit_index).size(); ++tick) {

	  double this_charge = hit_charge_v.at(hit_index).at(tick);
	  tmp_hit.sum_charge += this_charge;
	  if(tmp_hit.max_charge < this_charge) {
	    tmp_hit.max_charge = this_charge;

	    tmp_hit.peak_time  = hit_time_v.at(hit_index).at(tick);
	  }

	}
	tmp_hits_v[hit_shower_index_v.at(hit_index)][plane].push_back(tmp_hit);

	// Check if hit already exists or not & insert
	auto const tmp_hit_shower_iter = tmp_hit_m.find(hit_shower_index_v.at(hit_index));
	if(tmp_hit_shower_iter != tmp_hit_m.end()) {

	  auto const tmp_hit_tick_iter = (*tmp_hit_shower_iter).second.find(tmp_hit.peak_time);
	  if( tmp_hit_tick_iter != (*tmp_hit_shower_iter).second.end() ) {

	    // This hit already exists. Let's combine them.
	    TmpHitContainer_t prev_hit = (*tmp_hit_tick_iter).second;
	    tmp_hit.sum_charge += prev_hit.sum_charge;
	    tmp_hit.max_charge += prev_hit.max_charge;

	    tmp_hit.start_time  = ( tmp_hit.start_time < prev_hit.start_time ? tmp_hit.start_time : prev_hit.start_time );
	    tmp_hit.end_time    = ( tmp_hit.end_time   > prev_hit.end_time   ? tmp_hit.end_time   : prev_hit.end_time   );
	    
	    tmp_hit_m[hit_shower_index_v.at(hit_index)][tmp_hit.start_time] = tmp_hit;
	  }
	  // This hit does not exist yet
	  tmp_hit_m[hit_shower_index_v.at(hit_index)].insert(std::pair<unsigned int,TmpHitContainer_t>(tmp_hit.start_time,tmp_hit) );

	}else{
	  tmp_hit_m[hit_shower_index_v.at(hit_index)]=std::map<unsigned int,TmpHitContainer_t>();
	  tmp_hit_m[hit_shower_index_v.at(hit_index)].insert(std::pair<unsigned int,TmpHitContainer_t>(tmp_hit.start_time,tmp_hit) );
	}

      }	

      // Now insert hits into a storage container
      for(auto hit_shower_iter = tmp_hit_m.begin();
	  hit_shower_iter != tmp_hit_m.end();
	  ++hit_shower_iter) {

	unsigned int shower_index = (*hit_shower_iter).first;
	for(auto hit_tick_iter = (*hit_shower_iter).second.begin();
	    hit_tick_iter != (*hit_shower_iter).second.end();
	    ++hit_tick_iter) {

	  hits_v[shower_index][plane].push_back(recob::Hit(view,signal_type,wire_id,
							   (double)((*hit_tick_iter).second.start_time), 0,
							   (double)((*hit_tick_iter).second.end_time),   0,
							   (double)((*hit_tick_iter).second.peak_time),  0,
							   (*hit_tick_iter).second.sum_charge,        0,
							   (*hit_tick_iter).second.max_charge,        0,
							   1,1) );
	}
      }
    }// end looping over channels


    if(_debug_mode){
      for(size_t shower_index=0; shower_index<tmp_hits_v.size(); shower_index++) {
	for(size_t plane_index=0; plane_index<tmp_hits_v.at(shower_index).size(); ++plane_index)
	  
	  std::cout<<Form(" Plane %zu ... %zu hits...",plane_index,tmp_hits_v.at(shower_index).at(plane_index).size())<<std::endl;
	
      }
    }

    // Next, process cluster ingredients
    // Loop over showers
    for(size_t shower_index=0; shower_index<tmp_hits_v.size(); shower_index++) {
      
      double shower_charge=0;
      double shower_phi=0;
      double shower_theta=0;
      double shower_dcosvtx[3]={0.};
      double shower_dcosvtx_err[3]={0.};
      double shower_width_xy[2]={0.};
      double shower_length=0.;

      for(size_t plane_index=0; plane_index<geom->Nplanes(); plane_index++) {

	if(!tmp_hits_v.at(shower_index).at(plane_index).size()) {
	  if(_debug_mode) 
	    std::cout << Form("Encountered empty hits on plane %zu for shower %zu!!!",plane_index,shower_index)<<std::endl;
	  continue;
	}
	  //throw cet::exception(__FUNCTION__) << Form("Encountered empty hits on plane %zu for shower %zu!!!",plane_index,shower_index);

	std::vector<double> hit_vtx(3,0);
	double cluster_phi    = 0;
	double cluster_theta  = 0;
	double cluster_charge = 0;
	double start_time = 0;
	double end_time   = 0;
	unsigned int start_wire = 0;
	unsigned int end_wire   = 0;
	double dist,theta,phi;
	double min_dist=1e10;
	double max_dist=0;
	for(auto const tmp_hit : tmp_hits_v.at(shower_index).at(plane_index)) {
	  hit_vtx[0] = tmp_hit.x;
	  hit_vtx[1] = tmp_hit.y;
	  hit_vtx[2] = tmp_hit.z;
	  GetDistanceAndShowerAngle(shower_vtx.at(shower_index),hit_vtx,dist,theta,phi);

	  cluster_phi    += phi * tmp_hit.sum_charge;
	  cluster_theta  += theta * tmp_hit.sum_charge;
	  cluster_charge += tmp_hit.sum_charge;
	  if(dist < min_dist) {
	    min_dist = dist;
	    start_time = tmp_hit.peak_time;
	    start_wire = tmp_hit.wire_number;
	  }
	  if(dist > max_dist) {
	    max_dist = dist;
	    end_time = tmp_hit.peak_time;
	    end_wire = tmp_hit.wire_number;
	  }
	}
	shower_phi    += cluster_phi;
	shower_theta  += cluster_theta;
	shower_charge += cluster_charge;

	cluster_phi = cluster_phi / cluster_charge;
	clusters_v[shower_index].push_back(recob::Cluster((double)(start_wire), 0,
							  start_time, 0,
							  (double)(end_wire), 0,
							  end_time, 0,
							  shower_phi, 0,
							  0, 0,
							  cluster_charge,
							  (geo::View_t)(plane_index),
							  (shower_index * geom->Nplanes() + plane_index)));
      } // end looping over planes

      shower_theta = shower_theta / shower_charge;
      shower_phi   = shower_phi / shower_charge;
      shower_dcosvtx[0] = TMath::Sin(shower_theta) * TMath::Cos(shower_phi);
      shower_dcosvtx[1] = TMath::Sin(shower_theta) * TMath::Sin(shower_phi);
      shower_dcosvtx[2] = TMath::Cos(shower_phi);


      showers_v.push_back(recob::Shower(shower_dcosvtx,
					shower_dcosvtx_err,
					shower_width_xy,
					shower_length,
					shower_charge,
					showers_v.size()) );
    }
  }

  //----------------------------------------------------------------
  void UBMCShowerReco::ConstructGranularShower(const art::Event& evt)
  //----------------------------------------------------------------
  {

    art::Handle<std::vector<simb::MCParticle> > mcpArray;
    evt.getByLabel(fG4ModName,mcpArray);

    std::set<unsigned int> shower_parent_id;

    _track_id.reserve(mcpArray->size());
    _mother.reserve(mcpArray->size());
    _pdgcode.reserve(mcpArray->size());
    _start_vtx.reserve(mcpArray->size());
    _end_vtx.reserve(mcpArray->size());
    _daughters.reserve(mcpArray->size());

    // Read in all particles' information

    for(size_t i=0; i < mcpArray->size(); ++i) {

      const art::Ptr<simb::MCParticle> mcp_ptr(mcpArray,i);

      _track_index.insert(std::pair<unsigned int, unsigned int>(mcp_ptr->TrackId(), i));
      _track_id.push_back(mcp_ptr->TrackId());
      _mother.push_back(mcp_ptr->Mother());
      _pdgcode.push_back(mcp_ptr->PdgCode());

      std::vector<double> start(4,0);
      std::vector<double> end(4,0);
      start[0] = mcp_ptr->Vx();
      start[1] = mcp_ptr->Vy();
      start[2] = mcp_ptr->Vz();
      start[3] = mcp_ptr->T();
      end[0]   = mcp_ptr->EndX();
      end[1]   = mcp_ptr->EndY();
      end[2]   = mcp_ptr->EndZ();
      end[3]   = mcp_ptr->EndT();

      _start_vtx.push_back(start);
      _end_vtx.push_back(end);

      std::set<unsigned int> daughters;
      for(size_t i=0; i<(size_t)(mcp_ptr->NumberDaughters()); ++i)
	
	daughters.insert(mcp_ptr->Daughter(i));

      _daughters.push_back(daughters);

      if( mcp_ptr->PdgCode() == 22 ||
	  mcp_ptr->PdgCode() == 11 ||
	  mcp_ptr->PdgCode() == -11 ) {
	
	shower_parent_id.insert(mcp_ptr->TrackId());

      }

    }
    _shower_id.clear();
    _shower_id.resize(_start_vtx.size(),-1);

    if(!_mother.size()) return;

    // Construct granular showers (will combine later)

    // Commented out: this version combines all e-/e+/gamma daughters
    // Problem with this is that it does not differentiate 2 gammas
    // coming from the same pi0 (as they share same parent ID).
    // We define shower "mother" has to be e+/e-/gamma, so the
    // parent must have pdgcode 11/-11/22
    /*
    for(size_t i=0; i<_mother.size(); ++i) {

      unsigned int track_id   = _track_id.at(i);
      unsigned int parent_id  = _mother.at(i);
      double       time       = _start_vtx.at(i).at(3);

      while(1) {

	if( shower_parent_id.find(parent_id) == shower_parent_id.end() )
	  
	  break;
	  
	else {

	  auto grandma_iter = _track_index.find(parent_id);

	  if(grandma_iter == _track_index.end()) break;

	  parent_id = (*grandma_iter).second;

	}

      }

      auto shower_index_iter = _shower_index.find(i);
      size_t shower_index = 0;
      if(shower_index_iter == _shower_index.end()) {

	shower_index = _shower_index.size();
	_shower_index.insert(std::pair<unsigned int, unsigned int>(i,shower_index));
	_shower_daughters.push_back(std::map<double,unsigned int>());
      }
      else shower_index = (*shower_index_iter).second;
      
      _shower_daughters[shower_index].insert(std::pair<double,unsigned int>(time,i));
      _shower_id[i] = shower_index;
      // In case this particle has daughters, add them
      for(auto const daughter_track : _daughters.at(i)) {

	auto const daughter_index_iter = _track_index.find(daughter_track);
	if( daughter_index_iter == _track_index.end() ) continue;
	
	_shower_daughters[shower_index].insert(std::pair<double,unsigned int>(_start_vtx.at((*daughter_index_iter).second).at(3),
	                                                                      (*daughter_index_iter).second));
	_shower_id[i] = shower_index;
      }

    }
    */
    
    std::vector<std::map<Double_t,std::set<unsigned int> > > ordered_shower_daughters;

    for(size_t i=0; i<_mother.size(); ++i) {
      
      unsigned int parent_id  = 0;
      unsigned int grandma_id = _mother.at(i);
      double       time       = _start_vtx.at(i).at(3);

      while(1) {

	if( shower_parent_id.find(grandma_id) == shower_parent_id.end() )

	  break;
	  
	else {

	  parent_id = grandma_id;

	  auto grandma_iter = _track_index.find(grandma_id);

	  if(grandma_iter == _track_index.end()) break;

	  grandma_id = _mother.at((*grandma_iter).second);

	}

      }

      if(!parent_id) continue;
      unsigned int parent_index = (*(_track_index.find(parent_id))).second;

      auto shower_index_iter = _shower_index.find(parent_index);
      size_t shower_index = 0;
      if(shower_index_iter == _shower_index.end()) {

	shower_index = _shower_index.size();
	_shower_index.insert(std::pair<unsigned int, unsigned int>(parent_index,shower_index));
	ordered_shower_daughters.push_back(std::map<double,std::set<unsigned int> >());
      }
      else shower_index = (*shower_index_iter).second;

      // Add to a shower daughter... if there's already a particle @ same time T, shift by 1e-12
      if(ordered_shower_daughters.at(shower_index).find(time) == ordered_shower_daughters.at(shower_index).end() )

	ordered_shower_daughters[shower_index].insert(std::pair<double,std::set<unsigned int> >(time,std::set<unsigned int>()));

      ordered_shower_daughters[shower_index][time].insert(i);
      _shower_id[i]=shower_index;

      // In case this particle has daughters, add them
      for(auto const daughter_track : _daughters.at(i)) {

	auto const daughter_index_iter = _track_index.find(daughter_track);
	if( daughter_index_iter == _track_index.end() ) continue;
	
	ordered_shower_daughters[shower_index][time].insert((*daughter_index_iter).second);
	_shower_id[(*daughter_index_iter).second]=shower_index;
      }
    }

    // Store ordered daughters' list
    _shower_daughters.clear();
    _shower_daughters.reserve(ordered_shower_daughters.size());
    for(auto time_daughters : ordered_shower_daughters) {
      _shower_daughters.push_back(std::vector<unsigned int>());
      auto shower_daughters = _shower_daughters.rbegin();
      shower_daughters->reserve(time_daughters.size());
      std::set<unsigned int> unique_daughter_set;
      for(auto const time_daughter_pair : time_daughters){
	for(auto const daughter_index : time_daughter_pair.second){
	  if(unique_daughter_set.find(daughter_index)==unique_daughter_set.end()) {
	    unique_daughter_set.insert(daughter_index);
	    shower_daughters->push_back(daughter_index);
	  }
	}
      }
    }

    if(_debug_mode) {

      std::cout << std::endl << Form("Found %zu granular showers...",_shower_index.size())<<std::endl;
      
      for(auto mother_iter = _shower_index.begin();
	  mother_iter != _shower_index.end();
	  ++mother_iter) {
	
	unsigned int mother_index = (*mother_iter).first;
	unsigned int shower_index = (*mother_iter).second;
	unsigned int mother_track = _track_id.at(mother_index);
	unsigned int ndaughters = _shower_daughters.at(shower_index).size();
	double x,y,z,t;
	GetTrackStartInfo(mother_index,x,y,z,t);

	std::cout<<Form("Shower track ID = %d, PDG=%d,  @(%g, %g, %g, %g) ...  %d daughters ",
			mother_track,
			_pdgcode.at(mother_index),
			x,y,z,t,
			ndaughters)<<std::endl;
	
	unsigned int first_daughter_index = (*(_shower_daughters.at(shower_index).begin()));
	unsigned int daughter_pdgcode = _pdgcode.at(first_daughter_index);
	
	GetTrackStartInfo(first_daughter_index,x,y,z,t);
	std::cout << Form("  Daughter %d starting @ (%g, %g, %g, %g) =>",daughter_pdgcode,x,y,z,t);
	GetTrackEndInfo(first_daughter_index,x,y,z,t);
	std::cout << Form(" (%g, %g, %g, %g)",x,y,z,t)<<std::endl;
	
      }
      std::cout<<std::endl<<"End ConstructGranularShower ..."<<std::endl<<std::endl;
    }
  }

  void UBMCShowerReco::GetTrackInfo(const unsigned int &index,
				    double &start_x,
				    double &start_y,
				    double &start_z,
				    double &start_time,
				    double &end_x,
				    double &end_y,
				    double &end_z,
				    double &end_time)
  {

    if( index > _track_id.size() )

      throw cet::exception(__FUNCTION__) << Form("Particle index %d not found!",index);
    
    start_x = _start_vtx.at(index).at(0);
    start_y = _start_vtx.at(index).at(1);
    start_z = _start_vtx.at(index).at(2);
    start_time = _start_vtx.at(index).at(3);

    end_x = _end_vtx.at(index).at(0);
    end_y = _end_vtx.at(index).at(1);
    end_z = _end_vtx.at(index).at(2);
    end_time = _end_vtx.at(index).at(3);

  }

  void UBMCShowerReco::GetTrackStartInfo(const unsigned int &index,
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

  void UBMCShowerReco::GetTrackEndInfo(const unsigned int &index,
				      double &end_x,
				      double &end_y,
				      double &end_z,
				      double &end_time)
  {
    if(index > _track_id.size())

      throw cet::exception(__FUNCTION__) << Form("Particle index %d not found!",index);

    end_x = _end_vtx.at(index).at(0);
    end_y = _end_vtx.at(index).at(1);
    end_z = _end_vtx.at(index).at(2);
    end_time = _end_vtx.at(index).at(3);

  }

  void UBMCShowerReco::CombineGranularShower()
  {
    // Let's get start & end track ID of each shower
    std::vector<unsigned int> track_start_v;
    track_start_v.reserve(_shower_daughters.size());
    std::vector<unsigned int> track_end_v;
    track_end_v.reserve(_shower_daughters.size());

    for(auto shower_iter = _shower_index.begin();
	shower_iter != _shower_index.end();
	++shower_iter) {
      
      unsigned int shower_index  = (*shower_iter).second;
      unsigned int primary_index = (*shower_iter).first;

      track_start_v.push_back(primary_index);
      track_end_v.push_back  ((*(_shower_daughters.at(shower_index).rbegin())));

      if(_debug_mode) std::cout << Form("Shower Track ID %d => %d", 
					_track_id.at(primary_index),
					track_end_v.at(track_end_v.size()-1))
				<<std::endl;
    }
    
    std::vector<std::set<unsigned int> > combination_v;
    combination_v.reserve(track_start_v.size());
    
    for(size_t i=0; i<track_start_v.size(); ++i) {

      std::set<unsigned int> combination;
      double tstart_i, xstart_i, ystart_i, zstart_i;
      GetTrackStartInfo(track_start_v.at(i),
			xstart_i,ystart_i,zstart_i,tstart_i);

      for(size_t j=0; j<track_start_v.size(); ++j) {
	
	// Check if shower "i" belongs to "j".
	// Skip if "i"=="j"
	if(i==j) continue;

	// Skip if the combination (i,j) is already inspected and found to be combined
	if(combination_v.size()>j && combination_v.at(j).find(i) != combination_v.at(j).end())
	  continue;

	// Skip if "j" start time is not inside "i" shower time, which is 1st daughter to last daughter time
	double tstart_j, xstart_j, ystart_j, zstart_j;
	double tend_j, xend_j, yend_j, zend_j;
	GetTrackStartInfo( (*(_shower_daughters.at(j).begin())),
			   xstart_j,ystart_j,zstart_j,tstart_j);
	GetTrackStartInfo( (*(_shower_daughters.at(j).rbegin())),
			   xend_j,yend_j,zend_j,tend_j);
	
	if( tstart_j > (tstart_i + _comb_time_cut) || (tend_j + _comb_time_cut) < tstart_i )
	  continue;

	if(_debug_mode) std::cout<<Form("Shower %zu starts within %zu! Inspecting for merging...",i,j)<<std::endl;

	bool combine = false;
	double dT = 0;
	double dL2 = 0;
	if( tstart_i < tstart_j ){
	// Case 1: shower "i" is just before shower "j" ... use 1st daughter
	  dT = tstart_i - tstart_j;
	  double x,y,z,t;
	  auto const daughter_iter = (_shower_daughters.at(j).begin());
	  GetTrackStartInfo((*daughter_iter), x, y, z, t);
	  dL2 = (pow(xstart_i - x,2) + pow(ystart_i - y,2) + pow(zstart_i - z,2));
	  combine = dL2 < _comb_dist2_cut;
	}else if( tend_j < tstart_i){
	// Case 2: shower "i" is just after shower "j"
	  dT = tstart_i - tend_j;
	  double x,y,z,t;
	  auto const daughter_iter = (_shower_daughters.at(j).rbegin());
	  GetTrackStartInfo((*daughter_iter), x, y, z, t);
	  dL2 = (pow(xstart_i - x,2) + pow(ystart_i - y,2) + pow(zstart_i - z,2));
	  combine = dL2 < _comb_dist2_cut;
	}else{
	// Case 3: shower "i" is within shower "j"

	// Loop over "j"'s daughters and see if i's mother is generated nearby
	  double x,y,z,t;
	  for(auto daughter_iter = _shower_daughters.at(j).begin();
	      daughter_iter != _shower_daughters.at(j).end();
	      ++daughter_iter) {

	    double       step_time  = _start_vtx.at((*daughter_iter)).at(3);
	    unsigned int step_index = (*daughter_iter);
	    
	    if( (step_time < tstart_i) && (tstart_i - step_time) > _comb_time_cut ) continue;
	    if( (step_time > tstart_i) && (step_time - tstart_i) > _comb_time_cut ) break;

	    if(_debug_mode) std::cout<<Form(" Inspecting @ %d...",_track_id.at(step_index))<<std::endl;
	    
	    GetTrackStartInfo(step_index,x,y,z,t);
	    double dist2 = (pow(xstart_i - x,2) + pow(ystart_i - y,2) + pow(zstart_i - z,2));
	    combine = dist2 < _comb_dist2_cut;
	    if(combine) { dL2 = dist2; dT = tstart_i - step_time; break; }
	  }

	}

	if(combine) { 
	 
	  if(_debug_mode) std::cout<<Form("    Found merging point! dT = %g, dX^2 = %g", dT,dL2)<<std::endl;

	  combination.insert(j);
	}
      }
      combination_v.push_back(combination);
    }

    // Now make a combined shower index
    std::vector<std::set<unsigned int> > sorted_combination_v;

    for(size_t i=0; i<combination_v.size(); ++i) {

      int sorted_index = -1;

      // First, attempt to find "i" and associated indexes in sorted_combination_v contents.
      for(size_t j=0; j<sorted_combination_v.size(); ++j) {

	if(sorted_combination_v.at(j).find(i) != sorted_combination_v.at(j).end()) {
	  sorted_index = j; 
	  break;
	}
	for(auto const index : combination_v.at(i)) {

	  if(sorted_combination_v.at(j).find(index) != sorted_combination_v.at(j).end()) {
	    sorted_index = j;
	    break;
	  }
	}
	if(sorted_index>0) break;
      }
      // If not found, create a new set and insert
      if(sorted_index<0) {
	std::set<unsigned int> sorted_combination(combination_v.at(i));
	sorted_combination.insert(i);
	sorted_combination_v.push_back(sorted_combination);
	
      }
      // Else combine sets
      else{

	sorted_combination_v[sorted_index].insert(i);
	for(auto const index : combination_v.at(i))

	  sorted_combination_v[sorted_index].insert(index);
      }
    }
    
    
    std::map<unsigned int,unsigned int> combined_shower_index;
    std::vector<std::map<double,std::vector<unsigned int> > > combined_shower_daughters;
    for(size_t i=0; i<sorted_combination_v.size(); ++i) {

      int super_mother_index=-1;
      double super_mother_time=-1;
      std::map<double,std::vector<unsigned int> > daughters;
      for(auto const index : sorted_combination_v.at(i)) {

	// Combine daughters
	for(auto const daughter_index : _shower_daughters.at(index)) {
	  
	  double daughter_time = _start_vtx.at(daughter_index).at(3);
	  if(daughters.find(daughter_time)==daughters.end())
	    daughters.insert(std::pair<double,std::vector<unsigned int> >(daughter_time,
									  std::vector<unsigned int>(1,
												    daughter_index)
									  )
			     );
	  else
	    daughters[daughter_time].push_back(daughter_index);
	  
	}

	// Find super mother
	int mother_index = track_start_v.at(index);
	double mother_time = _start_vtx.at(mother_index).at(3);
	if(super_mother_index<0 || mother_time < super_mother_time) {
	  super_mother_index=mother_index;
	  super_mother_time=mother_time;
	}else if(mother_time == super_mother_time && mother_index < super_mother_index) {
	  super_mother_index = mother_index;
	  super_mother_time=mother_time;
	}
      }
      // Add mothers that is not super mother, to a daughter
      for(auto const index : sorted_combination_v.at(i)) {
	int mother_index = track_start_v.at(index);
	double mother_time = _start_vtx.at(mother_index).at(3);
	if(mother_index!=super_mother_index) {
	  if(daughters.find(mother_time)==daughters.end())
	    daughters.insert(std::pair<double,std::vector<unsigned int> >(mother_time,
									  std::vector<unsigned int>(1,
												    mother_index)
									  )
			     );
	  else
	    daughters[mother_time].push_back(mother_index);
	}
	else
	  combined_shower_index.insert(std::pair<unsigned int,unsigned int>(super_mother_index,i));
      }
      combined_shower_daughters.push_back(daughters);
    }

    // Re-set showers
    _shower_index = combined_shower_index; 

    // Update shower id
    _shower_daughters.clear();
    _shower_daughters.reserve(combined_shower_daughters.size());
    for(size_t i=0; i<_shower_id.size(); ++i)
      _shower_id[i]=-1;

    for(auto const daughters_map : combined_shower_daughters) {
      
      unsigned int shower_index = _shower_daughters.size();
      std::vector<unsigned int> daughters_v;
      daughters_v.reserve(daughters_map.size());
      for(auto daughters_map_iter = daughters_map.begin();
	  daughters_map_iter != daughters_map.end();
	  ++daughters_map_iter) {

	for(auto const index : (*daughters_map_iter).second) {
	  _shower_id[index] = shower_index;
	  daughters_v.push_back(index);
	}
      }
      _shower_daughters.push_back(daughters_v);
    }
    for(auto mother_iter = _shower_index.begin();
	mother_iter != _shower_index.end();
	++mother_iter)

      _shower_id[(*mother_iter).first] = (*mother_iter).second;
    
    if(_debug_mode) {
      std::vector<unsigned int> part_count_v(_shower_index.size(),0);
      unsigned int undefined=0;
      for(size_t i=0; i<_shower_id.size(); ++i) {
	if(_shower_id.at(i)<0) {
	  undefined++;
	  std::cout << Form("  Track %d (PDG=%d) @ (%g, %g, %g, %g) with %zu daughters does not belong to a shower!",
			    _track_id[i],
			    _pdgcode.at(i),
			    _start_vtx.at(i).at(0), _start_vtx.at(i).at(1), _start_vtx.at(i).at(2), _start_vtx.at(i).at(3),
			    _daughters.at(i).size()) 
		    << std::endl;
	}else if((size_t)(_shower_id.at(i))>=part_count_v.size())
	  
	  throw cet::exception(__FUNCTION__) << Form("Track %d PDG %d has ill-defined shower index %d!",
						     _track_id.at(i),_pdgcode.at(i),_shower_id.at(i));
	else
	  part_count_v[_shower_id.at(i)]++;
      }
      std::cout << Form("  %d tracks do not belong to shower...",undefined) << std::endl;
      
      for(size_t i=0; i<part_count_v.size(); ++i)
	
	std::cout << Form("  %d tracks belong to shower %zu...", part_count_v.at(i),i)<<std::endl;
      
    }
    
  } // namespace opdet
}
//  LocalWords:  ifndef
