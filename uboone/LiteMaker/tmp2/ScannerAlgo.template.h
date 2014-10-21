#ifndef SCANNERALGO_TEMPLATE_H
#define SCANNERALGO_TEMPLATE_H

/*
  This file defines certain specilization of templated functions.
  In particular it implements:

  ScannerAlgo::ScanData
  ScannerAlgo::GetPtrMap
  ScannerAlgo::LiteDataType

  One has to implement the relevant specilization for introducing a new data product!

 */

namespace larlite {

  //
  // ScanData definition
  //
  template <> void ScannerAlgo::ScanData(art::Handle<std::vector< ::simb::MCTruth> > const &dh,
					 ::larlite::event_base* lite_dh);
  template <> void ScannerAlgo::ScanData(art::Handle<std::vector< ::simb::GTruth> > const &dh,
					 ::larlite::event_base* lite_dh);
  template <> void ScannerAlgo::ScanData(art::Handle<std::vector< ::simb::MCParticle> > const &dh,
					 ::larlite::event_base* lite_dh);
  template <> void ScannerAlgo::ScanData(art::Handle<std::vector< ::simb::MCFlux> > const &dh,
					 ::larlite::event_base* lite_dh);
  template <> void ScannerAlgo::ScanData(art::Handle<std::vector< ::sim::MCShower> > const &dh,
					 ::larlite::event_base* lite_dh);
  template <> void ScannerAlgo::ScanData(art::Handle<std::vector< ::sim::SimChannel> > const &dh,
					 ::larlite::event_base* lite_dh);
  template <> void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::Wire> > const &dh,
					 ::larlite::event_base* lite_dh);
  template <> void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::Hit> > const &dh,
					 ::larlite::event_base* lite_dh);
  template <> void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::OpHit> > const &dh,
					 ::larlite::event_base* lite_dh);
  template <> void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::OpFlash> > const &dh,
					 ::larlite::event_base* lite_dh);
  template <> void ScannerAlgo::ScanData(art::Handle<std::vector< ::anab::CosmicTag> > const &dh,
					 ::larlite::event_base* lite_dh);
  template <> void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::Cluster> > const &dh,
					 ::larlite::event_base* lite_dh);
  template <> void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::Seed> > const &dh,
					 ::larlite::event_base* lite_dh);
  template <> void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::EndPoint2D> > const &dh,
					 ::larlite::event_base* lite_dh);
  template <> void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::SpacePoint> > const &dh,
					 ::larlite::event_base* lite_dh);
  template <> void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::Track> > const &dh,
					 ::larlite::event_base* lite_dh);
  template <> void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::Shower> > const &dh,
					 ::larlite::event_base* lite_dh);
  template <> void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::Vertex> > const &dh,
					 ::larlite::event_base* lite_dh);
  template <> void ScannerAlgo::ScanData(art::Handle<std::vector< ::anab::Calorimetry> > const &dh,
					 ::larlite::event_base* lite_dh);
  template <> void ScannerAlgo::ScanData(art::Handle<std::vector< ::anab::ParticleID> > const &dh,
					 ::larlite::event_base* lite_dh);
  template <> void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::PFParticle> > const &dh,
					 ::larlite::event_base* lite_dh);
  template <class T>
  void ScanData(art::Handle<std::vector<T> > const &dh,
		::larlite::event_base* lite_dh)
  { throw cet::exception(__PRETTY_FUNCTION__)<<"Not implemented!"; }
  
  //
  // Getter for associated data product pointer 
  //
  template <> const std::map<art::Ptr< ::simb::MCTruth>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const;
  template <> const std::map<art::Ptr< ::simb::GTruth>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const;
  template <> const std::map<art::Ptr< ::simb::MCFlux>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const;
  template <> const std::map<art::Ptr< ::simb::MCParticle>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const;
  template <> const std::map<art::Ptr< ::sim::SimChannel>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const;
  template <> const std::map<art::Ptr< ::sim::MCShower>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const;
  template <> const std::map<art::Ptr< ::recob::OpHit>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const;
  template <> const std::map<art::Ptr< ::recob::OpFlash>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const;
  template <> const std::map<art::Ptr< ::recob::Hit>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const;
  template <> const std::map<art::Ptr< ::recob::Wire>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const;
  template <> const std::map<art::Ptr< ::recob::Cluster>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const;
  template <> const std::map<art::Ptr< ::recob::Track>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const;
  template <> const std::map<art::Ptr< ::recob::Shower>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const;
  template <> const std::map<art::Ptr< ::recob::Vertex>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const;
  template <> const std::map<art::Ptr< ::recob::SpacePoint>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const;
  template <> const std::map<art::Ptr< ::recob::EndPoint2D>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const;
  template <> const std::map<art::Ptr< ::anab::CosmicTag>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const;
  template <> const std::map<art::Ptr< ::anab::Calorimetry>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const;
  template <> const std::map<art::Ptr< ::anab::ParticleID>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const;
  template <> const std::map<art::Ptr< ::recob::PFParticle>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const;

  template <class T>
  const std::map<art::Ptr<T>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap() const
  { throw cet::exception(__PRETTY_FUNCTION__)<<"Not implemented for a specified data product type..."; }

  //
  // Type identifier functions
  //
  // simb
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::simb::GTruth> () const;
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::simb::MCTruth> () const;
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::simb::MCParticle> () const;
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::simb::MCFlux> () const;
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::simb::MCTrajectory> () const;
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::simb::MCNeutrino> () const;
  // sim
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::sim::SimChannel> () const;
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::sim::MCShower> () const;
  // raw
  // recob
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::Wire> () const;
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::Hit> () const;
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::Cluster> () const;
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::SpacePoint> () const;
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::OpHit> () const;
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::OpFlash> () const;
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::Seed> () const;
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::Track> () const;
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::Shower> () const;
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::Vertex> () const;
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::EndPoint2D> () const;
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::PFParticle> () const;
  // anab
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::anab::CosmicTag> () const;
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::anab::Calorimetry> () const;
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::anab::ParticleID> () const;

  //
  // LocateLiteProduct implementation
  //
  template <class T>
  bool ScannerAlgo::LocateLiteProduct(art::Ptr<T> const ptr,
				      std::pair<size_t,size_t> &loc) const
  { 
    auto ptr_map = GetPtrMap<T>();
    auto id_iter = ptr_map.find(ptr);
    if(id_iter == ptr_map.end()) return false; 
    
    loc.first = (*id_iter).second.first;
    loc.second = (*id_iter).second.second;
    return true;
  }
  
  //
  // ScanAssociation implementation 
  //
  template <class T,class U>
  void ScannerAlgo::ScanAssociation(art::Event const& e,
				    art::Handle<std::vector<T> > &dh,
				    ::larlite::event_base* lite_dh) const
  { 

    art::FindManyP<U> ptr_coll_v(dh, e, lite_dh->name());

    auto ass_type = LiteDataType<U>();

    // Instantiate association container. length = # of producers for associated data type
    std::vector< ::larlite::AssSet_t > ass_set_v(fModuleLabel_v[(size_t)(ass_type)].size());

    // Return if there's no data products stored for associated data type
    if(!(ass_set_v.size())) return;

    std::pair<size_t,size_t> lite_location;

    // Loop over origin data product vector, find associated objects and store association info
    for(size_t i=0; i<dh->size(); ++i) {

      const std::vector<art::Ptr<U> > ptr_coll = ptr_coll_v.at(i);

      // Association vector: one per associated data product producers
      std::vector<larlite::AssUnit_t> ass_unit_v(ass_set_v.size());

      for(auto& art_ptr : ptr_coll) {

	if(!LocateLiteProduct(art_ptr,lite_location)) continue;

	ass_unit_v[lite_location.first].push_back(lite_location.second);
      }
      for(size_t i=0; i<ass_set_v.size(); ++i)

	ass_set_v[i].push_back(ass_unit_v[i]);
      
    } // end looping over origin data products
	
    // Store associations in larlite data products
    for(size_t i=0; i<ass_set_v.size(); ++i) {

      // Loop over in one association set, store if there's any
      for(auto const& ass_unit : ass_set_v[i]) {

	if(ass_unit.size()) {

	  auto ass_name = fModuleLabel_v[(size_t)ass_type][i];

	  larlite::product_id ass_id(ass_type,ass_name);
	  
	  lite_dh->set_association(ass_id,ass_set_v[i]);

	  break;
	}
      }// end looping over association set
    }// end looping over a vector of association set
  }

  template <> void ScannerAlgo::ScanAssociation <::recob::Cluster,::recob::Cluster>(art::Event const& e,
										art::Handle<std::vector<::recob::Cluster> > &dh,
										::larlite::event_base* lite_dh) const
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::recob::EndPoint2D,::recob::EndPoint2D>(art::Event const& e,
										      art::Handle<std::vector<::recob::EndPoint2D> > &dh,
										      ::larlite::event_base* lite_dh) const
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::recob::Vertex,::recob::Vertex>(art::Event const& e,
										      art::Handle<std::vector<::recob::Vertex> > &dh,
										      ::larlite::event_base* lite_dh) const
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::recob::SpacePoint,::recob::SpacePoint>(art::Event const& e,
										      art::Handle<std::vector<::recob::SpacePoint> > &dh,
										      ::larlite::event_base* lite_dh) const
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::recob::Track,::recob::Track>(art::Event const& e,
									    art::Handle<std::vector<::recob::Track> > &dh,
									    ::larlite::event_base* lite_dh) const
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::recob::Shower,::recob::Shower>(art::Event const& e,
									      art::Handle<std::vector<::recob::Shower> > &dh,
									      ::larlite::event_base* lite_dh) const
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::anab::Calorimetry,::anab::Calorimetry>(art::Event const& e,
										      art::Handle<std::vector<::anab::Calorimetry> > &dh,
										      ::larlite::event_base* lite_dh) const
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::anab::ParticleID,::anab::ParticleID>(art::Event const& e,
										    art::Handle<std::vector<::anab::ParticleID> > &dh,
										    ::larlite::event_base* lite_dh) const
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::simb::MCTruth,::simb::MCTruth>(art::Event const& e,
									      art::Handle<std::vector<::simb::MCTruth> > &dh,
									      ::larlite::event_base* lite_dh) const
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::simb::MCParticle,::simb::MCParticle>(art::Event const& e,
										    art::Handle<std::vector<::simb::MCParticle> > &dh,
										    ::larlite::event_base* lite_dh) const
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::recob::PFParticle,::recob::PFParticle>(art::Event const& e,
										      art::Handle<std::vector<::recob::PFParticle> > &dh,
										      ::larlite::event_base* lite_dh) const
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  //
  // LiteDataType
  //
  template <class T>
  const ::larlite::data::DataType_t ScannerAlgo::LiteDataType() const
  { throw cet::exception(__PRETTY_FUNCTION__)<<"Unsupported data type conversion attempted!";
    return ::larlite::data::kUndefined;
  }

  //
  // ProductID
  //
  template <class T>
  const ::larlite::product_id ScannerAlgo::ProductID(size_t name_index) const
  { auto data_type = LiteDataType<T>();
    if(fModuleLabel_v[(size_t)data_type].size() <= name_index)
      throw cet::exception(__PRETTY_FUNCTION__)
	<< "Length of registered products for data type " << ::larlite::data::kDATA_TREE_NAME[data_type].c_str()
	<< " is " << fModuleLabel_v[(size_t)data_type].size()
	<< " while you requested " << name_index;
    return ::larlite::product_id(data_type,fModuleLabel_v[(size_t)(data_type)][name_index]);
  }
}
#endif
