
#ifndef SCANNERALGO_H
#define SCANNERALGO_H

#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArLite include
#include "DataFormat/storage_manager.h"

// LArSoft includes
#include "uboone/MuCS/MuCSData.h"
#include "uboone/MuCS/MuCSRecoData.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "larsimobj/Simulation/SimChannel.h"
#include "larsimobj/Simulation/SimPhotons.h"
#include "larsimobj/Simulation/AuxDetSimChannel.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/OpticalDetectorData/FIFOChannel.h"
#include "lardataobj/OpticalDetectorData/OpticalTypes.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// std 
#include <vector>
#include <string>
#include <iostream>


#include <map>

namespace larlite {

  /**
     \class ScannerAlgo
     Algorithm class to scan LArSoft data products into LArLite data products.
     This class itself does not own any data product as it's just an algorithm.
     The design is to receive data container art::Handle<std::vector<T>> (LArSoft
     data product container) and resulting data container in LArLite and fill the
     latter.
     It also supports storing of associations.
   */
  class ScannerAlgo {

  public:
    /// default ctor
    ScannerAlgo() 
      : fModuleLabel_v    ((size_t)(::larlite::data::kDATA_TYPE_MAX),std::vector<std::string>())
      , fAssModuleLabel_v ((size_t)(::larlite::data::kDATA_TYPE_MAX),std::vector<std::string>())
      , fDataReadFlag_v   ((size_t)(::larlite::data::kDATA_TYPE_MAX),std::map<std::string,bool>())
    {}
    /// default dtor
    ~ScannerAlgo(){}

    /**
       Function to register a producer for a specified data product.
       This is used to internally keep track of recorded data products for association
       and also generating LArLite data product ID when requested.
    */
    void Register(std::string const& name,
		  ::larlite::data::DataType_t const data_type)
    {
      bool exist=false;
      for(auto const& label : fModuleLabel_v[data_type]) if(label==name) { exist=true; break;}
      if(!exist) {
	fModuleLabel_v[data_type].push_back(name);
	fDataReadFlag_v[data_type].insert(std::make_pair(name,false));
      }
      AssociationRegister(name,data_type);
    }

    /**
       Function to register a producer for associated data products.
    */
    void AssociationRegister(std::string const& name,
			     ::larlite::data::DataType_t const data_type)
    {
      bool exist=false;
      for(auto const& label : fAssModuleLabel_v[data_type]) if(label==name) { exist=true; break;}
      if(!exist) fAssModuleLabel_v[data_type].push_back(name);
    }

    /// Accessor to the list of registered producer module labels
    std::vector<std::vector<std::string> > const& ModuleLabels() const { return fModuleLabel_v; }

    /// Accessor to the list of registered associated products' producers' module labels
    std::vector<std::vector<std::string> > const& AssLabels() const { return fAssModuleLabel_v; }

    /// Function to be called @ end or beginning of each event
    void EventClear();

    /// Method to define a integer key per data product
    template <class T>
    void ProducePtrMapKey(const art::Ptr<T>& ptr, size_t& key1, size_t& key2);

    /// Method to generate a data product ID for a specified type (through template) and producer module label's index
    template <class T>
    const ::larlite::product_id ProductID(size_t name_index) const;

    /// Method to generate a association data product ID for a specified type (through template) and producer module label's index
    template <class T>
    const ::larlite::product_id AssProductID(size_t name_index) const;

    /// Core method: convert LArSoft data product (dh) to LArLite (lite_dh)
    template <class T>
    void ScanData(art::Handle<std::vector<T> > const &dh,
		  ::larlite::event_base* lite_dh);
    
    /// Core method: generate LArLite association data product and store (in lite_dh)
    template <class T, class U>
    void ScanAssociation(art::Event const& e,
			 art::Handle<std::vector<T> > &dh,
			 ::larlite::event_ass* lite_dh);

    /// Accessor to art::Ptr map ... used to locate associated data product location
    template <class T>
    std::map<art::Ptr<T>,std::pair<size_t,size_t> >& GetPtrMap(size_t key1=0, size_t key2=0);

    /// Utility function to link LArSoft data product type to LArLite enum
    template <class T>
    const ::larlite::data::DataType_t LiteDataType() const;

  private:
    
    /// Internal utility function to look up associated object locater in LArLite data holder
    template <class T>
    bool LocateLiteProduct(art::Ptr<T> const ptr,
			   std::pair<size_t,size_t> &loc);

    /// Internal utility function to look up the index number of given producer's name (for a specified type)
    size_t NameIndex(::larlite::data::DataType_t const data_type,
		     std::string const& name) const;

    /// Holder for producer module labels. Outer vector index = larlite::data::DataType_t
    std::vector<std::vector<std::string> > fModuleLabel_v;

    /// Holder for associated data producers' label. Outer vector index = larlite::data::DataType_t
    std::vector<std::vector<std::string> > fAssModuleLabel_v;

    /// Boolean holder to tell us whether specific producer's data is read or not
    std::vector<std::map<std::string,bool> > fDataReadFlag_v;

    // art::Ptr local storage. Value = index of data product & index of label
    std::vector< std::vector< std::map< art::Ptr<::simb::MCTruth>,     std::pair<size_t,size_t> > > > fPtrIndex_mctruth;
    std::vector< std::vector< std::map< art::Ptr<::simb::GTruth>,      std::pair<size_t,size_t> > > > fPtrIndex_gtruth;
    std::vector< std::vector< std::map< art::Ptr<::simb::MCFlux>,      std::pair<size_t,size_t> > > > fPtrIndex_mcflux;
    std::vector< std::vector< std::map< art::Ptr<::simb::MCParticle>,  std::pair<size_t,size_t> > > > fPtrIndex_mcpart;
    std::vector< std::vector< std::map< art::Ptr<::sim::SimChannel>,   std::pair<size_t,size_t> > > > fPtrIndex_simch;
    std::vector< std::vector< std::map< art::Ptr<::sim::AuxDetSimChannel>,   std::pair<size_t,size_t> > > > fPtrIndex_auxsimch;
    std::vector< std::vector< std::map< art::Ptr<::sim::MCShower>,     std::pair<size_t,size_t> > > > fPtrIndex_mcshower;
    std::vector< std::vector< std::map< art::Ptr<::sim::MCTrack>,      std::pair<size_t,size_t> > > > fPtrIndex_mctrack;
    std::vector< std::vector< std::map< art::Ptr<::raw::RawDigit>,     std::pair<size_t,size_t> > > > fPtrIndex_rawdigit;
    std::vector< std::vector< std::map< art::Ptr<::raw::OpDetWaveform>,std::pair<size_t,size_t> > > > fPtrIndex_opdigit;
    std::vector< std::vector< std::map< art::Ptr<::raw::Trigger>,      std::pair<size_t,size_t> > > > fPtrIndex_trigger;
    std::vector< std::vector< std::map< art::Ptr<::recob::Wire>,       std::pair<size_t,size_t> > > > fPtrIndex_wire;
    std::vector< std::vector< std::map< art::Ptr<::recob::Hit>,        std::pair<size_t,size_t> > > > fPtrIndex_hit;
    std::vector< std::vector< std::map< art::Ptr<::recob::OpHit>,      std::pair<size_t,size_t> > > > fPtrIndex_ophit;
    std::vector< std::vector< std::map< art::Ptr<::recob::OpFlash>,    std::pair<size_t,size_t> > > > fPtrIndex_opflash;
    std::vector< std::vector< std::map< art::Ptr<::recob::Cluster>,    std::pair<size_t,size_t> > > > fPtrIndex_cluster;
    std::vector< std::vector< std::map< art::Ptr<::recob::Shower>,     std::pair<size_t,size_t> > > > fPtrIndex_shower;
    std::vector< std::vector< std::map< art::Ptr<::recob::Vertex>,     std::pair<size_t,size_t> > > > fPtrIndex_vertex;
    std::vector< std::vector< std::map< art::Ptr<::recob::Track>,      std::pair<size_t,size_t> > > > fPtrIndex_track;
    std::vector< std::vector< std::map< art::Ptr<::anab::CosmicTag>,   std::pair<size_t,size_t> > > > fPtrIndex_cosmictag;
    std::vector< std::vector< std::map< art::Ptr<::anab::Calorimetry>, std::pair<size_t,size_t> > > > fPtrIndex_calo;
    std::vector< std::vector< std::map< art::Ptr<::recob::SpacePoint>, std::pair<size_t,size_t> > > > fPtrIndex_sps;
    std::vector< std::vector< std::map< art::Ptr<::recob::EndPoint2D>, std::pair<size_t,size_t> > > > fPtrIndex_end2d;
    std::vector< std::vector< std::map< art::Ptr<::recob::Seed>,       std::pair<size_t,size_t> > > > fPtrIndex_seed;
    std::vector< std::vector< std::map< art::Ptr<::anab::ParticleID>,  std::pair<size_t,size_t> > > > fPtrIndex_partid;
    std::vector< std::vector< std::map< art::Ptr<::recob::PFParticle>, std::pair<size_t,size_t> > > > fPtrIndex_pfpart;
    std::vector< std::vector< std::map< art::Ptr<::recob::PCAxis>,     std::pair<size_t,size_t> > > > fPtrIndex_pcaxis;
    std::vector< std::vector< std::map< art::Ptr<::anab::FlashMatch>,  std::pair<size_t,size_t> > > > fPtrIndex_fmatch;
  };
}

#include "ScannerAlgo.template.h"

#endif
