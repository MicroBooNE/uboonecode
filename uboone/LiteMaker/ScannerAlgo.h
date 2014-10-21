
#ifndef SCANNERALGO_H
#define SCANNERALGO_H

#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArLite include
#include "DataFormat/storage_manager.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "RawData/RawDigit.h"
#include "RecoBase/Wire.h"
#include "RecoBase/Hit.h"
#include "RecoBase/OpHit.h"
#include "RecoBase/OpFlash.h"
#include "RecoBase/Track.h"
#include "RecoBase/Seed.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Shower.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/EndPoint2D.h"
#include "RecoBase/PFParticle.h"
#include "AnalysisBase/ParticleID.h"
#include "AnalysisBase/Calorimetry.h"
#include "AnalysisBase/CosmicTag.h"
#include "Simulation/SimChannel.h"
#include "SimulationBase/MCFlux.h"
#include "SimulationBase/GTruth.h"
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"
#include "OpticalDetectorData/FIFOChannel.h"
#include "OpticalDetectorData/OpticalTypes.h"
#include "MCBase/MCShower.h"
#include "Utilities/LArProperties.h"
#include "Utilities/GeometryUtilities.h"
#include "Utilities/DetectorProperties.h"

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
      : fModuleLabel_v  ((size_t)(::larlite::data::kDATA_TYPE_MAX),std::vector<std::string>())
      , fDataReadFlag_v ((size_t)(::larlite::data::kDATA_TYPE_MAX),std::map<std::string,bool>())
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
      fModuleLabel_v[data_type].push_back(name);
      fDataReadFlag_v[data_type].insert(std::make_pair(name,false));
    }

    /// Accessor to the list of registered producer module labels
    std::vector<std::vector<std::string> > const& ModuleLabels() const { return fModuleLabel_v; }

    /// Function to be called @ end or beginning of each event
    void EventClear();

    /// Method to generate a data product ID for a specified type (through template) and producer module label's index
    template <class T>
    const ::larlite::product_id ProductID(size_t name_index) const;

    /// Core method: convert LArSoft data product (dh) to LArLite (lite_dh)
    template <class T>
    void ScanData(art::Handle<std::vector<T> > const &dh,
		  ::larlite::event_base* lite_dh);
    
    /// Core method: generate LArLite association data product and store (in lite_dh)
    template <class T, class U>
    void ScanAssociation(art::Event const& e,
			 art::Handle<std::vector<T> > &dh,
			 ::larlite::event_base* lite_dh) const;

  private:

    /// Internal utility function to locate art::Ptr map ... used to locate associated data product location
    template <class T>
    const std::map<art::Ptr<T>,std::pair<size_t,size_t> >& GetPtrMap() const;
    
    /// Internal utility function to look up associated object locater in LArLite data holder
    template <class T>
    bool LocateLiteProduct(art::Ptr<T> const ptr,
			   std::pair<size_t,size_t> &loc) const;

    /// Internal utility function to link LArSoft data product type to LArLite enum
    //template <class T>
    //const ::larlite::data::DataType_t LiteDataType() const;
    template <class T>
    const ::larlite::data::DataType_t LiteDataType() const;

    /// Internal utility function to look up the index number of given producer's name (for a specified type)
    size_t NameIndex(::larlite::data::DataType_t const data_type,
		     std::string const& name) const;

    /// Holder for producer module labels. Outer vector index = larlite::data::DataType_t
    std::vector<std::vector<std::string> > fModuleLabel_v;

    /// Boolean holder to tell us whether specific producer's data is read or not
    std::vector<std::map<std::string,bool> > fDataReadFlag_v;

    // art::Ptr local storage. Value = index of data product & index of label
    std::map<art::Ptr<::simb::MCTruth>,     std::pair<size_t,size_t> > fPtrIndex_mctruth;
    std::map<art::Ptr<::simb::GTruth>,      std::pair<size_t,size_t> > fPtrIndex_gtruth;
    std::map<art::Ptr<::simb::MCFlux>,      std::pair<size_t,size_t> > fPtrIndex_mcflux;
    std::map<art::Ptr<::simb::MCParticle>,  std::pair<size_t,size_t> > fPtrIndex_mcpart;
    std::map<art::Ptr<::sim::SimChannel>,   std::pair<size_t,size_t> > fPtrIndex_simch;
    std::map<art::Ptr<::sim::MCShower>,     std::pair<size_t,size_t> > fPtrIndex_mcshower;
    std::map<art::Ptr<::raw::RawDigit>,     std::pair<size_t,size_t> > fPtrIndex_rawdigit;
    std::map<art::Ptr<::recob::Wire>,       std::pair<size_t,size_t> > fPtrIndex_wire;
    std::map<art::Ptr<::recob::Hit>,        std::pair<size_t,size_t> > fPtrIndex_hit;
    std::map<art::Ptr<::recob::OpHit>,      std::pair<size_t,size_t> > fPtrIndex_ophit;
    std::map<art::Ptr<::recob::OpFlash>,    std::pair<size_t,size_t> > fPtrIndex_opflash;
    std::map<art::Ptr<::recob::Cluster>,    std::pair<size_t,size_t> > fPtrIndex_cluster;
    std::map<art::Ptr<::recob::Shower>,     std::pair<size_t,size_t> > fPtrIndex_shower;
    std::map<art::Ptr<::recob::Vertex>,     std::pair<size_t,size_t> > fPtrIndex_vertex;
    std::map<art::Ptr<::recob::Track>,      std::pair<size_t,size_t> > fPtrIndex_track;
    std::map<art::Ptr<::anab::CosmicTag>,   std::pair<size_t,size_t> > fPtrIndex_cosmictag;
    std::map<art::Ptr<::anab::Calorimetry>, std::pair<size_t,size_t> > fPtrIndex_calo;
    std::map<art::Ptr<::recob::SpacePoint>, std::pair<size_t,size_t> > fPtrIndex_sps;
    std::map<art::Ptr<::recob::EndPoint2D>, std::pair<size_t,size_t> > fPtrIndex_end2d;
    std::map<art::Ptr<::recob::Seed>,       std::pair<size_t,size_t> > fPtrIndex_seed;
    std::map<art::Ptr<::anab::ParticleID>,  std::pair<size_t,size_t> > fPtrIndex_partid;
    std::map<art::Ptr<::recob::PFParticle>, std::pair<size_t,size_t> > fPtrIndex_pfpart;

  };
}

#include "ScannerAlgo.template.h"

#endif
