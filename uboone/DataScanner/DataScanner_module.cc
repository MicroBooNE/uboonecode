/**
 * \file AnaProcess_module.cc
 *
 * \ingroup DataScanner
 * 
 * \brief Class definition file of DataScanner
 *
 * @author Kazu - Nevis 2013
 */

/** \addtogroup DataScanner

@{*/
#ifndef DataScanner_H
#define DataScanner_H

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "RawData/RawDigit.h"
#include "RecoBase/Wire.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Track.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Shower.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/EndPoint2D.h"
#include "AnalysisBase/Calorimetry.h"
#include "Simulation/SimChannel.h"
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"
#include "OpticalDetectorData/FIFOChannel.h"
#include "OpticalDetectorData/OpticalTypes.h"
#include "MCBase/MCShower.h"
//#include "RecoAlg/ClusterParamsAlg.h"
#include "Utilities/LArProperties.h"
#include "Utilities/GeometryUtilities.h"
#include "Utilities/DetectorProperties.h"

// ART includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"

// LArLight includes
#include <Base/Base-TypeDef.hh>
#include <DataFormat/DataFormat-TypeDef.hh>

// ROOT includes
#include "TTree.h"
#include "TPrincipal.h"

namespace datascanner {

  typedef enum lar_data_type {
    kLAR_MCTRUTH=0,
    kLAR_MCTRAJECTORY,
    kLAR_MCPARTICLE,
    kLAR_MCNEUTRINO,
    kLAR_SIMCHANNEL,
    kLAR_MCSHOWER,
    kLAR_WIRE,
    kLAR_FIFO,
    kLAR_PMTFIFO,
    kLAR_TPCFIFO,
    kLAR_HIT,
    kLAR_CLUSTER,
    kLAR_SHOWER,
    kLAR_SEED,
    kLAR_SPS,
    kLAR_TRACK,
    kLAR_CALO,
    kLAR_VERTEX,
    kLAR_END2D,
    kLAR_DATA_TYPE_MAX
  } lar_data_type_t ;

  /**
     \class DataScanner
     DataScanner module to copy LArSoft data contents into LArLight data formatted file
  */ 
  class DataScanner : public art::EDAnalyzer{

  public:

    /// Constructor
    DataScanner(const fhicl::ParameterSet&);

    /// Destructor
    virtual ~DataScanner(){}

    /// Function to be called before an event loop
    void beginJob();

    /// Function to be called per event
    void analyze (const art::Event&); 

  private:

    /// Function to read & store SimChannel data
    void ReadSimChannel(const art::Event& evt, const std::string mod_name, larlight::event_simch* data_ptr);

    /// Function to read & store MC truth shower object
    void ReadMCShower(const art::Event& evt, const std::string mod_name, larlight::event_mcshower* data_ptr);

    /// Function to read & store calibrated wire data
    void ReadWire(const art::Event& evt,  const std::string mod_name, larlight::event_wire* data_ptr);

    /// Function to read & store reconstructed hit data
    void ReadHit(const art::Event& evt,  const std::string mod_name, larlight::event_hit* data_ptr);

    /// Function to read & store reconstructed hit data
    void ReadCluster(const art::Event& evt,  const std::string mod_name, larlight::event_cluster* data_ptr, 
		     const std::vector<larlight::DATA::DATA_TYPE>& ass_types);

    /// Function to read & store PMT FIFO
    void ReadPMT(const art::Event& evt,  const std::string mod_name, larlight::event_pmtfifo* data_ptr);

    /// Function to read & store TPC FIFO
    void ReadTPC(const art::Event& evt,  const std::string mod_name, larlight::event_tpcfifo* data_ptr);

    /// Function to read & store spacepoints
    void ReadSPS(const art::Event& evt,  const std::string mod_name, larlight::event_sps* data_ptr,
		 const std::vector<larlight::DATA::DATA_TYPE>& ass_types);

    /// Function to read & store Tracking information
    void ReadTrack(const art::Event& evt,  const std::string mod_name, larlight::event_track* data_ptr,
		   const std::vector<larlight::DATA::DATA_TYPE>& ass_types);

    /// Function to read & store MCTruth information
    void ReadMCTruth(const art::Event& evt,  const std::string mod_name, larlight::event_mctruth* data_ptr);

    /// Function to read & store MCTruth information
    void ReadMCPartArray(const art::Event& evt,  const std::string mod_name, larlight::event_mcpart* data_ptr);

    /// Function to read & store Shower variables
    void ReadShower(const art::Event& evt,  const std::string mod_name, larlight::event_shower* data_ptr,
		    const std::vector<larlight::DATA::DATA_TYPE>& ass_types);

    /// Function to read & store Calorimetry variables
    void ReadCalorimetry(const art::Event& evt,  const std::string mod_name, larlight::event_calorimetry* data_ptr,
			 const std::vector<larlight::DATA::DATA_TYPE>& ass_types);

    /// Function to read & store Calorimetry variables
    void ReadVertex(const art::Event& evt,  const std::string mod_name, larlight::event_vertex* data_ptr,
		    const std::vector<larlight::DATA::DATA_TYPE>& ass_types);

    /// Function to read & store Calorimetry variables
    void ReadEndPoint2D(const art::Event& evt,  const std::string mod_name, larlight::event_endpoint2d* data_ptr,
			const std::vector<larlight::DATA::DATA_TYPE>& ass_types);

    /// Function to store user specific variables
    void StoreUserInfo(const art::Event& evt, larlight::event_user* data_ptr);

    /// Utility function for converting larlight::DATA::DATA_TYPE to lar_data_type_t
    datascanner::lar_data_type_t ConvertDataType(larlight::DATA::DATA_TYPE type);

    /// Utility function to parse module name string
    void ParseModuleName(std::vector<std::string> &mod_names, std::string name);

    /// Utility function to parse associated data types
    void ParseAssociationType(std::vector<std::vector<larlight::DATA::DATA_TYPE> >&ass_types, std::string name);

    /// Utility function to convert string to larlight::DATA::DATA_TYPE value
    larlight::DATA::DATA_TYPE Str2DataType(std::string name);

    /// Utility function to check if a particle is in the fiducial volume or not
    bool IsFV(Double_t x, Double_t y, Double_t z, Double_t t) const;

    /// One of main functions to read & copy association
    //void CopyAssociation(const art::Event& evt, const std::string mod_name, larlight::DATA::DATA_TYPE type, data_base* data_ptr)

    std::vector<TTree*>                    _trees;     ///< output data holder TTree
    std::vector<std::vector<std::string> > _mod_names; ///< input data production module names input from FCL file
    std::vector<larlight::event_base*>      _data_ptr;  ///< output data holder class object pointers
    
    std::vector<std::vector<std::vector<larlight::DATA::DATA_TYPE> > > _ass_types; ///< association data type to be specified via user

    std::vector<std::map<art::Ptr<recob::Hit>,        size_t > > _ass_map_hit;
    std::vector<std::map<art::Ptr<recob::Cluster>,    size_t > > _ass_map_cluster;
    std::vector<std::map<art::Ptr<recob::Shower>,     size_t > > _ass_map_shower;
    std::vector<std::map<art::Ptr<recob::Vertex>,     size_t > > _ass_map_vertex;
    std::vector<std::map<art::Ptr<recob::Track>,      size_t > > _ass_map_track;
    std::vector<std::map<art::Ptr<anab::Calorimetry>, size_t > > _ass_map_calo;
    std::vector<std::map<art::Ptr<recob::SpacePoint>, size_t > > _ass_map_sps;
    std::vector<std::map<art::Ptr<recob::EndPoint2D>, size_t > > _ass_map_end2d;

    typedef std::map<art::Ptr<recob::Hit>,        size_t >::const_iterator assmap_hit_citer;
    typedef std::map<art::Ptr<recob::Cluster>,    size_t >::const_iterator assmap_cluster_citer;
    typedef std::map<art::Ptr<recob::Shower>,     size_t >::const_iterator assmap_shower_citer;
    typedef std::map<art::Ptr<recob::Vertex>,     size_t >::const_iterator assmap_vertex_citer;
    typedef std::map<art::Ptr<recob::Track>,      size_t >::const_iterator assmap_track_citer;
    typedef std::map<art::Ptr<anab::Calorimetry>, size_t >::const_iterator assmap_calo_citer;
    typedef std::map<art::Ptr<recob::SpacePoint>, size_t >::const_iterator assmap_sps_citer;
    typedef std::map<art::Ptr<recob::EndPoint2D>, size_t >::const_iterator assmap_endp2d_citer;

    // ClusterParamsAlg module
    //cluster::ClusterParamsAlg fCParamsAlg;

    // Utility conversion factors
    double CONV_WIRE2CM; ///< wire -> length (cm) conversion factor
    double CONV_TIME2CM; ///< time -> length (cm) conversion factor

    Double_t _x_max; ///< Maximum X boundary of the detector
    Double_t _y_max; ///< Maximum X boundary of the detector
    Double_t _z_max; ///< Maximum X boundary of the detector
    Double_t _x_min; ///< Maximum X boundary of the detector
    Double_t _y_min; ///< Maximum X boundary of the detector
    Double_t _z_min; ///< Maximum X boundary of the detector
    Double_t _readout_startT; ///< Time at which readout window starts in G4 clock
    Double_t _readout_endT;   ///< Time at which readout window ends in G4 clock 
    Double_t _readout_freq;   ///< TPC sampling frequency in MHz
    Double_t _readout_size;   ///< Readout window size in readout time thick

    // Service modules
    art::ServiceHandle<util::DetectorProperties> _detp;
    art::ServiceHandle<util::LArProperties> _larp;
    art::ServiceHandle<geo::Geometry> _geo;

  };

} 

#endif//  DataScanner_H

// DataScanner.cc

namespace datascanner {
  DEFINE_ART_MODULE(DataScanner)
}

namespace datascanner {

  //#######################################################################################################
  DataScanner::DataScanner(fhicl::ParameterSet const& pset) : EDAnalyzer(pset),
							      _trees(larlight::DATA::DATA_TYPE_MAX,0), 
							      _data_ptr(larlight::DATA::DATA_TYPE_MAX,0)
  //#######################################################################################################
  {

    //fCParamsAlg.reconfigure(pset.get< fhicl::ParameterSet >("ClusterParamsAlg"));

    // Set detector boundaries
    art::ServiceHandle<geo::Geometry> geo;
    _y_max = geo->DetHalfHeight();
    _y_min = (-1.) * _y_max;
    _z_min = 0;
    _z_max = geo->DetLength();
    _x_min = 0;
    _x_max = 2.*(geo->DetHalfWidth());

    // These should be obtained from time service
    _readout_startT = -1.6e-3;
    _readout_endT   = 3.2e-3;
    _readout_freq   = 2.;
    _readout_size   = 3200.;

    // Initialize module name container
    for(size_t i=0; i<(size_t)(larlight::DATA::DATA_TYPE_MAX); i++) {
      
      _mod_names.push_back(std::vector<std::string>());

      _ass_map_hit.push_back     ( std::map<art::Ptr<recob::Hit>,size_t>()        );
      _ass_map_cluster.push_back ( std::map<art::Ptr<recob::Cluster>,size_t>()    );
      _ass_map_shower.push_back  ( std::map<art::Ptr<recob::Shower>,size_t>()     );
      _ass_map_vertex.push_back  ( std::map<art::Ptr<recob::Vertex>,size_t>()     );
      _ass_map_sps.push_back     ( std::map<art::Ptr<recob::SpacePoint>,size_t>() );
      _ass_map_track.push_back   ( std::map<art::Ptr<recob::Track>,size_t>()      );
      _ass_map_calo.push_back    ( std::map<art::Ptr<anab::Calorimetry>,size_t>() );
      _ass_map_end2d.push_back   ( std::map<art::Ptr<recob::EndPoint2D>,size_t>() );
    }

    // Obtain module names for input data
    // If a user set an empty string for these params, they are ignored for processing.
    ParseModuleName ( _mod_names[larlight::DATA::Bezier],               pset.get<std::string>("fModName_Bezier")               );
    ParseModuleName ( _mod_names[larlight::DATA::Kalman3DSPS],          pset.get<std::string>("fModName_Kalman3DSPS")          );
    ParseModuleName ( _mod_names[larlight::DATA::Kalman3DHit],          pset.get<std::string>("fModName_Kalman3DHit")          );
    ParseModuleName ( _mod_names[larlight::DATA::MCTruth],              pset.get<std::string>("fModName_MCTruth")              );
    ParseModuleName ( _mod_names[larlight::DATA::MCParticle],           pset.get<std::string>("fModName_MCParticle")           );
    ParseModuleName ( _mod_names[larlight::DATA::SpacePoint],           pset.get<std::string>("fModName_SpacePoint")           );
    ParseModuleName ( _mod_names[larlight::DATA::PMTFIFO],              pset.get<std::string>("fModName_FIFOChannel")          );
    ParseModuleName ( _mod_names[larlight::DATA::TPCFIFO],              pset.get<std::string>("fModName_RawDigit")             );
    ParseModuleName ( _mod_names[larlight::DATA::Wire],                 pset.get<std::string>("fModName_CalData")              );
    ParseModuleName ( _mod_names[larlight::DATA::CrawlerHit],           pset.get<std::string>("fModName_CrawlerHit")           );
    ParseModuleName ( _mod_names[larlight::DATA::GausHit],              pset.get<std::string>("fModName_GausHit")              );
    ParseModuleName ( _mod_names[larlight::DATA::APAHit],               pset.get<std::string>("fModName_APAHit")               );
    ParseModuleName ( _mod_names[larlight::DATA::FFTHit],               pset.get<std::string>("fModName_FFTHit")               );
    ParseModuleName ( _mod_names[larlight::DATA::RFFHit],               pset.get<std::string>("fModName_RFFHit")               );
    ParseModuleName ( _mod_names[larlight::DATA::CrawlerCluster],       pset.get<std::string>("fModName_CrawlerCluster")       );
    ParseModuleName ( _mod_names[larlight::DATA::RyanCluster],          pset.get<std::string>("fModName_RyanCluster")          );
    ParseModuleName ( _mod_names[larlight::DATA::DBCluster],            pset.get<std::string>("fModName_DBCluster")            );
    ParseModuleName ( _mod_names[larlight::DATA::FuzzyCluster],         pset.get<std::string>("fModName_FuzzyCluster")         );
    ParseModuleName ( _mod_names[larlight::DATA::HoughCluster],         pset.get<std::string>("fModName_HoughCluster")         );
    ParseModuleName ( _mod_names[larlight::DATA::ShowerAngleCluster],   pset.get<std::string>("fModName_ShowerAngleCluster")   );
    ParseModuleName ( _mod_names[larlight::DATA::Shower],               pset.get<std::string>("fModName_Shower")               );
    ParseModuleName ( _mod_names[larlight::DATA::FeatureVertex],        pset.get<std::string>("fModName_FeatureVertex")        );
    ParseModuleName ( _mod_names[larlight::DATA::HarrisVertex],         pset.get<std::string>("fModName_HarrisVertex")         );
    ParseModuleName ( _mod_names[larlight::DATA::FeatureEndPoint2D],    pset.get<std::string>("fModName_FeatureEndPoint2D")    );
    ParseModuleName ( _mod_names[larlight::DATA::HarrisEndPoint2D],     pset.get<std::string>("fModName_HarrisEndPoint2D")     );
    ParseModuleName ( _mod_names[larlight::DATA::Calorimetry],          pset.get<std::string>("fModName_Calorimetry")          );
    ParseModuleName ( _mod_names[larlight::DATA::MCShower],             pset.get<std::string>("fModName_MCShower")             );
    ParseModuleName ( _mod_names[larlight::DATA::SimChannel],           pset.get<std::string>("fModName_SimChannel")           );

    // Initialize association type input flags with correct length
    for(size_t i=0; i<(size_t)(larlight::DATA::DATA_TYPE_MAX); i++){

      std::vector<std::vector<larlight::DATA::DATA_TYPE> > ass_per_type;
      ass_per_type.reserve(_mod_names[i].size());

      for(size_t j=0; j<_mod_names[i].size(); ++j)

	ass_per_type.push_back(std::vector<larlight::DATA::DATA_TYPE>());

      _ass_types.push_back(ass_per_type);

    }

    ParseAssociationType ( _ass_types[larlight::DATA::Bezier],               pset.get<std::string>("fAssType_Bezier")               );
    ParseAssociationType ( _ass_types[larlight::DATA::Kalman3DSPS],          pset.get<std::string>("fAssType_Kalman3DSPS")          );
    ParseAssociationType ( _ass_types[larlight::DATA::Kalman3DHit],          pset.get<std::string>("fAssType_Kalman3DHit")          );
    ParseAssociationType ( _ass_types[larlight::DATA::MCTruth],              pset.get<std::string>("fAssType_MCTruth")              );
    ParseAssociationType ( _ass_types[larlight::DATA::MCParticle],           pset.get<std::string>("fAssType_MCParticle")           );
    ParseAssociationType ( _ass_types[larlight::DATA::SpacePoint],           pset.get<std::string>("fAssType_SpacePoint")           );
    ParseAssociationType ( _ass_types[larlight::DATA::PMTFIFO],              pset.get<std::string>("fAssType_FIFOChannel")          );
    ParseAssociationType ( _ass_types[larlight::DATA::TPCFIFO],              pset.get<std::string>("fAssType_RawDigit")             );
    ParseAssociationType ( _ass_types[larlight::DATA::Wire],                 pset.get<std::string>("fAssType_CalData")              );
    ParseAssociationType ( _ass_types[larlight::DATA::CrawlerHit],           pset.get<std::string>("fAssType_CrawlerHit")           );
    ParseAssociationType ( _ass_types[larlight::DATA::GausHit],              pset.get<std::string>("fAssType_GausHit")              );
    ParseAssociationType ( _ass_types[larlight::DATA::APAHit],               pset.get<std::string>("fAssType_APAHit")               );
    ParseAssociationType ( _ass_types[larlight::DATA::FFTHit],               pset.get<std::string>("fAssType_FFTHit")               );
    ParseAssociationType ( _ass_types[larlight::DATA::RFFHit],               pset.get<std::string>("fAssType_RFFHit")               );    
    ParseAssociationType ( _ass_types[larlight::DATA::CrawlerCluster],       pset.get<std::string>("fAssType_CrawlerCluster")       );
    ParseAssociationType ( _ass_types[larlight::DATA::RyanCluster],          pset.get<std::string>("fAssType_RyanCluster")          );
    ParseAssociationType ( _ass_types[larlight::DATA::DBCluster],            pset.get<std::string>("fAssType_DBCluster")            );
    ParseAssociationType ( _ass_types[larlight::DATA::FuzzyCluster],         pset.get<std::string>("fAssType_FuzzyCluster")         );
    ParseAssociationType ( _ass_types[larlight::DATA::HoughCluster],         pset.get<std::string>("fAssType_HoughCluster")         );
    ParseAssociationType ( _ass_types[larlight::DATA::ShowerAngleCluster],   pset.get<std::string>("fAssType_ShowerAngleCluster")   );
    ParseAssociationType ( _ass_types[larlight::DATA::Shower],               pset.get<std::string>("fAssType_Shower")               );
    ParseAssociationType ( _ass_types[larlight::DATA::FeatureVertex],        pset.get<std::string>("fAssType_FeatureVertex")        );
    ParseAssociationType ( _ass_types[larlight::DATA::HarrisVertex],         pset.get<std::string>("fAssType_HarrisVertex")         );
    ParseAssociationType ( _ass_types[larlight::DATA::FeatureEndPoint2D],    pset.get<std::string>("fAssType_FeatureEndPoint2D")    );
    ParseAssociationType ( _ass_types[larlight::DATA::HarrisEndPoint2D],     pset.get<std::string>("fAssType_HarrisEndPoint2D")     );
    ParseAssociationType ( _ass_types[larlight::DATA::Calorimetry],          pset.get<std::string>("fAssType_Calorimetry")          );
    ParseAssociationType ( _ass_types[larlight::DATA::MCShower],             pset.get<std::string>("fAssType_MCShower")             );
    ParseAssociationType ( _ass_types[larlight::DATA::SimChannel],           pset.get<std::string>("fAssType_SimChannel")           );

    // Next we make storage data class objects for those data types specified in fcl files.
    art::ServiceHandle<art::TFileService>  fileService;

    for(size_t i=0; i<(int)(larlight::DATA::DATA_TYPE_MAX); i++){

      larlight::DATA::DATA_TYPE type = (larlight::DATA::DATA_TYPE)i;

      // Check if a user provided an input module name for this data type.
      //if(_mod_names[i].size() || type==larlight::DATA::UserInfo){
      if(_mod_names[i].size()) {
	// Create TTree
	_trees[i] = fileService->make<TTree>(Form("%s_tree",larlight::DATA::DATA_TREE_NAME[i].c_str()),"");

	// Next, create data class objects
	switch(type){
	case larlight::DATA::MCTruth:
	  _data_ptr[i]=(larlight::event_base*)(new larlight::event_mctruth(type));
	  break;
	case larlight::DATA::MCParticle:
	  _data_ptr[i]=(larlight::event_base*)(new larlight::event_mcpart(type));
	  break;
	case larlight::DATA::Track:
	case larlight::DATA::Kalman3DSPS:
	case larlight::DATA::Kalman3DHit:
	case larlight::DATA::Bezier:
	  _data_ptr[i]=(larlight::event_base*)(new larlight::event_track(type));
	  break;
	case larlight::DATA::SpacePoint:
	  _data_ptr[i]=(larlight::event_base*)(new larlight::event_sps(type));
	  break;

	case larlight::DATA::TPCFIFO:
	  _data_ptr[i]=(larlight::event_base*)(new larlight::event_tpcfifo(type));
	  break;
	case larlight::DATA::PMTFIFO:
	  _data_ptr[i]=(larlight::event_base*)(new larlight::event_pmtfifo(type));
	  break;
	case larlight::DATA::Hit:
	case larlight::DATA::MCShowerHit:
	case larlight::DATA::CrawlerHit:
	case larlight::DATA::GausHit:
	case larlight::DATA::APAHit:
	case larlight::DATA::FFTHit:
	case larlight::DATA::RFFHit:
	  _data_ptr[i]=(larlight::event_base*)(new larlight::event_hit(type));
	  break;
	case larlight::DATA::Wire:
	  _data_ptr[i]=(larlight::event_base*)(new larlight::event_wire(type));
	  break;
	case larlight::DATA::Cluster:
	case larlight::DATA::MCShowerCluster:
	case larlight::DATA::RyanCluster:
	case larlight::DATA::CrawlerCluster:
	case larlight::DATA::DBCluster:
	case larlight::DATA::FuzzyCluster:
	case larlight::DATA::HoughCluster:
	case larlight::DATA::ShowerAngleCluster:
	  _data_ptr[i]=(larlight::event_base*)(new larlight::event_cluster(type));
	  break;
	case larlight::DATA::UserInfo:
	  _data_ptr[i]=(larlight::event_base*)(new larlight::event_user(type));
	  break;
	case larlight::DATA::RyanShower:
	case larlight::DATA::Shower:
	  _data_ptr[i]=(larlight::event_base*)(new larlight::event_shower(type));
	  break;
	case larlight::DATA::Calorimetry:
	  _data_ptr[i]=(larlight::event_base*)(new larlight::event_calorimetry(type));
	  break;
	case larlight::DATA::Vertex:
	case larlight::DATA::FeatureVertex:
	case larlight::DATA::HarrisVertex:
	  _data_ptr[i]=(larlight::event_base*)(new larlight::event_vertex(type));
	  break;
	case larlight::DATA::EndPoint2D:
	case larlight::DATA::FeatureEndPoint2D:
	case larlight::DATA::HarrisEndPoint2D:
	  _data_ptr[i]=(larlight::event_base*)(new larlight::event_endpoint2d(type));
	  break;
	case larlight::DATA::MCShower:
	  _data_ptr[i]=(larlight::event_base*)(new larlight::event_mcshower(type));
	  break;
	case larlight::DATA::SimChannel:
	  _data_ptr[i]=(larlight::event_base*)(new larlight::event_simch(type));
	  break;
	case larlight::DATA::MCTrajectory:
	case larlight::DATA::MCNeutrino:
	case larlight::DATA::FIFO:
	case larlight::DATA::Event:
	case larlight::DATA::Seed:
	case larlight::DATA::Pulse:
	case larlight::DATA::PMTPulse_FixedWin:
	case larlight::DATA::PMTPulse_ThresWin:
	case larlight::DATA::TPCPulse_FixedWin:
	case larlight::DATA::TPCPulse_ThresWin:
	case larlight::DATA::Trigger:
	case larlight::DATA::DATA_TYPE_MAX:
	  mf::LogError("DataScanner")<<Form("Data type %d not supported!",type);
	  break;
	}

	if(_data_ptr[i]) {
	  // Set TTree branch to the created data class object's address
	  _trees[i]->Branch(Form("%s_branch",larlight::DATA::DATA_TREE_NAME[i].c_str()),
			    _data_ptr[i]->GetName(),
			    &(_data_ptr[i]));
	}else
	  _trees[i]=0;
      }
    }
  }

  //**********************************************************************###
  datascanner::lar_data_type_t DataScanner::ConvertDataType(larlight::DATA::DATA_TYPE larlight_type)
  //**********************************************************************###
  {
    datascanner::lar_data_type_t type=kLAR_DATA_TYPE_MAX;
    
    switch(larlight_type){
    case larlight::DATA::MCTruth:
      type=kLAR_MCTRUTH; break;
    case larlight::DATA::MCTrajectory:
      type=kLAR_MCTRAJECTORY; break;
    case larlight::DATA::MCNeutrino:
      type=kLAR_MCNEUTRINO; break;
    case larlight::DATA::MCParticle:
      type=kLAR_MCPARTICLE; break;
    case larlight::DATA::Track:
    case larlight::DATA::Kalman3DSPS:
    case larlight::DATA::Kalman3DHit:
    case larlight::DATA::Bezier:
      type=kLAR_TRACK; break;
    case larlight::DATA::SpacePoint:
      type=kLAR_SPS; break;
      break;
    case larlight::DATA::FIFO:
      type=kLAR_FIFO; break;
    case larlight::DATA::TPCFIFO:
      type=kLAR_TPCFIFO; break;
    case larlight::DATA::PMTFIFO:
      type=kLAR_PMTFIFO; break;
    case larlight::DATA::Hit:
    case larlight::DATA::MCShowerHit:
    case larlight::DATA::CrawlerHit:
    case larlight::DATA::GausHit:
    case larlight::DATA::APAHit:
    case larlight::DATA::FFTHit:
    case larlight::DATA::RFFHit:
      type=kLAR_HIT; break;
    case larlight::DATA::Wire:
      type=kLAR_WIRE; break;
    case larlight::DATA::Cluster:
    case larlight::DATA::MCShowerCluster:
    case larlight::DATA::RyanCluster:
    case larlight::DATA::CrawlerCluster:
    case larlight::DATA::DBCluster:
    case larlight::DATA::FuzzyCluster:
    case larlight::DATA::HoughCluster:
    case larlight::DATA::ShowerAngleCluster:
      type=kLAR_CLUSTER; break;
    case larlight::DATA::RyanShower:
    case larlight::DATA::Shower:
      type=kLAR_SHOWER; break;
    case larlight::DATA::Calorimetry:
      type=kLAR_CALO; break;
    case larlight::DATA::Vertex:
    case larlight::DATA::FeatureVertex:
    case larlight::DATA::HarrisVertex:
      type=kLAR_VERTEX; break;
    case larlight::DATA::EndPoint2D:
    case larlight::DATA::FeatureEndPoint2D:
    case larlight::DATA::HarrisEndPoint2D:
      type=kLAR_END2D; break;
    case larlight::DATA::MCShower:
      type=kLAR_MCSHOWER; break;
    case larlight::DATA::SimChannel:
      type=kLAR_SIMCHANNEL; break;
    case larlight::DATA::UserInfo:
    case larlight::DATA::Event:
    case larlight::DATA::Seed:
    case larlight::DATA::Pulse:
    case larlight::DATA::PMTPulse_FixedWin:
    case larlight::DATA::PMTPulse_ThresWin:
    case larlight::DATA::TPCPulse_FixedWin:
    case larlight::DATA::TPCPulse_ThresWin:
    case larlight::DATA::Trigger:
    case larlight::DATA::DATA_TYPE_MAX:
    default:
      type=kLAR_DATA_TYPE_MAX;
    }
    
    return type;
    
  }
  
  //#######################################################################################################
  void DataScanner::ParseModuleName(std::vector<std::string> &mod_names, std::string name) 
  //#######################################################################################################
  {
    while(1){

      size_t pos = name.find(":");

      if(pos>=name.size()) break;

      mod_names.push_back(std::string(name.substr(0,pos)));

      name = name.substr(pos+1);
      
    }

    if(name.size()) mod_names.push_back(name);

  }

  //#######################################################################################################################
  void DataScanner::ParseAssociationType(std::vector<std::vector<larlight::DATA::DATA_TYPE> > &ass_types, std::string name)
  //#######################################################################################################################
  {
    std::vector<std::string> tmp_str_v;
    ParseModuleName(tmp_str_v,name);

    if(ass_types.size() < tmp_str_v.size()) {
      
      std::ostringstream msg;

      msg << std::endl
	  << Form("Problem with the association type specifier string: \"%s\".", name.c_str())
	  << std::endl
	  << Form("The sring contains more input modules than what is specified by producer module specifier string!")
	  << std::endl
	  << "This association won't be stored (ignored) ..." 
	  << std::endl
	  << std::endl;
	
      mf::LogError(__FUNCTION__) << msg.str();

      return;
    }
      

    for(size_t i=0; i<tmp_str_v.size(); ++i) {

      while(1){

	size_t pos = tmp_str_v.at(i).find(",");
	
	if(pos>=tmp_str_v.at(i).size()) break;

	std::string ass_name = tmp_str_v.at(i).substr(0,pos);

	tmp_str_v[i] = tmp_str_v.at(i).substr(pos+1);

	if(ass_name.empty()) continue;
	
	larlight::DATA::DATA_TYPE ass_type = Str2DataType(ass_name);

	if(ass_type!=larlight::DATA::DATA_TYPE_MAX) 

	  ass_types[i].push_back(ass_type);

	else

	  mf::LogError(__FUNCTION__) << Form("Invalid association type string provided: %s ... ignored!",ass_name.c_str());
	
      }

      if(tmp_str_v[i].size()) {

	larlight::DATA::DATA_TYPE ass_type = Str2DataType(tmp_str_v[i]);

	if(ass_type!=larlight::DATA::DATA_TYPE_MAX)

	  ass_types[i].push_back(ass_type);

	else

	  mf::LogError(__FUNCTION__) << Form("Invalid association type string provided: %s ... ignored!",tmp_str_v[i].c_str());

      }
      
    }
      
  }

  //#######################################################################################################################
  larlight::DATA::DATA_TYPE DataScanner::Str2DataType(std::string name)
  //#######################################################################################################################
  {

    // Dirty function to map-out string data type name to enum value

    larlight::DATA::DATA_TYPE type=larlight::DATA::DATA_TYPE_MAX;

    if( name=="Event" )
      type=larlight::DATA::Event; 
    else if( name=="mcshower" )
      type=larlight::DATA::MCShower;
    else if( name=="simch" )
      type=larlight::DATA::SimChannel;
    else if( name=="MCTruth" )
      type=larlight::DATA::MCTruth; 
    else if( name=="MCParticle" )
      type=larlight::DATA::MCParticle; 
    else if( name=="Wire" )
      type=larlight::DATA::Wire; 
    else if( name=="Hit" )
      type=larlight::DATA::Hit; 
    else if( name=="MCShowerHit" )
      type=larlight::DATA::MCShowerHit;
    else if( name=="CrawlerHit" )
      type=larlight::DATA::CrawlerHit; 
    else if( name=="GausHit" )
      type=larlight::DATA::GausHit; 
    else if( name=="APAHit" )
      type=larlight::DATA::APAHit; 
    else if( name=="FFTHit" )
      type=larlight::DATA::FFTHit; 
    else if( name=="RFFHit" )
      type=larlight::DATA::RFFHit; 
    else if( name=="Cluster" )
      type=larlight::DATA::Cluster; 
    else if( name=="MCShowerCluster" )
      type=larlight::DATA::MCShowerCluster;
    else if( name=="RyanCluster" )
      type=larlight::DATA::RyanCluster; 
    else if( name=="FuzzyCluster" )
      type=larlight::DATA::FuzzyCluster; 
    else if( name=="DBCluster" )
      type=larlight::DATA::DBCluster; 
    else if( name=="CrawlerCluster" )
      type=larlight::DATA::CrawlerCluster; 
    else if( name=="HoughCluster" )
      type=larlight::DATA::HoughCluster; 
    else if( name=="ShowerAngleCluster" )
      type=larlight::DATA::ShowerAngleCluster; 
    else if( name=="Seed" )
      type=larlight::DATA::Seed; 
    else if( name=="SpacePoint" )
      type=larlight::DATA::SpacePoint; 
    else if( name=="Track" )
      type=larlight::DATA::Track; 
    else if( name=="Bezier" )
      type=larlight::DATA::Bezier; 
    else if( name=="Kalman3DSPS" )
      type=larlight::DATA::Kalman3DSPS; 
    else if( name=="Kalman3DHit" )
      type=larlight::DATA::Kalman3DHit; 
    else if( name=="Shower" )
      type=larlight::DATA::Shower; 
    else if( name=="RyanShower" )
      type=larlight::DATA::RyanShower;
    else if( name=="Vertex" )
      type=larlight::DATA::Vertex;
    else if( name=="FeatureVertex" )
      type=larlight::DATA::FeatureVertex;
    else if( name=="HarrisVertex" )
      type=larlight::DATA::HarrisVertex;
    else if( name=="EndPoint2D" )
      type=larlight::DATA::EndPoint2D;
    else if( name=="FeatureEndPoint2D")
      type=larlight::DATA::FeatureEndPoint2D;
    else if( name=="HarrisEndPoint2D")
      type=larlight::DATA::HarrisEndPoint2D;
    else if( name=="Calorimetry" )
      type=larlight::DATA::Calorimetry; 
    else if( name=="FIFO" )
      type=larlight::DATA::FIFO; 
    else if( name=="PMTFIFO" )
      type=larlight::DATA::PMTFIFO; 
    else if( name=="TPCFIFO" )
      type=larlight::DATA::TPCFIFO; 
    else if( name=="Trigger" )
      type=larlight::DATA::Trigger; 
    else {
      mf::LogError(__FUNCTION__) << Form("Unsupported type name: %s",name.c_str());
      type=larlight::DATA::DATA_TYPE_MAX; 
    }
    
    return type;

  }

  //****************************************************************************
  bool DataScanner::IsFV(Double_t x, Double_t y, Double_t z,Double_t t)  const 
  //****************************************************************************
  {
    // x, y, z should be in mm, time in ns

    // Hard volume cut
    if( x > _x_max || x < _x_min || z > _z_max || z < _z_min || y > _y_max || y < _y_min ) return false;

    // Now compare x-coordinate to make sure this point is inside the readout window.
    Double_t drift_v = (_x_max - _x_min) / (_readout_size / _readout_freq);  // unit = mm / us 
    Double_t drift_t = (x - _x_min) / drift_v; // unit = us

    // Charge arrives @ wire-plane @ drift_t + t ... compute this in second scale
    Double_t arrival_t = drift_t*1.e-6 + t*1.e-9;

    return ( _readout_startT < arrival_t && arrival_t < _readout_endT);
  }

  //#######################################################################################################
  void DataScanner::beginJob()
  //#######################################################################################################
  {
    // Some detector constants
    CONV_WIRE2CM = _geo->WirePitch(0,1,0);    //wire pitch in cm
    CONV_TIME2CM = (_detp->SamplingRate()/1000.) * _larp->DriftVelocity(_larp->Efield(),_larp->Temperature());
  }

  //#######################################################################################################
  void DataScanner::analyze(const art::Event& evt) 
  //#######################################################################################################
  {

    // Loop over data type to initialize storage
    for(size_t i=0; i<(size_t)(larlight::DATA::DATA_TYPE_MAX); i++){

      if(!(_trees[i])) continue;

      // Reset data
      _data_ptr[i]->clear_data();

      // Clear association pointer address map
      _ass_map_hit[i].clear();
      _ass_map_cluster[i].clear();
      _ass_map_shower[i].clear();
      _ass_map_sps[i].clear();
      _ass_map_track[i].clear();
      _ass_map_calo[i].clear();

    }
    

    // Loop over data type to store contents
    for(size_t i=0; i<(size_t)(larlight::DATA::DATA_TYPE_MAX); i++){
      
      // If data pointer is not set, we don't have to fill this data type
      if(!(_trees[i])) continue;
      
      // Fill common variables such as run, subrun, and event id
      _data_ptr[i]->set_run      ( evt.id().run()    );
      _data_ptr[i]->set_subrun   ( evt.id().subRun() );
      _data_ptr[i]->set_event_id ( evt.id().event()  );

      // Handle different kind of data (class wise)
      larlight::DATA::DATA_TYPE type = (larlight::DATA::DATA_TYPE)i;
      switch(type){

      case larlight::DATA::Track:
      case larlight::DATA::Kalman3DSPS:
      case larlight::DATA::Kalman3DHit:
      case larlight::DATA::Bezier:
	// Data types to be stored in event_track class
	for(size_t j=0; j<_mod_names[i].size(); ++j)
	  ReadTrack(evt, _mod_names[i][j], (larlight::event_track*)(_data_ptr[i]),_ass_types[i][j]);;
	break;

      case larlight::DATA::MCTruth:
	// Data type to be stored in event_mctruth class
	for(size_t j=0; j<_mod_names[i].size(); ++j)
	  ReadMCTruth     (evt, _mod_names[i][j], (larlight::event_mctruth*)(_data_ptr[i]));
	break;
      case larlight::DATA::MCParticle:
	// Data type to be stored in event_mctruth class
	for(size_t j=0; j<_mod_names[i].size(); ++j)
	  ReadMCPartArray (evt, _mod_names[i][j], (larlight::event_mcpart*)(_data_ptr[i]));
	break;
      case larlight::DATA::SpacePoint:
 	// Data type to be stored in event_sps class
	for(size_t j=0; j<_mod_names[i].size(); ++j)
	  ReadSPS(evt,_mod_names[i][j], (larlight::event_sps*)(_data_ptr[i]),_ass_types[i][j]);
	break;

      case larlight::DATA::FIFO:
      case larlight::DATA::PMTFIFO:
	// Data type to be stored in event_pmtfifo class
	for(size_t j=0; j<_mod_names[i].size(); ++j)
	  ReadPMT(evt,_mod_names[i][j], (larlight::event_pmtfifo*)(_data_ptr[i]));
	break;

      case larlight::DATA::TPCFIFO:
	// Data type to be stored in event_tpcfifo class
	for(size_t j=0; j<_mod_names[i].size(); ++j)
	  ReadTPC(evt,_mod_names[i][j], (larlight::event_tpcfifo*)(_data_ptr[i]));
	break;	

      case larlight::DATA::Wire:
	// Data type to be stored in event_wire class
	for(size_t j=0; j<_mod_names[i].size(); ++j)
	  ReadWire(evt,_mod_names[i][j],(larlight::event_wire*)(_data_ptr[i]));
	break;

      case larlight::DATA::Hit:
      case larlight::DATA::MCShowerHit:
      case larlight::DATA::CrawlerHit:
      case larlight::DATA::GausHit:
      case larlight::DATA::APAHit:
      case larlight::DATA::FFTHit:
      case larlight::DATA::RFFHit:
	// Data type to be stored in event_hit class
	for(size_t j=0; j<_mod_names[i].size(); ++j) 
	  ReadHit(evt,_mod_names[i][j],(larlight::event_hit*)(_data_ptr[i]));
	break;

      case larlight::DATA::Cluster:
      case larlight::DATA::MCShowerCluster:
      case larlight::DATA::RyanCluster:
      case larlight::DATA::DBCluster:
      case larlight::DATA::FuzzyCluster:
      case larlight::DATA::HoughCluster:
      case larlight::DATA::CrawlerCluster:
      case larlight::DATA::ShowerAngleCluster:
	// Data type to be stored in event_cluster class
	for(size_t j=0; j<_mod_names[i].size(); ++j) 
	  ReadCluster(evt, _mod_names[i][j], (larlight::event_cluster*)(_data_ptr[i]), _ass_types[i][j]);
	break;
	
      case larlight::DATA::UserInfo:
	// Data type to be stored in user_info class
	//StoreUserInfo(evt,(larlight::event_user*)(_data_ptr[i]));
	break;

      case larlight::DATA::RyanShower:
      case larlight::DATA::Shower:
	for(size_t j=0; j<_mod_names[i].size(); ++j)
	  ReadShower(evt,_mod_names[i][j],(larlight::event_shower*)(_data_ptr[i]),_ass_types[i][j]);
	break;

      case larlight::DATA::Calorimetry:
	for(size_t j=0; j<_mod_names[i].size(); ++j)
	  ReadCalorimetry(evt,_mod_names[i][j],(larlight::event_calorimetry*)(_data_ptr[i]),_ass_types[i][j]);
	break;

      case larlight::DATA::Vertex:
      case larlight::DATA::FeatureVertex:
      case larlight::DATA::HarrisVertex:
	for(size_t j=0; j<_mod_names[i].size(); ++j)
	  ReadVertex(evt,_mod_names[i][j],(larlight::event_vertex*)(_data_ptr[i]),_ass_types[i][j]);
	break;
	
      case larlight::DATA::EndPoint2D:
      case larlight::DATA::FeatureEndPoint2D:
      case larlight::DATA::HarrisEndPoint2D:
	for(size_t j=0; j<_mod_names[i].size(); ++j)
	  ReadEndPoint2D(evt,_mod_names[i][j],(larlight::event_endpoint2d*)(_data_ptr[i]),_ass_types[i][j]);
	break;
      case larlight::DATA::MCShower:
	for(size_t j=0; j<_mod_names[i].size(); ++j)
	  ReadMCShower(evt,_mod_names[i][j],(larlight::event_mcshower*)(_data_ptr[i]));
	break;
      case larlight::DATA::SimChannel:
	for(size_t j=0; j<_mod_names[i].size(); ++j)
	  ReadSimChannel(evt,_mod_names[i][j],(larlight::event_simch*)(_data_ptr[i]));
	break;
      case larlight::DATA::MCTrajectory:
      case larlight::DATA::MCNeutrino:
      case larlight::DATA::Event:
      case larlight::DATA::Seed:
      case larlight::DATA::Pulse:
      case larlight::DATA::PMTPulse_FixedWin:
      case larlight::DATA::PMTPulse_ThresWin:
      case larlight::DATA::TPCPulse_FixedWin:
      case larlight::DATA::TPCPulse_ThresWin:
      case larlight::DATA::Trigger:
      case larlight::DATA::DATA_TYPE_MAX:
	break;
      }

      // Fill this TTree
      _trees[i]->Fill(); 
    }
  }

  //##################################################################################################################
  void DataScanner::ReadSimChannel(const art::Event& evt, const std::string mod_name, larlight::event_simch* data_ptr)
  //##################################################################################################################
  {
    art::Handle<std::vector<sim::SimChannel> > schArray;
    evt.getByLabel(mod_name,schArray);

    if(!schArray.isValid()) return;

    for(size_t i=0; i<schArray->size(); ++i ) {

      const art::Ptr<sim::SimChannel> sch_ptr(schArray,i);

      std::map<unsigned short,std::vector<sim::IDE> > sch_map(sch_ptr->TDCIDEMap());

      larlight::simch light_sch;
      light_sch.set_channel(sch_ptr->Channel());

      for(auto sch_iter = sch_map.begin(); sch_iter!=sch_map.end(); ++sch_iter) {
	
	unsigned short tdc = (*sch_iter).first;

	for(auto const this_ide : (*sch_iter).second) {

	  larlight::ide light_ide;
	  light_ide.trackID = this_ide.trackID;
	  light_ide.numElectrons = this_ide.numElectrons;
	  light_ide.energy = this_ide.energy;
	  light_ide.x = this_ide.x;
	  light_ide.y = this_ide.y;
	  light_ide.z = this_ide.z;

	  light_sch.add_ide(tdc,light_ide);
	}
      }
      
      data_ptr->push_back(light_sch);
    }
  }
  
  //###################################################################################################################
  void DataScanner::ReadMCShower(const art::Event& evt, const std::string mod_name, larlight::event_mcshower* data_ptr)
  //###################################################################################################################
  {

    art::Handle<std::vector<sim::MCShower> > mcsHandle;
    evt.getByLabel(mod_name,mcsHandle);
    if(!mcsHandle.isValid()) return; 

    for(auto const& mcs : *mcsHandle) {

      larlight::mcshower light_prof;
      light_prof.SetMotherID(mcs.MotherPDGID(), mcs.MotherTrackID());
      
      //light_prof.SetMotherAngles(mcs.phiMother, mcs.thetaMother);
      //mcs.uAngleMother, mcs.vAngleMother, mcs.wAngleMother);

      light_prof.SetMotherPoint(mcs.MotherPosition());

      light_prof.SetMotherProcess(mcs.MotherCreationProcess());

      light_prof.SetMotherMomentum(mcs.MotherMomentum());

      light_prof.SetDaughterTrackList(mcs.DaughterTrackID());

      //light_prof.SetDaughterAngles(mcs.phiDaughter, mcs.thetaDaughter);
      //mcs.uAngleDaughter, mcs.vAngleDaughter, mcs.wAngleDaughter);

      light_prof.SetDaughterMomentum(mcs.DaughterMomentum());
      light_prof.SetDaughterPosition(mcs.DaughterPosition());

      light_prof.SetPlaneCharge(mcs.Charge());

      //light_prof.SetEdepVtx(mcs.vtxEdep);

      data_ptr->push_back(light_prof);
    }
    
  }
  

  //#######################################################################################################
  void DataScanner::ReadWire(const art::Event& evt, 
			     const std::string mod_name, 
			     larlight::event_wire* data_ptr){
  //#######################################################################################################

    art::Handle< std::vector<recob::Wire> > wireArray;
    evt.getByLabel(mod_name,wireArray);

    if(!wireArray.isValid()) return;

    for(size_t i=0; i<wireArray->size(); i++){

      const art::Ptr<recob::Wire> wire_ptr(wireArray,i);
      
      larlight::wire wire_light(wire_ptr->Signal(),
				wire_ptr->Channel(),
				(larlight::GEO::View_t)(wire_ptr->View()),
				(larlight::GEO::SigType_t)(wire_ptr->SignalType()));
      
      data_ptr->push_back(wire_light);
      
    }

  }

  //#######################################################################################################
  void DataScanner::ReadHit(const art::Event& evt, 
			    const std::string mod_name, 
			    larlight::event_hit* data_ptr){
  //#######################################################################################################

    art::Handle< std::vector<recob::Hit> > hitArray;
    evt.getByLabel(mod_name,hitArray);

    if(!hitArray.isValid()) return;

    for(size_t i=0; i<hitArray->size(); i++){
      
      art::Ptr<recob::Hit> hit_ptr(hitArray,i);

      larlight::hit hit_light(data_ptr->data_type());

      hit_light.set_waveform(hit_ptr->fHitSignal);
      hit_light.set_times(hit_ptr->StartTime(),
			  hit_ptr->PeakTime(),
			  hit_ptr->EndTime());
      hit_light.set_times_err(hit_ptr->SigmaStartTime(),
			      hit_ptr->SigmaPeakTime(),
			      hit_ptr->SigmaEndTime());
      hit_light.set_charge(hit_ptr->Charge(),hit_ptr->Charge(true));
      hit_light.set_charge_err(hit_ptr->SigmaCharge(),hit_ptr->SigmaCharge(true));
      hit_light.set_multiplicity(hit_ptr->Multiplicity());
      hit_light.set_channel(hit_ptr->Channel());
      hit_light.set_wire(hit_ptr->WireID().Wire);
      hit_light.set_fit_goodness(hit_ptr->GoodnessOfFit());
      hit_light.set_view((larlight::GEO::View_t)(hit_ptr->View()));
      hit_light.set_sigtype((larlight::GEO::SigType_t)(hit_ptr->SignalType()));

      // Store address map for downstream association
      _ass_map_hit[data_ptr->data_type()][hit_ptr] = data_ptr->size();

      data_ptr->push_back(hit_light);
      
    }

  }

  //#######################################################################################################
  void DataScanner::ReadCluster(const art::Event& evt, 
				const std::string mod_name, 
				larlight::event_cluster* data_ptr,
				const std::vector<larlight::DATA::DATA_TYPE>& ass_types){
  //#######################################################################################################

    art::Handle<std::vector<recob::Cluster> > clusterArray;
    evt.getByLabel(mod_name, clusterArray);
    if(!clusterArray.isValid()) return;

    art::FindManyP<recob::Hit> hit_m(clusterArray, evt, mod_name);

    for(size_t i=0; i<clusterArray->size(); ++i) {

      const art::Ptr<recob::Cluster> cluster_ptr(clusterArray,i);
      
      larlight::cluster cluster_light(data_ptr->data_type());
      cluster_light.set_charge(cluster_ptr->Charge());
      cluster_light.set_dtdw(cluster_ptr->dTdW());
      cluster_light.set_dqdw(cluster_ptr->dQdW());
      cluster_light.set_dtdw_err(cluster_ptr->SigmadTdW());
      cluster_light.set_dqdw_err(cluster_ptr->SigmadQdW());
      cluster_light.set_id(cluster_ptr->ID());
      cluster_light.set_view((larlight::GEO::View_t)(cluster_ptr->View()));
      cluster_light.set_start_vtx(cluster_ptr->StartPos());
      cluster_light.set_end_vtx(cluster_ptr->EndPos());
      cluster_light.set_start_vtx_err(cluster_ptr->SigmaStartPos());
      cluster_light.set_end_vtx_err(cluster_ptr->SigmaEndPos());

      // Store associated hits
      const std::vector<art::Ptr<recob::Hit> > hit_v = hit_m.at(i);
      
      for(auto const type : ass_types) {
	
	lar_data_type_t lar_type = ConvertDataType(type);
	
	if(lar_type != kLAR_HIT) {
	  
	  mf::LogError(__FUNCTION__) 
	    << Form("Storage of association from \"%s\" to \"%s\" not implemented!",
		    larlight::DATA::DATA_TREE_NAME[data_ptr->data_type()].c_str(),
		    larlight::DATA::DATA_TREE_NAME[type].c_str());
	  
	  continue;
	}
	
	// Find them
	std::vector<unsigned int> ass_index;
	
	ass_index.reserve(hit_v.size());
	
	for(auto const hit_ptr : hit_v){
	  
	  assmap_hit_citer iter = _ass_map_hit[type].find(hit_ptr);
	  
	  if(iter!=_ass_map_hit[type].end())
	    
	    ass_index.push_back((*iter).second);
	  
	}
	
	// Store them
	cluster_light.add_association(type,ass_index);
      }
      
      // Store address map for downstream association
      _ass_map_cluster[data_ptr->data_type()][cluster_ptr]=data_ptr->size();

      // Process the case of association among itself
      

      data_ptr->push_back(cluster_light);
    }

  }

  //#######################################################################################################
  void DataScanner::ReadPMT(const art::Event& evt, 
			    const std::string mod_name, 
			    larlight::event_pmtfifo* data_ptr){
  //#######################################################################################################

    art::Handle<std::vector<optdata::FIFOChannel> > pmtArray;
    evt.getByLabel(mod_name, pmtArray);
    if(!pmtArray.isValid()) return;
    
    for(size_t i=0; i<pmtArray->size(); ++i) {

      const art::Ptr<optdata::FIFOChannel> fifo_ptr(pmtArray,i);

      optdata::Optical_Category_t cat(fifo_ptr->Category());
      larlight::FEM::DISCRIMINATOR disc_id = larlight::FEM::DISC_MAX;
      switch(cat){
      case optdata::kFEMCosmicHighGain:
      case optdata::kFEMCosmicLowGain:
      case optdata::kCosmicPMTTrigger:
	disc_id = larlight::FEM::COSMIC_DISC;
	break;
      case optdata::kFEMBeamHighGain:
      case optdata::kFEMBeamLowGain:
      case optdata::kBeamPMTTrigger:
	disc_id = larlight::FEM::BEAM_DISC;
	break;
      case optdata::kUndefined:
      case optdata::kHighGain:
      case optdata::kLowGain:
      default:
	disc_id = larlight::FEM::DISC_MAX;
      }

      larlight::pmtfifo fifo_light(fifo_ptr->ChannelNumber(),
				0,
				0,
				fifo_ptr->Frame(),
				fifo_ptr->TimeSlice(),
				disc_id,
				data_ptr->data_type(),
				*fifo_ptr);
      
      data_ptr->push_back(fifo_light);
    }
    
  }

  //#######################################################################################################
  void DataScanner::ReadTPC(const art::Event& evt,
			    const std::string mod_name,
			    larlight::event_tpcfifo* data_ptr){
  //#######################################################################################################

    art::Handle<std::vector<raw::RawDigit> > tpcArray;
    evt.getByLabel(mod_name, tpcArray);
    if(!tpcArray.isValid()) return;


    larlight::tpcfifo light_wf;

    for(size_t i=0; i<tpcArray->size(); ++i) {

      const art::Ptr<raw::RawDigit> wf(tpcArray,i);

      light_wf.clear_data();

      // Copy ADC samples
      for(unsigned int i=0; i<wf->NADC(); ++i)

	light_wf.push_back(wf->ADC(i));
      
      // Copy channel-wise parameters
      light_wf.set_channel_number(wf->Channel());
      geo::SigType_t sigtype = _geo->SignalType(wf->Channel());
      if (sigtype == geo::kInduction)
	light_wf.set_signal(larlight::GEO::kInduction);
      else if (sigtype == geo::kCollection)
	light_wf.set_signal(larlight::GEO::kCollection);
      else
	light_wf.set_signal(larlight::GEO::kMysteryType);
      geo::View_t viewtype = _geo->View(wf->Channel());
      if (viewtype == geo::kU)
	light_wf.set_plane(larlight::GEO::kU);
      else if (viewtype == geo::kV)
	light_wf.set_plane(larlight::GEO::kV);
      else if (viewtype == geo::kZ)
	light_wf.set_plane(larlight::GEO::kZ);
      else
	light_wf.set_plane(larlight::GEO::kUnknown);

      data_ptr->push_back(light_wf);
  
    }

  }


  //#######################################################################################################
  void DataScanner::ReadSPS(const art::Event& evt, 
			    const std::string mod_name, 
			    larlight::event_sps* data_ptr,
			    const std::vector<larlight::DATA::DATA_TYPE>& ass_types){
  //#######################################################################################################

    art::Handle<std::vector<recob::SpacePoint> > spsArray;
    evt.getByLabel(mod_name, spsArray);
    if(!spsArray.isValid()) return;

    for(size_t i=0; i<spsArray->size(); ++i){

      const art::Ptr<recob::SpacePoint> sps_ptr(spsArray,i);
      
      larlight::spacepoint sps_light(sps_ptr->ID(),
				     sps_ptr->XYZ()[0],    sps_ptr->XYZ()[1],    sps_ptr->XYZ()[2],
				     sps_ptr->ErrXYZ()[0], sps_ptr->ErrXYZ()[1], sps_ptr->ErrXYZ()[2],
				     sps_ptr->Chisq());

      // Store associated hits/clusters
      for(auto const type : ass_types) {

	// Find them
	lar_data_type_t lar_type = ConvertDataType(type);
	if(lar_type != kLAR_HIT && lar_type!=kLAR_CLUSTER) {

	  mf::LogError(__FUNCTION__) 
	    << Form("Storage of association from \"%s\" to \"%s\" not implemented!",
		    larlight::DATA::DATA_TREE_NAME[data_ptr->data_type()].c_str(),
		    larlight::DATA::DATA_TREE_NAME[type].c_str());

	  continue;
	}

	std::vector<unsigned int> ass_index;

	if(lar_type == kLAR_CLUSTER) {

	  const std::vector<art::Ptr<recob::Cluster> > cluster_v = art::FindManyP<recob::Cluster>(spsArray,evt,mod_name).at(i);

	  ass_index.reserve(cluster_v.size());

	  for(auto const cluster_ptr : cluster_v){
	    
	    assmap_cluster_citer iter = _ass_map_cluster[type].find(cluster_ptr);
	    
	    if(iter!=_ass_map_cluster[type].end())
	      
	      ass_index.push_back((*iter).second);
	  }

	}else if(lar_type == kLAR_HIT) {

	  const std::vector<art::Ptr<recob::Hit> > hit_v = art::FindManyP<recob::Hit>(spsArray,evt,mod_name).at(i);

	  ass_index.reserve(hit_v.size());

	  for(auto const hit_ptr : hit_v){
	    
	    assmap_hit_citer iter = _ass_map_hit[type].find(hit_ptr);
	    
	    if(iter!=_ass_map_hit[type].end())
	      
	      ass_index.push_back((*iter).second);
	  }

	}

	// Store them
	sps_light.add_association(type,ass_index);	
      }

      // Store address map for downstream association      
      _ass_map_sps[data_ptr->data_type()][sps_ptr] = data_ptr->size();

      data_ptr->push_back(sps_light);
    }

  }

  //#######################################################################################################
  void DataScanner::ReadMCPartArray(const art::Event& evt, 
				    const std::string mod_name, 
				    larlight::event_mcpart* data_ptr){
  //#######################################################################################################

    art::Handle<std::vector<simb::MCParticle> > mcpArray;
    evt.getByLabel(mod_name, mcpArray);
    if(!mcpArray.isValid()) return;

    for(size_t i=0; i < mcpArray->size(); ++i) {

      const art::Ptr<simb::MCParticle> mcp_ptr(mcpArray,i);

      larlight::mcpart light_part(mcp_ptr->TrackId(), mcp_ptr->PdgCode(), mcp_ptr->Process(), 
				  mcp_ptr->Mother(), mcp_ptr->Mass(), mcp_ptr->StatusCode());

      light_part.SetPolarization(mcp_ptr->Polarization());
      light_part.SetRescatter(mcp_ptr->Rescatter());
      light_part.SetWeight(mcp_ptr->Weight());
      
      for(size_t k=0; k<(size_t)(mcp_ptr->NumberDaughters()); ++k)
	
	light_part.AddDaughter(mcp_ptr->Daughter(k));
      
      // Add trajectory points
      larlight::mctrajectory light_track(larlight::DATA::MCTrajectory);
      light_track.reserve(mcp_ptr->NumberTrajectoryPoints());
      
      bool   loopInFV=false;
      size_t start_FV=0;
      for(size_t l=0; l<(size_t)mcp_ptr->NumberTrajectoryPoints(); ++l) {
	
	light_track.push_back(mcp_ptr->Position(l),mcp_ptr->Momentum(l));
	
	// Record fiducial volume tracking
	bool inFV = IsFV(mcp_ptr->Vx(l),mcp_ptr->Vy(l),mcp_ptr->Vz(l),mcp_ptr->T(l));
	
	if(!loopInFV) {
	  
	  if(inFV) { loopInFV=true; start_FV=l; }
	  
	}else if(!inFV) {
	  light_part.AddFiducialTrack(start_FV,l-1);
	  loopInFV=false;
	}
      }
      if(loopInFV){ light_part.AddFiducialTrack(start_FV,mcp_ptr->NumberTrajectoryPoints()-1); }
      
      light_part.SetTrajectory(light_track);
      
      // Save
      data_ptr->push_back(light_part);
    }

  }

  //#######################################################################################################
  void DataScanner::ReadMCTruth(const art::Event& evt, 
				const std::string mod_name, 
				larlight::event_mctruth* data_ptr){
  //#######################################################################################################

    art::Handle<std::vector<simb::MCTruth> > mctArray;
    evt.getByLabel(mod_name, mctArray);
    if(!mctArray.isValid()) return;

    for(size_t i=0; i < mctArray->size(); ++i) {

      const art::Ptr<simb::MCTruth> mct_ptr(mctArray,i);
      
      larlight::mctruth light_mct(data_ptr->data_type());

      // Generator (Origin) ID
      light_mct.SetOrigin((larlight::MC::Origin_t) mct_ptr->Origin() );

      // Particle Information
      for(size_t j=0; j<(size_t)(mct_ptr->NParticles()); ++j) {

	const simb::MCParticle lar_part(mct_ptr->GetParticle(j));

	larlight::mcpart light_part(lar_part.TrackId(), lar_part.PdgCode(), lar_part.Process(), 
				    lar_part.Mother(), lar_part.Mass(), lar_part.StatusCode());

	light_part.SetPolarization(lar_part.Polarization());
	light_part.SetRescatter(lar_part.Rescatter());
	light_part.SetWeight(lar_part.Weight());

	for(size_t k=0; k<(size_t)(lar_part.NumberDaughters()); ++k)

	  light_part.AddDaughter(lar_part.Daughter(k));

	// Add trajectory points
	larlight::mctrajectory light_track(larlight::DATA::MCTrajectory);
	light_track.reserve(lar_part.NumberTrajectoryPoints());
	
	bool   loopInFV=false;
	size_t start_FV=0;
	for(size_t l=0; l<(size_t)lar_part.NumberTrajectoryPoints(); ++l) {
	  
	  light_track.push_back(lar_part.Position(l),lar_part.Momentum(l));

	  // Record fiducial volume tracking
	  bool inFV = IsFV(lar_part.Vx(l),lar_part.Vy(l),lar_part.Vz(l),lar_part.T(l));
 
	  if(!loopInFV) {

	    if(inFV) { loopInFV=true; start_FV=l; }

	  }else if(!inFV) {
	      light_part.AddFiducialTrack(start_FV,l-1);
	      loopInFV=false;
	  }
	}
	if(loopInFV){ light_part.AddFiducialTrack(start_FV,lar_part.NumberTrajectoryPoints()-1); }

	light_part.SetTrajectory(light_track);
	
	light_mct.Add(light_part);
      }

      // Neutrino Information
      const simb::MCNeutrino lar_nu(mct_ptr->GetNeutrino());
      
      light_mct.SetNeutrino( lar_nu.CCNC(),
			     lar_nu.Mode(),
			     lar_nu.InteractionType(),
			     lar_nu.Target(),
			     lar_nu.HitNuc(),
			     lar_nu.HitQuark(),
			     lar_nu.W(),
			     lar_nu.X(),
			     lar_nu.Y(),
			     lar_nu.QSqr() );
      
      // Save
      data_ptr->push_back(light_mct);
    }
  }

  //#######################################################################################################
  void DataScanner::ReadShower(const art::Event& evt, 
			       const std::string mod_name, 
			       larlight::event_shower* data_ptr,
			       const std::vector<larlight::DATA::DATA_TYPE>& ass_types){
  //#######################################################################################################
    
    art::Handle< std::vector<recob::Shower> > showerArray;
    evt.getByLabel(mod_name, showerArray);
    if(!showerArray.isValid()) return;

    art::FindManyP<recob::Cluster> cluster_m(showerArray, evt, mod_name);

    for(size_t i=0; i < showerArray->size(); ++i) {

      const art::Ptr<recob::Shower> shower_ptr(showerArray,i);

      larlight::shower light_shower(data_ptr->data_type());

      light_shower.set_id(shower_ptr->ID());
      //light_shower.set_total_charge(shower_ptr->TotalCharge());
      light_shower.set_direction(shower_ptr->Direction());
      light_shower.set_direction_err(shower_ptr->DirectionErr());
      //light_shower.set_max_width(shower_ptr->MaxWidthX(),shower_ptr->MaxWidthY());
      //light_shower.set_distance_max_width(shower_ptr->DistanceMaxWidth());

      // Store associated clusters
      const std::vector<art::Ptr<recob::Cluster> > cluster_v = cluster_m.at(i);

      for(auto const type : ass_types) {

	lar_data_type_t lar_type = ConvertDataType(type);
	if(lar_type!=kLAR_CLUSTER) {

	  mf::LogError(__FUNCTION__) 
	    << Form("Storage of association from \"%s\" to \"%s\" not implemented!",
		    larlight::DATA::DATA_TREE_NAME[data_ptr->data_type()].c_str(),
		    larlight::DATA::DATA_TREE_NAME[type].c_str());
	  continue;
	}

	// Find them
	std::vector<unsigned int> ass_index;

	ass_index.reserve(cluster_v.size());

	for(auto const cluster_ptr : cluster_v){
	  
	  assmap_cluster_citer iter = _ass_map_cluster[type].find(cluster_ptr);
	  
	  if(iter!=_ass_map_cluster[type].end())
	    
	    ass_index.push_back((*iter).second);
	  
	}
	
	// Store them
	light_shower.add_association(type,ass_index);	
      }

      // Store address map for downstream association
      _ass_map_shower[data_ptr->data_type()][shower_ptr]=data_ptr->size();
      
      data_ptr->push_back(light_shower);
    }

  }

  //##########################################################################################################################
  void DataScanner::ReadCalorimetry(const art::Event& evt, 
				    const std::string mod_name, 
				    larlight::event_calorimetry* data_ptr,
				    const std::vector<larlight::DATA::DATA_TYPE>& ass_types){
  //##########################################################################################################################

    art::Handle< std::vector<anab::Calorimetry> > caloArray;
    evt.getByLabel(mod_name, caloArray);
    if(!caloArray.isValid()) return;

    for(size_t i=0; i < caloArray->size(); ++i) {

      const art::Ptr<anab::Calorimetry> calo_ptr(caloArray,i);
      
      larlight::calorimetry light_calo(data_ptr->data_type());

      light_calo.set_dedx(calo_ptr->dEdx());
      light_calo.set_dqdx(calo_ptr->dQdx());
      light_calo.set_residual_range(calo_ptr->ResidualRange());
      light_calo.set_deadwire_range(calo_ptr->DeadWireResRC());
      light_calo.set_kinetic_energy(calo_ptr->KineticEnergy());
      light_calo.set_range(calo_ptr->Range());
      light_calo.set_track_pitch(calo_ptr->TrkPitchVec());

      // Store associated shower/track

      for(auto const type : ass_types) {

	// Find them
        lar_data_type_t lar_type = ConvertDataType(type);
        if(lar_type != kLAR_TRACK && lar_type!=kLAR_SHOWER) {

	  mf::LogError(__FUNCTION__)
            << Form("Storage of association from \"%s\" to \"%s\" not implemented!",
                    larlight::DATA::DATA_TREE_NAME[data_ptr->data_type()].c_str(),
                    larlight::DATA::DATA_TREE_NAME[type].c_str());

	  continue;
        }

	std::vector<unsigned int> ass_index;

	if(lar_type == kLAR_TRACK) {

	  const std::vector<art::Ptr<recob::Track> > track_v = art::FindManyP<recob::Track>(caloArray,evt,mod_name).at(i);

	  ass_index.reserve(track_v.size());

          for(auto const track_ptr : track_v){

            assmap_track_citer iter = _ass_map_track[type].find(track_ptr);

            if(iter!=_ass_map_track[type].end())

	      ass_index.push_back((*iter).second);
          }

	}else if(lar_type == kLAR_SHOWER) {

	  const std::vector<art::Ptr<recob::Shower> > shower_v = art::FindManyP<recob::Shower>(caloArray,evt,mod_name).at(i);

	  ass_index.reserve(shower_v.size());
	  
          for(auto const shower_ptr : shower_v){

            assmap_shower_citer iter = _ass_map_shower[type].find(shower_ptr);

            if(iter!=_ass_map_shower[type].end())

	      ass_index.push_back((*iter).second);
          }
	}
	// Store them
        light_calo.add_association(type,ass_index);
      }

      // Store address map for downstream association
      _ass_map_calo[data_ptr->data_type()][calo_ptr]=data_ptr->size();

      data_ptr->push_back(light_calo);

    }

  }

  //##########################################################################################################################
  void DataScanner::ReadVertex(const art::Event& evt, 
			       const std::string mod_name, 
			       larlight::event_vertex* data_ptr,
			       const std::vector<larlight::DATA::DATA_TYPE>& ass_types){
  //##########################################################################################################################

    art::Handle< std::vector<recob::Vertex> > vtxArray;
    evt.getByLabel(mod_name, vtxArray);
    if(!vtxArray.isValid()) return;


    Double_t xyz[3]={0};
    
    for(size_t i=0; i < vtxArray->size(); ++i) {

      const art::Ptr<recob::Vertex> vtx_ptr(vtxArray,i);

      vtx_ptr->XYZ(xyz);

      // get vertex info
      larlight::vertex light_vtx(xyz,
				 vtx_ptr->ID(),
				 data_ptr->data_type());
      // store associations
      for(auto const type : ass_types) {

	// Find them
        lar_data_type_t lar_type = ConvertDataType(type);
        if(lar_type != kLAR_HIT && lar_type!=kLAR_TRACK && lar_type!=kLAR_SHOWER) {

	  mf::LogError(__FUNCTION__)
            << Form("Storage of association from \"%s\" to \"%s\" not implemented!",
                    larlight::DATA::DATA_TREE_NAME[data_ptr->data_type()].c_str(),
                    larlight::DATA::DATA_TREE_NAME[type].c_str());

	  continue;
        }

	std::vector<unsigned int> ass_index;

	if(lar_type == kLAR_HIT) {

	  // Store associated hit, cluster, and spacepoint
	  const std::vector<art::Ptr<recob::Hit> > hit_v = art::FindManyP<recob::Hit>(vtxArray, evt, mod_name).at(i);

	  ass_index.reserve(hit_v.size());

          for(auto const hit_ptr : hit_v){

            assmap_hit_citer iter = _ass_map_hit[type].find(hit_ptr);

            if(iter!=_ass_map_hit[type].end())

	      ass_index.push_back((*iter).second);
          }

	}else if(lar_type == kLAR_TRACK) {

	  const std::vector<art::Ptr<recob::Track> > track_v = art::FindManyP<recob::Track>(vtxArray, evt, mod_name).at(i);

	  ass_index.reserve(track_v.size());
	  
          for(auto const track_ptr : track_v){

            assmap_track_citer iter = _ass_map_track[type].find(track_ptr);

            if(iter!=_ass_map_track[type].end())

	      ass_index.push_back((*iter).second);
          }

	}else if(lar_type == kLAR_SHOWER) {

	  const std::vector<art::Ptr<recob::Shower> > shower_v = art::FindManyP<recob::Shower>(vtxArray,evt,mod_name).at(i);

	  ass_index.reserve(shower_v.size());

          for(auto const shower_ptr : shower_v){

            assmap_shower_citer iter = _ass_map_shower[type].find(shower_ptr);

            if(iter!=_ass_map_shower[type].end())

	      ass_index.push_back((*iter).second);
          }
	}
	// Store them
	light_vtx.add_association(type,ass_index);
      }

      // Store address map for downstream association
      _ass_map_vertex[data_ptr->data_type()][vtx_ptr]=data_ptr->size();

      // Store data
      data_ptr->push_back(light_vtx);
      
    }
    
  }

  //##########################################################################################################################
  void DataScanner::ReadEndPoint2D(const art::Event& evt, 
				   const std::string mod_name, 
				   larlight::event_endpoint2d* data_ptr,
				   const std::vector<larlight::DATA::DATA_TYPE>& ass_types){
  //##########################################################################################################################

    art::Handle< std::vector<recob::EndPoint2D> > end2dArray;
    evt.getByLabel(mod_name, end2dArray);
    if(!end2dArray.isValid()) return;

    for(size_t i=0; i < end2dArray->size(); ++i) {

      const art::Ptr<recob::EndPoint2D> end2d_ptr(end2dArray,i);

      // get vertex info
      larlight::endpoint2d light_end2d(end2d_ptr->DriftTime(),
				       end2d_ptr->WireID().Wire,
				       end2d_ptr->Strength(),
				       end2d_ptr->ID(),
				       (larlight::GEO::View_t)(end2d_ptr->View()),
				       end2d_ptr->Charge(),
				       data_ptr->data_type());

      // store associations
      for(auto const type : ass_types) {

	// Find them
        lar_data_type_t lar_type = ConvertDataType(type);
        if(lar_type != kLAR_HIT) {

	  mf::LogError(__FUNCTION__)
            << Form("Storage of association from \"%s\" to \"%s\" not implemented!",
                    larlight::DATA::DATA_TREE_NAME[data_ptr->data_type()].c_str(),
                    larlight::DATA::DATA_TREE_NAME[type].c_str());

	  continue;
        }

	std::vector<unsigned int> ass_index;

	if(lar_type == kLAR_HIT) {

	  // Store associated hit
	  const std::vector<art::Ptr<recob::Hit> > hit_v = art::FindManyP<recob::Hit>(end2dArray, evt, mod_name).at(i);

	  ass_index.reserve(hit_v.size());

          for(auto const hit_ptr : hit_v){

            assmap_hit_citer iter = _ass_map_hit[type].find(hit_ptr);

            if(iter!=_ass_map_hit[type].end())

	      ass_index.push_back((*iter).second);
          }

	}
	// Store them
	light_end2d.add_association(type,ass_index);
      }

      // Store address map for downstream association
      _ass_map_end2d[data_ptr->data_type()][end2d_ptr]=data_ptr->size();

      // Store data
      data_ptr->push_back(light_end2d);
      
    }
    
  }



  //#######################################################################################################
  void DataScanner::ReadTrack(const art::Event& evt, 
			      const std::string mod_name, 
			      larlight::event_track* data_ptr,
			      const std::vector<larlight::DATA::DATA_TYPE>& ass_types){
  //#######################################################################################################

    art::Handle< std::vector<recob::Track> > trackArray;
    evt.getByLabel(mod_name, trackArray);
    if(!trackArray.isValid()) return;

    for(size_t i=0; i < trackArray->size(); ++i){

      // Obtain recob::Track object pointer
      const art::Ptr<recob::Track> track_ptr(trackArray,i);

      // Prepare storage track object
      larlight::track track_light(data_ptr->data_type());
      
      //
      // Start copying data
      //
      
      // ID
      track_light.set_track_id ( track_ptr->ID()   );

      // Direction & points
      for(size_t i=0; i<track_ptr->NumberTrajectoryPoints(); i++) {

	track_light.add_vertex     (track_ptr->LocationAtPoint(i));

	track_light.add_direction  (track_ptr->DirectionAtPoint(i));

      }

      // Covariance
      for(size_t i=0; i<track_ptr->NumberCovariance(); i++)

	track_light.add_covariance (track_ptr->CovarianceAtPoint(i));

      // Momentum
      for(size_t i=0; i<track_ptr->NumberFitMomentum(); i++)

	track_light.add_momentum   (track_ptr->MomentumAtPoint(i));
      
      for(auto const type : ass_types) {

	// Find them
        lar_data_type_t lar_type = ConvertDataType(type);
        if(lar_type != kLAR_HIT && lar_type!=kLAR_CLUSTER && lar_type!=kLAR_SPS) {

	  mf::LogError(__FUNCTION__)
            << Form("Storage of association from \"%s\" to \"%s\" not implemented!",
                    larlight::DATA::DATA_TREE_NAME[data_ptr->data_type()].c_str(),
                    larlight::DATA::DATA_TREE_NAME[type].c_str());

	  continue;
        }

	std::vector<unsigned int> ass_index;

	if(lar_type == kLAR_HIT) {

	  // Store associated hit, cluster, and spacepoint
	  const std::vector<art::Ptr<recob::Hit> > hit_v = art::FindManyP<recob::Hit>(trackArray, evt, mod_name).at(i);

	  ass_index.reserve(hit_v.size());

          for(auto const hit_ptr : hit_v){

            assmap_hit_citer iter = _ass_map_hit[type].find(hit_ptr);

            if(iter!=_ass_map_hit[type].end())

	      ass_index.push_back((*iter).second);
          }

	}else if(lar_type == kLAR_CLUSTER) {

	  const std::vector<art::Ptr<recob::Cluster> > cluster_v = art::FindManyP<recob::Cluster>(trackArray, evt, mod_name).at(i);

	  ass_index.reserve(cluster_v.size());
	  
          for(auto const cluster_ptr : cluster_v){

            assmap_cluster_citer iter = _ass_map_cluster[type].find(cluster_ptr);

            if(iter!=_ass_map_cluster[type].end())

	      ass_index.push_back((*iter).second);
          }

	}else if(lar_type == kLAR_SPS) {

	  const std::vector<art::Ptr<recob::SpacePoint> > sps_v = art::FindManyP<recob::SpacePoint>(trackArray,evt,mod_name).at(i);

	  ass_index.reserve(sps_v.size());

          for(auto const sps_ptr : sps_v){

            assmap_sps_citer iter = _ass_map_sps[type].find(sps_ptr);

            if(iter!=_ass_map_sps[type].end())

	      ass_index.push_back((*iter).second);
          }
	}
	// Store them
        track_light.add_association(type,ass_index);
      }

      // Store address map for downstream association
      _ass_map_track[data_ptr->data_type()][track_ptr]=data_ptr->size();      

      // Store this track
      data_ptr->push_back(track_light);

    }

    // Done
  }

  //#########################################################################
  void DataScanner::StoreUserInfo(const art::Event& evt, larlight::event_user* data_ptr)
  //###########################################################################
  {
    // This function does any necessary operation to store extra variables in user_info.

    larlight::user_info my_ui(data_ptr->data_type());

    // Block to fill ShowerAngleCluster related variables
    /*
    if(_data_ptr[larlight::DATA::ShowerAngleCluster]) {

      for(auto const mod_name : _mod_names[larlight::DATA::ShowerAngleCluster]) {

	std::vector<const recob::Cluster*> clusterArray;
	try{
	  evt.getView(mod_name,clusterArray);
	}catch (art::Exception const& e) {
	  if (e.categoryCode() != art::errors::ProductNotFound ) throw;
	  continue;
	}

	// Loop over individual clusters. In this loop, we fill:
	// (1) variables calculated by ClusterParamsAlg::isShower() function	
	art::Handle< std::vector<recob::Cluster> > clusterListHandle;
	evt.getByLabel(mod_name,clusterListHandle);
	art::FindManyP<recob::Hit> fmh(clusterListHandle, evt, mod_name);
	for(size_t i=0; i < clusterListHandle->size(); i++){

	  art::Ptr<recob::Cluster> iCluster (clusterListHandle,i);
	  std::vector< art::Ptr<recob::Hit> > iHitList = fmh.at(i);

	  // Prepare output bector
	  double lineslope  = iCluster->dQdW();
	  double start_wire = iCluster->StartPos()[0];
	  double start_time = iCluster->StartPos()[1];
	  double end_wire   = iCluster->EndPos()[0];
	  double end_time   = iCluster->EndPos()[1];

	  // Do an equivalent calculation as ClusterParamsAlg::isShower()
	  
	  double length = TMath::Sqrt( pow((start_wire - end_wire),2) * CONV_WIRE2CM + 
				       pow((start_time - end_time),2) * CONV_TIME2CM );

	  double xangle = fCParamsAlg.Get2DAngleForHit(start_wire, start_time, iHitList);

	  if(xangle >  90) xangle -= 180;
	  if(xangle < -90) xangle += 180;

	  double HighBin,LowBin,invHighBin,invLowBin,altWeight;
	  double PrincipalEigenvalue=1.,ModPV=-900;
	  double multihit=0; 
	  double ModHitDensity=0;

	  // OffAxisHits calculation
	  fCParamsAlg.FindDirectionWeights(lineslope, 
					   start_wire, start_time,
					   end_wire, end_time,
					   iHitList,
					   HighBin, LowBin, invHighBin, invLowBin, &altWeight);
	  altWeight /= length;
	  
	  // Principal value calculation
	  TPrincipal mypc(2,"D");
	  fCParamsAlg.GetPrincipal(iHitList,&mypc);
	  PrincipalEigenvalue = (*mypc.GetEigenValues())[0];
	  ModPV = TMath::Log(1-PrincipalEigenvalue);

	  // MultiHit calculation
	  multihit = fCParamsAlg.MultiHitWires(iHitList) / length;

	  // Modified hit density calculation
	  double HitDensity = iHitList.size()/length;
	  ModHitDensity = HitDensity - (1.82 * cosh(3.14/180*abs(xangle)-1.24) - 1.56);

	  // Store
	  my_ui.append("vModHitDensity",ModHitDensity);
	  my_ui.append("vMultiHit",multihit);
	  my_ui.append("vPrincipalHD",ModPV);
	  my_ui.append("vOffAxisHit",altWeight);
	} // end of cluster loop

	// Loop over a set of clusters. In this loop, we store:
	// (1) matched cluster combination indexes
	// (2) spacepoint associated with matching
	art::Handle< std::vector<art::PtrVector<recob::Cluster> > >    clusterAssociationHandle;	
	art::Handle< std::vector<art::PtrVector<recob::SpacePoint> > > spsAssociationHandle;
	evt.getByLabel(mod_name,clusterAssociationHandle);
	evt.getByLabel(mod_name,spsAssociationHandle);

	bool sps_exist = spsAssociationHandle.isValid();

	// Check both vector have the same length
	if(sps_exist &&
	   clusterAssociationHandle->size() && 
	   (spsAssociationHandle->size() != clusterAssociationHandle->size())) {
	  
	  mf::LogError(__PRETTY_FUNCTION__) << "Stored cluster & SPS vector size for shower different!";
	  
	  continue;
	}

	// Store number of showers
	my_ui.store("numShower",(int)(clusterAssociationHandle->size()));
	for(size_t i=0; i < clusterAssociationHandle->size(); ++i){

	  const art::PtrVector<recob::Cluster>    cluster_set = clusterAssociationHandle->at(i);

	  // Store cluster ID
	  for(size_t j=0; j<cluster_set.size(); ++j) {

	    const art::Ptr<recob::Cluster> cluster = cluster_set.at(j);

	    my_ui.append("vShowerClusterID",(int)(cluster->ID()));
	    
	  }

	  // Store # of SPS in THIS shower, then loop over SPS and append to a grand vector
	  if(sps_exist){
	    const art::PtrVector<recob::SpacePoint> sps_set     = spsAssociationHandle->at(i);
	    my_ui.append("vNumShowerSPS",(int)(sps_set.size()));
	    for(size_t j=0; j<sps_set.size(); ++j) {
	      
	      const art::Ptr<recob::SpacePoint> sps = sps_set.at(j);
	      
	      my_ui.append("vShowerSPS_X",(double)(sps->XYZ()[0]));
	      my_ui.append("vShowerSPS_Y",(double)(sps->XYZ()[1]));
	      my_ui.append("vShowerSPS_Z",(double)(sps->XYZ()[2]));
	      
	      my_ui.append("vShowerSPS_Xe",(double)(sps->ErrXYZ()[0]));
	      my_ui.append("vShowerSPS_Ye",(double)(sps->ErrXYZ()[1]));
	      my_ui.append("vShowerSPS_Ze",(double)(sps->ErrXYZ()[2]));

	      my_ui.append("vShowerSPS_Chi2",(double)(sps->Chisq()));
	    }
	  }
	}

      } // end of input module loop
    }
    */
    // store user_info
    data_ptr->push_back(my_ui);
  }

  
} // namespace 

/** @} */ // end of doxygen group
