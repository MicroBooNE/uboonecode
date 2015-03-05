////////////////////////////////////////////////////////////////////////
// $Id: DBSCANfinderAna.cxx,v 1.36 2010/09/15  bpage Exp $
//
// module to create a TTree for analysis
//
// \author tjyang@fnal.gov, sowjanyag@phys.ksu.edu
//
////////////////////////////////////////////////////////////////////////
// To reduce memory usage:
// [x] create the data structure connected to the tree only when needed
// [x] reduce the size of the elemental items (Double_t => Float_t could damage precision)
// [x] create a different structure for each tracker, allocate only what needed
// [x] use variable size array buffers for each tracker datum instead of [kMaxTrack]
// [x] turn the truth/GEANT information into vectors
// [ ] move hit_trkid into the track information, remove kMaxTrackers
// [ ] turn the hit information into vectors (~1 MB worth), remove kMaxHits
// [ ] fill the tree branch by branch
// 
// Current implementation:
// There is one tree only, with one set of branches for each tracking algorithm.
// The data structure which hosts the addresses of the tree branches is
// dynamically allocated on demand, and it can be optionally destroyed at the
// end of each event.
// The data structure (AnalysisTreeDataStruct) firectly contains the truth and
// simulation information as C arrays. The data from tracking algorithms is the
// largest, and it is contained in a C++ vector of structures (TrackDataStruct),
// one per algorithm. These structures can also be allocated on demand.
// Each of these structures is connected to a set of branches, one branch per
// data member. Data members are vectors of numbers or vectors of fixed-size
// C arrays. The vector index represents the tracks reconstructed by the
// algorithm, and each has a fixed size pool for hits (do ROOT t3rees support
// branches with more than one dimension with variable size?).
// The data structures can assign default values to their data, connect to a
// ROOT tree (creating the branches they need) and resize.
// The AnalysisTreeDataStruct is constructed with as many tracking algorithms as
// there are named in the module configuration (even if they are not backed by
// any available tracking data).
// By default construction, TrackDataStruct is initialized in a state which does
// not allow any track (maximum tracks number is zero), and in such state trying
// to connect to a tree has no effect. This is done so that the
// AnalysisTreeDataStruct can be initialized first (and with unusable track data
// structures), and then the TrackDataStruct instances are initialized one by
// one when the number of tracks needed is known.
// A similar mechanism is implemented for the truth information.
// 
// The "UseBuffers: false" mode assumes that on each event a new
// AnalysisTreeDataStruct is created with unusable tracker data, connected to
// the ROOT tree (the addresses of the available branches are assigned), then
// each of the tracking algorithm data is resized to host the correct number
// of reconstructed tracks and connected to the tree. Then the normal process of
// filling the event data and then the tree take place. Finally, the whole
// data structure is freed and the tree is left in a invalid state (branch
// addresses are invalid). It could be possible to make the tree in a valid
// state by resetting the addresses, but there is no advantage in that.
// 
// The "UseBuffers: true" mode assumes that on the first event a new
// AnalysisTreeDataStruct is created and used just as in the other mode
// described above. At the end of the first event, the data structure is left
// around (and the tree is in a valid state). On the next event, all the
// addresses are checked, then for each tracker the data is resized to
// accomodate the right number of tracks for tis event. If the memory is
// increased, the address will be changed. All the branches are reconnected to
// the data structure, and the procedure goes on as normal.
// 
// Note that reducing the maximum number of tracks in a TrackDataStruct does not
// necessarily make memory available, because of how std::vector::resize()
// works; that feature can be implemented, but it currently has not been.
// 
// The BoxedArray<> class is a wrapper around a normal C array; it is needed
// to be able to include such structure in a std::vector. This container
// requires its objects to be default-constructable and copy-constructable,
// and a C array is neither. BoxedArray<> is: the default construction leaves it
// uninitialized (for speed reasons) while the copy construction is performed
// as in a Plain Old Data structure (memcpy; really!).
// 
////////////////////////////////////////////////////////////////////////

#ifndef ANALYSISTREE_H
#define ANALYSISTREE_H

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/FindMany.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geometry/Geometry.h"
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCFlux.h"
#include "Simulation/SimChannel.h"
#include "Simulation/AuxDetSimChannel.h"
#include "AnalysisBase/Calorimetry.h"
#include "AnalysisBase/ParticleID.h"
#include "RawData/RawDigit.h"
#include "RawData/BeamInfo.h"
#include "Utilities/LArProperties.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"
#include "SummaryData/POTSummary.h"
#include "MCCheater/BackTracker.h"
#include "RecoBase/Track.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "RecoBase/EndPoint2D.h"
#include "RecoBase/Vertex.h"
#include "SimpleTypesAndConstants/geo_types.h"
#include "RecoObjects/BezierTrack.h"
#include "RecoAlg/TrackMomentumCalculator.h"
#include "AnalysisBase/CosmicTag.h"
#include "AnalysisBase/FlashMatch.h"
	

#include <cstring> // std::memcpy()
#include <vector>
#include <map>
#include <iterator> // std::begin(), std::end()
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional> // std::mem_fun_ref
#include <typeinfo>

#include "TTree.h"
#include "TTimeStamp.h"

constexpr int kNplanes       = 3;     //number of wire planes
constexpr int kMaxHits       = 25000; //maximum number of hits;
constexpr int kMaxTrackHits  = 2000;  //maximum number of hits on a track
constexpr int kMaxTrackers   = 15;    //number of trackers passed into fTrackModuleLabel
constexpr unsigned short kMaxVertices   = 100;    //max number of 3D vertices
constexpr unsigned short kMaxAuxDets = 4; ///< max number of auxiliary detector cells per MC particle

/// total_extent\<T\>::value has the total number of elements of an array
template <typename T>
struct total_extent {
  using value_type = size_t;
  static constexpr value_type value
    = sizeof(T) / sizeof(typename std::remove_all_extents<T>::type);
}; // total_extent<>


namespace microboone {

  /// Data structure with all the tree information.
  /// 
  /// Can connect to a tree, clear its fields and resize its data.
  class AnalysisTreeDataStruct {
      public:
    
    /// A wrapper to a C array (needed to embed an array into a vector)
    template <typename Array_t>
    class BoxedArray {
        protected:
      Array_t array; // actual data
      
        public:
      using This_t = BoxedArray<Array_t>;
      typedef typename std::remove_all_extents<Array_t>::type Data_t;
      
      BoxedArray() {} // no initialization
      BoxedArray(const This_t& from)
        { std::memcpy((char*) &(data()), (char*) &(from.data()), sizeof(Array_t)); }
      
      Array_t& data() { return array; }
      const Array_t& data() const { return array; }
      
      //@{
      /// begin/end interface
      static constexpr size_t size() { return total_extent<Array_t>::value; }
      Data_t* begin() { return reinterpret_cast<Data_t*>(&array); }
      const Data_t* begin() const { return reinterpret_cast<const Data_t*>(&array); }
      Data_t* end() { return begin() + size(); }
      const Data_t* end() const { return begin() + size(); }
      //@}
      
      //@{
      /// Array interface
      auto operator[] (size_t index) -> decltype(*array) { return array[index]; }
      auto operator[] (size_t index) const -> decltype(*array) { return array[index]; }
      auto operator+ (ptrdiff_t index) -> decltype(&*array) { return array + index; }
      auto operator+ (ptrdiff_t index) const -> decltype(&*array) { return array + index; }
      auto operator- (ptrdiff_t index) -> decltype(&*array) { return array - index; }
      auto operator- (ptrdiff_t index) const -> decltype(&*array) { return array - index; }
      auto operator* () -> decltype(*array) { return *array; }
      auto operator* () const -> decltype(*array) { return *array; }
      
      operator decltype(&array[0]) () { return &array[0]; }
      operator decltype(&array[0]) () const { return &array[0]; }
      //@}
      
    }; // BoxedArray
    
    /// Tracker algorithm result
    /// 
    /// Can connect to a tree, clear its fields and resize its data.
    class TrackDataStruct {
        public:
      /* Data structure size:
       *
       * TrackData_t<Short_t>                    :  2  bytes/track
       * TrackData_t<Float_t>                    :  4  bytes/track
       * PlaneData_t<Float_t>, PlaneData_t<Int_t>: 12  bytes/track
       * HitData_t<Float_t>                      : 24k bytes/track
       * HitCoordData_t<Float_t>                 : 72k bytes/track
       */
      template <typename T>
      using TrackData_t = std::vector<T>;
      template <typename T>
      using PlaneData_t = std::vector<BoxedArray<T[kNplanes]>>;
      template <typename T>
      using HitData_t = std::vector<BoxedArray<T[kNplanes][kMaxTrackHits]>>;
      template <typename T>
      using HitCoordData_t = std::vector<BoxedArray<T[kNplanes][kMaxTrackHits][3]>>;
      
      size_t MaxTracks; ///< maximum number of storable tracks
      
      Short_t  ntracks;             //number of reconstructed tracks
      PlaneData_t<Float_t>    trkke;
      PlaneData_t<Float_t>    trkrange;
      PlaneData_t<Int_t>      trkidtruth;  //true geant trackid
      PlaneData_t<Short_t>    trkorigin;   //_ev_origin 0: unknown, 1: cosmic, 2: neutrino, 3: supernova, 4: singles
      PlaneData_t<Int_t>      trkpdgtruth; //true pdg code
      PlaneData_t<Float_t>    trkefftruth; //completeness
      PlaneData_t<Float_t>    trksimIDEenergytruth;
      PlaneData_t<Float_t>    trksimIDExtruth;
      PlaneData_t<Float_t>    trksimIDEytruth;
      PlaneData_t<Float_t>    trksimIDEztruth;
      PlaneData_t<Float_t>    trkpurtruth; //purity of track
      PlaneData_t<Float_t>    trkpitchc;
      PlaneData_t<Short_t>    ntrkhits;
      HitData_t<Float_t>      trkdedx;
      HitData_t<Float_t>      trkdqdx;
      HitData_t<Float_t>      trkresrg;
      HitCoordData_t<Float_t> trkxyz;

      // more track info
      TrackData_t<Short_t> trkId;
      TrackData_t<Short_t> trkncosmictags_tagger;
      TrackData_t<Float_t> trkcosmicscore_tagger;
      TrackData_t<Short_t> trkcosmictype_tagger;
      TrackData_t<Short_t> trkncosmictags_flashmatch;
      TrackData_t<Float_t> trkcosmicscore_flashmatch;
      TrackData_t<Short_t> trkcosmictype_flashmatch;
      TrackData_t<Float_t> trkstartx;     // starting x position.
      TrackData_t<Float_t> trkstarty;     // starting y position.
      TrackData_t<Float_t> trkstartz;     // starting z position.
      TrackData_t<Float_t> trkstartd;     // starting distance to boundary.
      TrackData_t<Float_t> trkendx;       // ending x position.
      TrackData_t<Float_t> trkendy;       // ending y position.
      TrackData_t<Float_t> trkendz;       // ending z position.
      TrackData_t<Float_t> trkendd;       // ending distance to boundary.
      TrackData_t<Float_t> trktheta;      // theta.
      TrackData_t<Float_t> trkphi;        // phi.
      TrackData_t<Float_t> trkstartdcosx;
      TrackData_t<Float_t> trkstartdcosy;
      TrackData_t<Float_t> trkstartdcosz;
      TrackData_t<Float_t> trkenddcosx;
      TrackData_t<Float_t> trkenddcosy;
      TrackData_t<Float_t> trkenddcosz;
      TrackData_t<Float_t> trkthetaxz;    // theta_xz.
      TrackData_t<Float_t> trkthetayz;    // theta_yz.
      TrackData_t<Float_t> trkmom;        // momentum.
      TrackData_t<Float_t> trklen;        // length.
      TrackData_t<Float_t> trkmomrange;    // track momentum from range using CSDA tables
      TrackData_t<Float_t> trkmommschi2;   // track momentum from multiple scattering Chi2 method
      TrackData_t<Float_t> trkmommsllhd;   // track momentum from multiple scattering LLHD method
      TrackData_t<Short_t> trksvtxid;     // Vertex ID associated with the track start
      TrackData_t<Short_t> trkevtxid;     // Vertex ID associated with the track end
      PlaneData_t<Int_t> trkpidpdg;       // particle PID pdg code
      PlaneData_t<Float_t> trkpidchi;
      PlaneData_t<Float_t> trkpidchipr;   // particle PID chisq for proton
      PlaneData_t<Float_t> trkpidchika;   // particle PID chisq for kaon
      PlaneData_t<Float_t> trkpidchipi;   // particle PID chisq for pion
      PlaneData_t<Float_t> trkpidchimu;   // particle PID chisq for muon
      PlaneData_t<Float_t> trkpidpida;    // particle PIDA
      TrackData_t<Short_t> trkpidbestplane; // this is defined as the plane with most hits     
       
      /// Creates an empty tracker data structure
      TrackDataStruct(): MaxTracks(0) { Clear(); }
      /// Creates a tracker data structure allowing up to maxTracks tracks
      TrackDataStruct(size_t maxTracks): MaxTracks(maxTracks) { Clear(); }
      void Clear();
      void SetMaxTracks(size_t maxTracks)
        { MaxTracks = maxTracks; Resize(MaxTracks); }
      void Resize(size_t nTracks);
      void SetAddresses(TTree* pTree, std::string tracker);
      
      size_t GetMaxTracks() const { return MaxTracks; }
      size_t GetMaxPlanesPerTrack(int /* iTrack */ = 0) const
        { return (size_t) kNplanes; }
      size_t GetMaxHitsPerTrack(int /* iTrack */ = 0, int /* ipl */ = 0) const
        { return (size_t) kMaxTrackHits; }
      
    }; // class TrackDataStruct
    
 
    enum DataBits_t: unsigned int {
      tdAuxDet = 0x01,
      tdCry = 0x02,
      tdGenie = 0x04,
      tdGeant = 0x08,
      tdHit = 0x10,
      tdTrack = 0x20,
      tdVtx = 0x40,
      tdDefault = 0
    }; // DataBits_t
    
/*    /// information from the run
    struct RunData_t {
        public:
      RunData_t() { Clear(); }
      void Clear() {}
    }; // struct RunData_t
*/
    /// information from the subrun
    struct SubRunData_t {
      SubRunData_t() { Clear(); }
      void Clear() { pot = -99999.; }
      Double_t pot; //protons on target
    }; // struct SubRunData_t

//    RunData_t    RunData; ///< run data collected at begin of run
    SubRunData_t SubRunData; ///< subrun data collected at begin of subrun

    //run information
    Int_t      run;                  //run number
    Int_t      subrun;               //subrun number
    Int_t      event;                //event number
    Double_t   evttime;              //event time in sec
    Double_t   beamtime;             //beam time
  //  Double_t   pot;                  //protons on target moved in subrun data
    Double_t   taulife;              //electron lifetime
    Char_t     isdata;               //flag, 0=MC 1=data

    // hit information (non-resizeable, 45x kMaxHits = 900k bytes worth)
    Int_t    no_hits;                  //number of hits
    Short_t  hit_plane[kMaxHits];      //plane number
    Short_t  hit_wire[kMaxHits];       //wire number
    Short_t  hit_channel[kMaxHits];    //channel ID
    Float_t  hit_peakT[kMaxHits];      //peak time
    Float_t  hit_charge[kMaxHits];     //charge (area)
    Float_t  hit_ph[kMaxHits];         //amplitude
    Float_t  hit_startT[kMaxHits];     //hit start time
    Float_t  hit_endT[kMaxHits];       //hit end time
    Float_t  hit_nelec[kMaxHits];     //hit number of electrons
    Float_t  hit_energy[kMaxHits];       //hit energy
    Short_t  hit_trkid[kMaxHits];      //is this hit associated with a reco track?

    // vertex information
    Short_t  nvtx;                     //number of vertices
    Float_t  vtx[kMaxVertices][3];     //vtx[3]  

    //track information
    Char_t   kNTracker;
    std::vector<TrackDataStruct> TrackData;
    
    //mctruth information
    Int_t     mcevts_truth;    //number of neutrino Int_teractions in the spill
    Int_t     nuPDG_truth;     //neutrino PDG code
    Int_t     ccnc_truth;      //0=CC 1=NC
    Int_t     mode_truth;      //0=QE/El, 1=RES, 2=DIS, 3=Coherent production
    Float_t  enu_truth;       //true neutrino energy
    Float_t  Q2_truth;        //Momentum transfer squared
    Float_t  W_truth;         //hadronic invariant mass
    Int_t     hitnuc_truth;    //hit nucleon
    Float_t  nuvtxx_truth;    //neutrino vertex x
    Float_t  nuvtxy_truth;    //neutrino vertex y
    Float_t  nuvtxz_truth;    //neutrino vertex z
    Float_t  nu_dcosx_truth;  //neutrino dcos x
    Float_t  nu_dcosy_truth;  //neutrino dcos y
    Float_t  nu_dcosz_truth;  //neutrino dcos z
    Float_t  lep_mom_truth;   //lepton momentum
    Float_t  lep_dcosx_truth; //lepton dcos x
    Float_t  lep_dcosy_truth; //lepton dcos y
    Float_t  lep_dcosz_truth; //lepton dcos z

    //flux information
    Float_t  tpx_flux;        //Px of parent particle leaving BNB target
    Float_t  tpy_flux;        //Py of parent particle leaving BNB target
    Float_t  tpz_flux;        //Pz of parent particle leaving BNB target
    Int_t     tptype_flux;     //Type of parent particle leaving BNB target

    //genie information
    size_t MaxGeniePrimaries = 0;
    Int_t     genie_no_primaries;
    std::vector<Int_t>     genie_primaries_pdg;
    std::vector<Float_t>  genie_Eng;
    std::vector<Float_t>  genie_Px;
    std::vector<Float_t>  genie_Py;
    std::vector<Float_t>  genie_Pz;
    std::vector<Float_t>  genie_P;
    std::vector<Int_t>     genie_status_code;
    std::vector<Float_t>  genie_mass;
    std::vector<Int_t>     genie_trackID;
    std::vector<Int_t>     genie_ND;
    std::vector<Int_t>     genie_mother;
    
    //cosmic cry information
    Int_t     mcevts_truthcry;    //number of neutrino Int_teractions in the spill
    Int_t     cry_no_primaries;
    std::vector<Int_t>    cry_primaries_pdg;
    std::vector<Float_t>  cry_Eng;
    std::vector<Float_t>  cry_Px;
    std::vector<Float_t>  cry_Py;
    std::vector<Float_t>  cry_Pz;
    std::vector<Float_t>  cry_P;
    std::vector<Float_t>  cry_StartPointx;
    std::vector<Float_t>  cry_StartPointy;
    std::vector<Float_t>  cry_StartPointz;
    std::vector<Int_t>    cry_status_code;
    std::vector<Float_t>  cry_mass;
    std::vector<Int_t>    cry_trackID;
    std::vector<Int_t>    cry_ND;
    std::vector<Int_t>    cry_mother;
    
    //geant information
    size_t MaxGEANTparticles = 0; ///! how many particles there is currently room for
    Int_t     no_primaries;      //number of primary geant particles
    Int_t     geant_list_size;  //number of all geant particles
    Int_t     geant_list_size_in_tpcAV;
    std::vector<Int_t>    pdg;
    std::vector<Int_t>    status;    
    std::vector<Float_t>  Eng;
    std::vector<Float_t>  EndE;
    std::vector<Float_t>  Mass;
    std::vector<Float_t>  Px;
    std::vector<Float_t>  Py;
    std::vector<Float_t>  Pz;
    std::vector<Float_t>  P;
    std::vector<Float_t>  StartPointx;
    std::vector<Float_t>  StartPointy;
    std::vector<Float_t>  StartPointz;
    std::vector<Float_t>  StartT;  
    std::vector<Float_t>  EndT;          
    std::vector<Float_t>  EndPointx;
    std::vector<Float_t>  EndPointy;
    std::vector<Float_t>  EndPointz;
    std::vector<Float_t>  theta;    
    std::vector<Float_t>  phi;    
    std::vector<Float_t>  theta_xz;    
    std::vector<Float_t>  theta_yz;    
    std::vector<Float_t>  pathlen;    
    std::vector<Int_t>    inTPCActive;    
    std::vector<Float_t>  StartPointx_tpcAV;
    std::vector<Float_t>  StartPointy_tpcAV;
    std::vector<Float_t>  StartPointz_tpcAV;
    std::vector<Float_t>  EndPointx_tpcAV;
    std::vector<Float_t>  EndPointy_tpcAV;
    std::vector<Float_t>  EndPointz_tpcAV;
    std::vector<Int_t>    NumberDaughters;
    std::vector<Int_t>    TrackId;
    std::vector<Int_t>    Mother;
    std::vector<Int_t>    process_primary;
    std::vector<std::string> processname;
    std::vector<Int_t>    MergedId; //geant track segments, which belong to the same particle, get the same
    
    // Auxiliary detector variables saved for each geant track
    // This data is saved as a vector (one item per GEANT particle) of C arrays
    // (wrapped in a BoxedArray for technical reasons), one item for each
    // affected detector cell (which one is saved in AuxDetID
    template <typename T>
    using AuxDetMCData_t = std::vector<BoxedArray<T[kMaxAuxDets]>>;
    
    std::vector<UShort_t> NAuxDets;         ///< Number of AuxDets crossed by this particle
    AuxDetMCData_t<Short_t> AuxDetID;       ///< Which AuxDet this particle went through
    AuxDetMCData_t<Float_t> entryX;         ///< Entry X position of particle into AuxDet
    AuxDetMCData_t<Float_t> entryY;         ///< Entry Y position of particle into AuxDet
    AuxDetMCData_t<Float_t> entryZ;         ///< Entry Z position of particle into AuxDet
    AuxDetMCData_t<Float_t> entryT;         ///< Entry T position of particle into AuxDet
    AuxDetMCData_t<Float_t> exitX;          ///< Exit X position of particle out of AuxDet
    AuxDetMCData_t<Float_t> exitY;          ///< Exit Y position of particle out of AuxDet
    AuxDetMCData_t<Float_t> exitZ;          ///< Exit Z position of particle out of AuxDet
    AuxDetMCData_t<Float_t> exitT;          ///< Exit T position of particle out of AuxDet
    AuxDetMCData_t<Float_t> exitPx;         ///< Exit x momentum of particle out of AuxDet
    AuxDetMCData_t<Float_t> exitPy;         ///< Exit y momentum of particle out of AuxDet
    AuxDetMCData_t<Float_t> exitPz;         ///< Exit z momentum of particle out of AuxDet
    AuxDetMCData_t<Float_t> CombinedEnergyDep; ///< Sum energy of all particles with this trackID (+ID or -ID) in AuxDet
    
    unsigned int bits; ///< complementary information
 
    /// Returns whether we have auxiliary detector data
    bool hasAuxDetector() const { return bits & tdAuxDet; }
    
    /// Returns whether we have Cry data
    bool hasCryInfo() const { return bits & tdCry; }
    
    /// Returns whether we have Genie data
    bool hasGenieInfo() const { return bits & tdGenie; }
    
    /// Returns whether we have Hit data
    bool hasHitInfo() const { return bits & tdHit; }

    /// Returns whether we have Track data
    bool hasTrackInfo() const { return bits & tdTrack; }
    
    /// Returns whether we have Vertex data
    bool hasVertexInfo() const { return bits & tdVtx; }
    
    /// Returns whether we have Geant data
    bool hasGeantInfo() const { return bits & tdGeant; }

    /// Sets the specified bits
    void SetBits(unsigned int setbits, bool unset = false)
      { if (unset) bits &= ~setbits; else bits |= setbits; }
      
    /// Constructor; clears all fields
    AnalysisTreeDataStruct(size_t nTrackers = 0): bits(tdDefault) 
      { SetTrackers(nTrackers); Clear(); }

    TrackDataStruct& GetTrackerData(size_t iTracker)
      { return TrackData.at(iTracker); }
    const TrackDataStruct& GetTrackerData(size_t iTracker) const
      { return TrackData.at(iTracker); }
    
    
    /// Clear all fields if this object (not the tracker algorithm data)
    void ClearLocalData();
    
    /// Clear all fields
    void Clear();
    
    
    /// Allocates data structures for the given number of trackers (no Clear())
    void SetTrackers(size_t nTrackers) { TrackData.resize(nTrackers); }
    
    /// Resize the data strutcure for GEANT particles
    void ResizeGEANT(int nParticles);
    
    /// Resize the data strutcure for Genie primaries
    void ResizeGenie(int nPrimaries);
    
    /// Resize the data strutcure for Cry primaries
    void ResizeCry(int nPrimaries);
    
    /// Connect this object with a tree
    void SetAddresses(TTree* pTree, const std::vector<std::string>& trackers);
    
    
    /// Returns the number of trackers for which data structures are allocated
    size_t GetNTrackers() const { return TrackData.size(); }
    
    /// Returns the number of hits for which memory is allocated
    size_t GetMaxHits() const { return kMaxHits; }
    
    /// Returns the number of trackers for which memory is allocated
    size_t GetMaxTrackers() const { return TrackData.capacity(); }
    
    /// Returns the number of GEANT particles for which memory is allocated
    size_t GetMaxGEANTparticles() const { return MaxGEANTparticles; }
        
    /// Returns the number of GENIE primaries for which memory is allocated
    size_t GetMaxGeniePrimaries() const { return MaxGeniePrimaries; }
    
    
  private:
    /// Little helper functor class to create or reset branches in a tree
    class BranchCreator {
        public:
      TTree* pTree; ///< the tree to be worked on
      BranchCreator(TTree* tree): pTree(tree) {}

      //@{
      /// Create a branch if it does not exist, and set its address
      void operator()
        (std::string name, void* address, std::string leaflist /*, int bufsize = 32000 */)
        {
          if (!pTree) return;
          TBranch* pBranch = pTree->GetBranch(name.c_str());
          if (!pBranch) {
            pTree->Branch(name.c_str(), address, leaflist.c_str() /*, bufsize */);
            LOG_DEBUG("AnalysisTreeStructure")
              << "Creating branch '" << name << " with leaf '" << leaflist << "'";
          }
          else if (pBranch->GetAddress() != address) {
            pBranch->SetAddress(address);
            LOG_DEBUG("AnalysisTreeStructure")
              << "Reassigning address to branch '" << name << "'";
          }
          else {
            LOG_DEBUG("AnalysisTreeStructure")
              << "Branch '" << name << "' is fine";
          }
        } // operator()
      void operator()
        (std::string name, void* address, const std::stringstream& leaflist /*, int bufsize = 32000 */)
        { return this->operator() (name, address, leaflist.str() /*, int bufsize = 32000 */); }
      template <typename T>
      void operator()
        (std::string name, std::vector<T>& data, std::string leaflist /*, int bufsize = 32000 */)
        { return this->operator() (name, (void*) data.data(), leaflist /*, int bufsize = 32000 */); }

      template <typename T>
      void operator() (std::string name, std::vector<T>& data)
        {
          // overload for a generic object expressed directly by reference
          // (as opposed to a generic object expressed by a pointer or
          // to a simple leaf sequence specification);
          // TTree::Branch(name, T* obj, Int_t bufsize, splitlevel) and
          // TTree::SetObject() are used.
          if (!pTree) return;
          TBranch* pBranch = pTree->GetBranch(name.c_str());
          if (!pBranch) {
            pTree->Branch(name.c_str(), &data);
            // ROOT needs a TClass definition for T in order to create a branch,
            // se we are sure that at this point the TClass exists
            LOG_DEBUG("AnalysisTreeStructure")
              << "Creating object branch '" << name
              << " with " << TClass::GetClass(typeid(T))->ClassName();
          }
          else if
            (*(reinterpret_cast<std::vector<T>**>(pBranch->GetAddress())) != &data)
          {
            // when an object is provided directly, the address of the object
            // is assigned in TBranchElement::fObject (via TObject::SetObject())
            // and the address itself is set to the address of the fObject
            // member. Here we check that the address of the object in fObject
            // is the same as the address of our current data type
            pBranch->SetObject(&data);
            LOG_DEBUG("AnalysisTreeStructure")
              << "Reassigning object to branch '" << name << "'";
          }
          else {
            LOG_DEBUG("AnalysisTreeStructure")
              << "Branch '" << name << "' is fine";
          }
        } // operator()
      //@}
    }; // class BranchCreator

  }; // class AnalysisTreeDataStruct
  
  
  /**
   * @brief Creates a simple ROOT tree with tracking and calorimetry information
   * 
   * <h2>Configuration parameters</h2>
   * - <b>UseBuffers</b> (default: false): if enabled, memory is allocated for
   *   tree data for all the run; otherwise, it's allocated on each event, used
   *   and freed; use "true" for speed, "false" to save memory
   * - <b>SaveAuxDetInfo</b> (default: false): if enabled, auxiliary detector
   *   data will be extracted and included in the tree
   */
  class AnalysisTree : public art::EDAnalyzer {

  public:

    explicit AnalysisTree(fhicl::ParameterSet const& pset);
    virtual ~AnalysisTree();

    /// read access to event
    void analyze(const art::Event& evt);
  //  void beginJob() {}
    void beginSubRun(const art::SubRun& sr);

  private:

    void   HitsPurity(std::vector< art::Ptr<recob::Hit> > const& hits, Int_t& trackid, Float_t& purity, double& maxe);
    double length(const recob::Track& track);
    double length(const simb::MCParticle& part, TVector3& start, TVector3& end);
    double bdist(const TVector3& pos);

    TTree* fTree;

    // event information is huge and dynamic;
    // run information is much smaller and we still store it statically
    // in the event
    AnalysisTreeDataStruct* fData;
//    AnalysisTreeDataStruct::RunData_t RunData;
    AnalysisTreeDataStruct::SubRunData_t SubRunData;

    std::string fDigitModuleLabel;
    std::string fHitsModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fCalDataModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fCryGenModuleLabel;
    std::string fG4ModuleLabel;
    std::string fVertexModuleLabel;
    std::vector<std::string> fTrackModuleLabel;
    std::vector<std::string> fCalorimetryModuleLabel;
    std::vector<std::string> fParticleIDModuleLabel;
    std::string fPOTModuleLabel;
    bool fUseBuffer; ///< whether to use a permanent buffer (faster, huge memory)    
    bool fSaveAuxDetInfo; ///< whether to extract and save auxiliary detector data
    bool fSaveCryInfo; ///whether to extract and save CRY particle data
    bool fSaveGenieInfo; ///whether to extract and save Genie information
    bool fSaveGeantInfo; ///whether to extract and save Geant information
    bool fSaveHitInfo; ///whether to extract and save Hit information
    bool fSaveTrackInfo; ///whether to extract and save Track information
    bool fSaveVertexInfo; ///whether to extract and save Vertex information
    
    std::vector<std::string> fCosmicTaggerAssocLabel;
    std::vector<std::string> fFlashMatchAssocLabel;

    /// Returns the number of trackers configured
    size_t GetNTrackers() const { return fTrackModuleLabel.size(); }
       
    /// Creates the structure for the tree data; optionally initializes it
    void CreateData(bool bClearData = false)
      {
        if (!fData) {
          fData = new AnalysisTreeDataStruct(GetNTrackers());
          fData->SetBits(AnalysisTreeDataStruct::tdAuxDet, !fSaveAuxDetInfo);
	  fData->SetBits(AnalysisTreeDataStruct::tdCry, !fSaveCryInfo);	  
	  fData->SetBits(AnalysisTreeDataStruct::tdGenie, !fSaveGenieInfo);
	  fData->SetBits(AnalysisTreeDataStruct::tdGeant, !fSaveGeantInfo); 
        }
        else {
          fData->SetBits(AnalysisTreeDataStruct::tdHit, !fSaveHitInfo);	
	  fData->SetBits(AnalysisTreeDataStruct::tdTrack, !fSaveTrackInfo);	
	  fData->SetBits(AnalysisTreeDataStruct::tdVtx, !fSaveVertexInfo);	  	  	    	  	    	  	    	  
	  fData->SetTrackers(GetNTrackers());
          if (bClearData) fData->Clear();
        }
      } // CreateData()
    
    /// Sets the addresses of all the tree branches, creating the missing ones
    void SetAddresses()
      {
        CheckData("SetAddress()"); CheckTree("SetAddress()");
        fData->SetAddresses(fTree, fTrackModuleLabel);
      } // SetAddresses()
    
    /// Sets the addresses of all the tree branches of the specified tracking algo,
    /// creating the missing ones
    void SetTrackerAddresses(size_t iTracker)
      {
        CheckData("SetTrackerAddresses()"); CheckTree("SetTrackerAddresses()");
        if (iTracker >= fData->GetNTrackers()) {
          throw art::Exception(art::errors::LogicError)
            << "AnalysisTree::SetTrackerAddresses(): no tracker #" << iTracker
            << " (" << fData->GetNTrackers() << " available)";
        }
        fData->GetTrackerData(iTracker) \
          .SetAddresses(fTree, fTrackModuleLabel[iTracker]);
      } // SetTrackerAddresses()
    
    /// Create the output tree and the data structures, if needed
    void CreateTree(bool bClearData = false);
    
    /// Destroy the local buffers (existing branches will point to invalid address!)
    void DestroyData() { if (fData) { delete fData; fData = nullptr; } }
    
    /// Helper function: throws if no data structure is available
    void CheckData(std::string caller) const
      {
        if (fData) return;
        throw art::Exception(art::errors::LogicError)
          << "AnalysisTree::" << caller << ": no data";
      } // CheckData()
    /// Helper function: throws if no tree is available
    void CheckTree(std::string caller) const
      {
        if (fTree) return;
        throw art::Exception(art::errors::LogicError)
          << "AnalysisTree::" << caller << ": no tree";
      } // CheckData()
  }; // class microboone::AnalysisTree
} // namespace microboone


namespace { // local namespace
  /// Simple stringstream which empties its buffer on operator() call
  class AutoResettingStringSteam: public std::ostringstream {
      public:
    AutoResettingStringSteam& operator() () { str(""); return *this; }
  }; // class AutoResettingStringSteam

  /// Fills a sequence of TYPE elements
  template <typename ITER, typename TYPE>
  inline void FillWith(ITER from, ITER to, TYPE value)
    { std::fill(from, to, value); }

  /// Fills a sequence of TYPE elements
  template <typename ITER, typename TYPE>
  inline void FillWith(ITER from, size_t n, TYPE value)
    { std::fill(from, from + n, value); }

  /// Fills a container with begin()/end() interface
  template <typename CONT, typename V>
  inline void FillWith(CONT& data, const V& value)
    { FillWith(std::begin(data), std::end(data), value); }

} // local namespace


//------------------------------------------------------------------------------
//---  AnalysisTreeDataStruct::TrackDataStruct
//---

void microboone::AnalysisTreeDataStruct::TrackDataStruct::Resize(size_t nTracks)
{
  MaxTracks = nTracks;
  
  trkId.resize(MaxTracks);
  trkncosmictags_tagger.resize(MaxTracks);
  trkcosmicscore_tagger.resize(MaxTracks);
  trkcosmictype_tagger.resize(MaxTracks);
  trkncosmictags_flashmatch.resize(MaxTracks);
  trkcosmicscore_flashmatch.resize(MaxTracks);
  trkcosmictype_flashmatch.resize(MaxTracks);
  trkstartx.resize(MaxTracks);
  trkstarty.resize(MaxTracks);
  trkstartz.resize(MaxTracks);
  trkstartd.resize(MaxTracks);
  trkendx.resize(MaxTracks);
  trkendy.resize(MaxTracks);
  trkendz.resize(MaxTracks);
  trkendd.resize(MaxTracks);
  trktheta.resize(MaxTracks);
  trkphi.resize(MaxTracks);
  trkstartdcosx.resize(MaxTracks);
  trkstartdcosy.resize(MaxTracks);
  trkstartdcosz.resize(MaxTracks);
  trkenddcosx.resize(MaxTracks);
  trkenddcosy.resize(MaxTracks);
  trkenddcosz.resize(MaxTracks);
  trkthetaxz.resize(MaxTracks);
  trkthetayz.resize(MaxTracks);
  trkmom.resize(MaxTracks);
  trkmomrange.resize(MaxTracks);
  trkmommschi2.resize(MaxTracks);
  trkmommsllhd.resize(MaxTracks);  
  trklen.resize(MaxTracks);
  trksvtxid.resize(MaxTracks);
  trkevtxid.resize(MaxTracks);
  // PID variables
  trkpidpdg.resize(MaxTracks);
  trkpidchi.resize(MaxTracks);
  trkpidchipr.resize(MaxTracks);
  trkpidchika.resize(MaxTracks);
  trkpidchipi.resize(MaxTracks);
  trkpidchimu.resize(MaxTracks);
  trkpidpida.resize(MaxTracks);
  trkpidbestplane.resize(MaxTracks);
  
  trkke.resize(MaxTracks);
  trkrange.resize(MaxTracks);
  trkidtruth.resize(MaxTracks);
  trkorigin.resize(MaxTracks);
  trkpdgtruth.resize(MaxTracks);
  trkefftruth.resize(MaxTracks);
  trksimIDEenergytruth.resize(MaxTracks);
  trksimIDExtruth.resize(MaxTracks);
  trksimIDEytruth.resize(MaxTracks);
  trksimIDEztruth.resize(MaxTracks);
  trkpurtruth.resize(MaxTracks);
  trkpitchc.resize(MaxTracks);
  ntrkhits.resize(MaxTracks);
  
  trkdedx.resize(MaxTracks);
  trkdqdx.resize(MaxTracks);
  trkresrg.resize(MaxTracks);
  trkxyz.resize(MaxTracks);
  
} // microboone::AnalysisTreeDataStruct::TrackDataStruct::Resize()

void microboone::AnalysisTreeDataStruct::TrackDataStruct::Clear() {
  Resize(MaxTracks);
  ntracks = 0;
  
  FillWith(trkId        , -9999  );
  FillWith(trkncosmictags_tagger, -9999  );
  FillWith(trkcosmicscore_tagger, -99999.);
  FillWith(trkcosmictype_tagger, -9999  );
  FillWith(trkncosmictags_flashmatch, -9999  );
  FillWith(trkcosmicscore_flashmatch, -99999.);
  FillWith(trkcosmictype_flashmatch, -9999  );
  FillWith(trkstartx    , -99999.);
  FillWith(trkstarty    , -99999.);
  FillWith(trkstartz    , -99999.);
  FillWith(trkstartd    , -99999.);
  FillWith(trkendx      , -99999.);
  FillWith(trkendy      , -99999.);
  FillWith(trkendz      , -99999.);
  FillWith(trkendd      , -99999.);
  FillWith(trktheta     , -99999.);
  FillWith(trkphi       , -99999.);
  FillWith(trkstartdcosx, -99999.);
  FillWith(trkstartdcosy, -99999.);
  FillWith(trkstartdcosz, -99999.);
  FillWith(trkenddcosx  , -99999.);
  FillWith(trkenddcosy  , -99999.);
  FillWith(trkenddcosz  , -99999.);
  FillWith(trkthetaxz   , -99999.);
  FillWith(trkthetayz   , -99999.);
  FillWith(trkmom       , -99999.);
  FillWith(trkmomrange  , -99999.);  
  FillWith(trkmommschi2 , -99999.);  
  FillWith(trkmommsllhd , -99999.);  
  FillWith(trklen       , -99999.);
  FillWith(trksvtxid    , -1);
  FillWith(trkevtxid    , -1);
  FillWith(trkpidbestplane, -1); 
 
  for (size_t iTrk = 0; iTrk < MaxTracks; ++iTrk){
    
    // the following are BoxedArray's;
    // their iterators traverse all the array dimensions
    FillWith(trkke[iTrk]      , -99999.);
    FillWith(trkrange[iTrk]   , -99999.);
    FillWith(trkidtruth[iTrk] , -99999 );
    FillWith(trkorigin[iTrk]  , -1 );
    FillWith(trkpdgtruth[iTrk], -99999 );
    FillWith(trkefftruth[iTrk], -99999.);
    FillWith(trksimIDEenergytruth[iTrk], -99999.);
    FillWith(trksimIDExtruth[iTrk], -99999.);
    FillWith(trksimIDEytruth[iTrk], -99999.);
    FillWith(trksimIDEztruth[iTrk], -99999.);
    FillWith(trkpurtruth[iTrk], -99999.);
    FillWith(trkpitchc[iTrk]  , -99999.);
    FillWith(ntrkhits[iTrk]   ,  -9999 );
    
    FillWith(trkdedx[iTrk], 0.);
    FillWith(trkdqdx[iTrk], 0.);
    FillWith(trkresrg[iTrk], 0.);
    
    FillWith(trkxyz[iTrk], 0.);
 
    FillWith(trkpidpdg[iTrk]    , -1);
    FillWith(trkpidchi[iTrk]    , -99999.);
    FillWith(trkpidchipr[iTrk]  , -99999.);
    FillWith(trkpidchika[iTrk]  , -99999.);
    FillWith(trkpidchipi[iTrk]  , -99999.);
    FillWith(trkpidchimu[iTrk]  , -99999.);
    FillWith(trkpidpida[iTrk]   , -99999.);
  } // for track
  
} // microboone::AnalysisTreeDataStruct::TrackDataStruct::Clear()


void microboone::AnalysisTreeDataStruct::TrackDataStruct::SetAddresses(
  TTree* pTree, std::string tracker
) {
  if (MaxTracks == 0) return; // no tracks, no tree!
  
  microboone::AnalysisTreeDataStruct::BranchCreator CreateBranch(pTree);

  AutoResettingStringSteam sstr;
  sstr() << kMaxTrackHits;
  std::string MaxTrackHitsIndexStr("[" + sstr.str() + "]");
  
  std::string TrackLabel = tracker;
  std::string BranchName;

  BranchName = "ntracks_" + TrackLabel;
  CreateBranch(BranchName, &ntracks, BranchName + "/S");
  std::string NTracksIndexStr = "[" + BranchName + "]";
  
  BranchName = "trkId_" + TrackLabel;
  CreateBranch(BranchName, trkId, BranchName + NTracksIndexStr + "/S");
  
  BranchName = "trkncosmictags_tagger_" + TrackLabel;
  CreateBranch(BranchName, trkncosmictags_tagger, BranchName + NTracksIndexStr + "/S");
  
  BranchName = "trkcosmicscore_tagger_" + TrackLabel;
  CreateBranch(BranchName, trkcosmicscore_tagger, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkcosmictype_tagger_" + TrackLabel;
  CreateBranch(BranchName, trkcosmictype_tagger, BranchName + NTracksIndexStr + "/S");

  BranchName = "trkncosmictags_flashmatch_" + TrackLabel;
  CreateBranch(BranchName, trkncosmictags_flashmatch, BranchName + NTracksIndexStr + "/S");
  
  BranchName = "trkcosmicscore_flashmatch_" + TrackLabel;
  CreateBranch(BranchName, trkcosmicscore_flashmatch, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkcosmictype_flashmatch_" + TrackLabel;
  CreateBranch(BranchName, trkcosmictype_flashmatch, BranchName + NTracksIndexStr + "/S");
  
  BranchName = "trkke_" + TrackLabel;
  CreateBranch(BranchName, trkke, BranchName + NTracksIndexStr + "[3]/F");
  
  BranchName = "trkrange_" + TrackLabel;
  CreateBranch(BranchName, trkrange, BranchName + NTracksIndexStr + "[3]/F");
   
  BranchName = "trkidtruth_" + TrackLabel;
  CreateBranch(BranchName, trkidtruth, BranchName + NTracksIndexStr + "[3]/I");

  BranchName = "trkorigin_" + TrackLabel;
  CreateBranch(BranchName, trkorigin, BranchName + NTracksIndexStr + "[3]/S");
  
  BranchName = "trkpdgtruth_" + TrackLabel;
  CreateBranch(BranchName, trkpdgtruth, BranchName + NTracksIndexStr + "[3]/I");
  
  BranchName = "trkefftruth_" + TrackLabel;
  CreateBranch(BranchName, trkefftruth, BranchName + NTracksIndexStr + "[3]/F");
 
  BranchName = "trksimIDEenergytruth_" + TrackLabel;
  CreateBranch(BranchName, trksimIDEenergytruth, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trksimIDExtruth_" + TrackLabel;
  CreateBranch(BranchName, trksimIDExtruth, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trksimIDEytruth_" + TrackLabel;
  CreateBranch(BranchName, trksimIDEytruth, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trksimIDEztruth_" + TrackLabel;
  CreateBranch(BranchName, trksimIDEztruth, BranchName + NTracksIndexStr + "[3]/F");
 
  BranchName = "trkpurtruth_" + TrackLabel;
  CreateBranch(BranchName, trkpurtruth, BranchName + NTracksIndexStr + "[3]/F");
  
  BranchName = "trkpitchc_" + TrackLabel;
  CreateBranch(BranchName, trkpitchc, BranchName + NTracksIndexStr + "[3]/F");
  
  BranchName = "ntrkhits_" + TrackLabel;
  CreateBranch(BranchName, ntrkhits, BranchName + NTracksIndexStr + "[3]/S");
  
  BranchName = "trkdedx_" + TrackLabel;
  CreateBranch(BranchName, trkdedx, BranchName + NTracksIndexStr + "[3]" + MaxTrackHitsIndexStr + "/F");
  
  BranchName = "trkdqdx_" + TrackLabel;
  CreateBranch(BranchName, trkdqdx, BranchName + NTracksIndexStr + "[3]" + MaxTrackHitsIndexStr + "/F");
  
  BranchName = "trkresrg_" + TrackLabel;
  CreateBranch(BranchName, trkresrg, BranchName + NTracksIndexStr + "[3]" + MaxTrackHitsIndexStr + "/F");
  
  BranchName = "trkxyz_" + TrackLabel;
  CreateBranch(BranchName, trkxyz, BranchName + NTracksIndexStr + "[3]" + MaxTrackHitsIndexStr + "/F");
  
  BranchName = "trkstartx_" + TrackLabel;
  CreateBranch(BranchName, trkstartx, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkstarty_" + TrackLabel;
  CreateBranch(BranchName, trkstarty, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkstartz_" + TrackLabel;
  CreateBranch(BranchName, trkstartz, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkstartd_" + TrackLabel;
  CreateBranch(BranchName, trkstartd, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkendx_" + TrackLabel;
  CreateBranch(BranchName, trkendx, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkendy_" + TrackLabel;
  CreateBranch(BranchName, trkendy, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkendz_" + TrackLabel;
  CreateBranch(BranchName, trkendz, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkendd_" + TrackLabel;
  CreateBranch(BranchName, trkendd, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trktheta_" + TrackLabel;
  CreateBranch(BranchName, trktheta, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkphi_" + TrackLabel;
  CreateBranch(BranchName, trkphi, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkstartdcosx_" + TrackLabel;
  CreateBranch(BranchName, trkstartdcosx, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkstartdcosy_" + TrackLabel;
  CreateBranch(BranchName, trkstartdcosy, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkstartdcosz_" + TrackLabel;
  CreateBranch(BranchName, trkstartdcosz, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkenddcosx_" + TrackLabel;
  CreateBranch(BranchName, trkenddcosx, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkenddcosy_" + TrackLabel;
  CreateBranch(BranchName, trkenddcosy, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkenddcosz_" + TrackLabel;
  CreateBranch(BranchName, trkenddcosz, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkthetaxz_" + TrackLabel;
  CreateBranch(BranchName, trkthetaxz, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkthetayz_" + TrackLabel;
  CreateBranch(BranchName, trkthetayz, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkmom_" + TrackLabel;
  CreateBranch(BranchName, trkmom, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkmomrange_" + TrackLabel;
  CreateBranch(BranchName, trkmomrange, BranchName + NTracksIndexStr + "/F");

  BranchName = "trkmommschi2_" + TrackLabel;
  CreateBranch(BranchName, trkmommschi2, BranchName + NTracksIndexStr + "/F");

  BranchName = "trkmommsllhd_" + TrackLabel;
  CreateBranch(BranchName, trkmommsllhd, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trklen_" + TrackLabel;
  CreateBranch(BranchName, trklen, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trksvtxid_" + TrackLabel;
  CreateBranch(BranchName, trksvtxid, BranchName + NTracksIndexStr + "/S");
  
  BranchName = "trkevtxid_" + TrackLabel;
  CreateBranch(BranchName, trkevtxid, BranchName + NTracksIndexStr + "/S");

  BranchName = "trkpidpdg_" + TrackLabel;
  CreateBranch(BranchName, trkpidpdg, BranchName + NTracksIndexStr + "[3]/I");

  BranchName = "trkpidchi_" + TrackLabel;
  CreateBranch(BranchName, trkpidchi, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trkpidchipr_" + TrackLabel;
  CreateBranch(BranchName, trkpidchipr, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trkpidchika_" + TrackLabel;
  CreateBranch(BranchName, trkpidchika, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trkpidchipi_" + TrackLabel;
  CreateBranch(BranchName, trkpidchipi, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trkpidchimu_" + TrackLabel;
  CreateBranch(BranchName, trkpidchimu, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trkpidpida_" + TrackLabel;
  CreateBranch(BranchName, trkpidpida, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trkpidbestplane_" + TrackLabel;
  CreateBranch(BranchName, trkpidbestplane, BranchName + NTracksIndexStr + "/S");

} // microboone::AnalysisTreeDataStruct::TrackDataStruct::SetAddresses()

//------------------------------------------------------------------------------
//---  AnalysisTreeDataStruct
//---

void microboone::AnalysisTreeDataStruct::ClearLocalData() {

//  RunData.Clear();
  SubRunData.Clear();

  run = -99999;
  subrun = -99999;
  event = -99999;
  evttime = -99999;
  beamtime = -99999;
  isdata = -99;
  taulife = -99999;

  no_hits = 0;
  
  std::fill(hit_plane, hit_plane + sizeof(hit_plane)/sizeof(hit_plane[0]), -9999);
  std::fill(hit_wire, hit_wire + sizeof(hit_wire)/sizeof(hit_wire[0]), -9999);
  std::fill(hit_channel, hit_channel + sizeof(hit_channel)/sizeof(hit_channel[0]), -9999);
  std::fill(hit_peakT, hit_peakT + sizeof(hit_peakT)/sizeof(hit_peakT[0]), -99999.);
  std::fill(hit_charge, hit_charge + sizeof(hit_charge)/sizeof(hit_charge[0]), -99999.);
  std::fill(hit_ph, hit_ph + sizeof(hit_ph)/sizeof(hit_ph[0]), -99999.);
  std::fill(hit_startT, hit_startT + sizeof(hit_startT)/sizeof(hit_startT[0]), -99999.);
  std::fill(hit_endT, hit_endT + sizeof(hit_endT)/sizeof(hit_endT[0]), -99999.);
  std::fill(hit_trkid, hit_trkid + sizeof(hit_trkid)/sizeof(hit_trkid[0]), -9999);
  std::fill(hit_nelec, hit_nelec + sizeof(hit_nelec)/sizeof(hit_nelec[0]), -99999.);
  std::fill(hit_energy, hit_energy + sizeof(hit_energy)/sizeof(hit_energy[0]), -99999.);

  nvtx = 0;
  for (size_t ivtx = 0; ivtx < kMaxVertices; ++ivtx) {
    std::fill(vtx[ivtx], vtx[ivtx]+3, -99999.);
  }

  mcevts_truth = -99999;
  mcevts_truthcry = -99999;
  nuPDG_truth = -99999;
  ccnc_truth = -99999;
  mode_truth = -99999;
  enu_truth = -99999;
  Q2_truth = -99999;
  W_truth = -99999;
  hitnuc_truth = -99999;
  nuvtxx_truth = -99999;
  nuvtxy_truth = -99999;
  nuvtxz_truth = -99999;
  nu_dcosx_truth = -99999;
  nu_dcosy_truth = -99999;
  nu_dcosz_truth = -99999;
  lep_mom_truth = -99999;
  lep_dcosx_truth = -99999;
  lep_dcosy_truth = -99999;
  lep_dcosz_truth = -99999;
  tpx_flux = -99999;
  tpy_flux = -99999;
  tpz_flux = -99999;
  tptype_flux = -99999;

  genie_no_primaries = 0;
  cry_no_primaries = 0;
  no_primaries = 0;
  geant_list_size=0;
  geant_list_size_in_tpcAV = 0;
  
  FillWith(pdg, -99999);
  FillWith(status, -99999);
  FillWith(Mass, -99999.);
  FillWith(Eng, -99999.);
  FillWith(EndE, -99999.);
  FillWith(Px, -99999.);
  FillWith(Py, -99999.);
  FillWith(Pz, -99999.);
  FillWith(P, -99999.);
  FillWith(StartPointx, -99999.);
  FillWith(StartPointy, -99999.);
  FillWith(StartPointz, -99999.);
  FillWith(StartT, -99999.);
  FillWith(EndT, -99999.);    
  FillWith(EndPointx, -99999.);
  FillWith(EndPointy, -99999.);
  FillWith(EndPointz, -99999.);
  FillWith(EndT, -99999.);
  FillWith(theta, -99999.);
  FillWith(phi, -99999.);
  FillWith(theta_xz, -99999.);
  FillWith(theta_yz, -99999.);
  FillWith(pathlen, -99999.);
  FillWith(inTPCActive, -99999);
  FillWith(StartPointx_tpcAV, -99999.);
  FillWith(StartPointy_tpcAV, -99999.);
  FillWith(StartPointz_tpcAV, -99999.);
  FillWith(EndPointx_tpcAV, -99999.);
  FillWith(EndPointy_tpcAV, -99999.);
  FillWith(EndPointz_tpcAV, -99999.);  
  FillWith(NumberDaughters, -99999);
  FillWith(Mother, -99999);
  FillWith(TrackId, -99999);
  FillWith(process_primary, -99999);
  FillWith(processname, "noname");
  FillWith(MergedId, -99999);
  FillWith(genie_primaries_pdg, -99999);
  FillWith(genie_Eng, -99999.);
  FillWith(genie_Px, -99999.);
  FillWith(genie_Py, -99999.);
  FillWith(genie_Pz, -99999.);
  FillWith(genie_P, -99999.);
  FillWith(genie_status_code, -99999);
  FillWith(genie_mass, -99999.);
  FillWith(genie_trackID, -99999);
  FillWith(genie_ND, -99999);
  FillWith(genie_mother, -99999);
  FillWith(cry_primaries_pdg, -99999);
  FillWith(cry_Eng, -99999.);
  FillWith(cry_Px, -99999.);
  FillWith(cry_Py, -99999.);
  FillWith(cry_Pz, -99999.);
  FillWith(cry_P, -99999.);
  FillWith(cry_StartPointx, -99999.);
  FillWith(cry_StartPointy, -99999.);
  FillWith(cry_StartPointz, -99999.);  
  FillWith(cry_status_code, -99999);
  FillWith(cry_mass, -99999.);
  FillWith(cry_trackID, -99999);
  FillWith(cry_ND, -99999);
  FillWith(cry_mother, -99999);
  
  // auxiliary detector information;
  FillWith(NAuxDets, 0);
  // - set to -9999 all the values of each of the arrays in AuxDetID;
  //   this auto is BoxedArray<Short_t>
  for (auto& partInfo: AuxDetID) FillWith(partInfo, -9999);
  // - pythonish C++: as the previous line, for each one in a list of containers
  //   of the same type (C++ is not python yet), using pointers to avoid copy;
  for (AuxDetMCData_t<Float_t>* cont: {
   &entryX, &entryY, &entryZ, &entryT,
   &exitX , &exitY , &exitZ, &exitT, &exitPx, &exitPy, &exitPz,
   &CombinedEnergyDep
   })
  {
    // this auto is BoxedArray<Float_t>
    for (auto& partInfo: *cont) FillWith(partInfo, -99999.);
  } // for container
  
} // microboone::AnalysisTreeDataStruct::ClearLocalData()


void microboone::AnalysisTreeDataStruct::Clear() {
  ClearLocalData();
  std::for_each
    (TrackData.begin(), TrackData.end(), std::mem_fun_ref(&TrackDataStruct::Clear));
} // microboone::AnalysisTreeDataStruct::Clear()


void microboone::AnalysisTreeDataStruct::ResizeGEANT(int nParticles) {

  // minimum size is 1, so that we always have an address
  MaxGEANTparticles = (size_t) std::max(nParticles, 1);
  
  pdg.resize(MaxGEANTparticles);
  status.resize(MaxGEANTparticles);
  Mass.resize(MaxGEANTparticles);  
  Eng.resize(MaxGEANTparticles);
  EndE.resize(MaxGEANTparticles);
  Px.resize(MaxGEANTparticles);
  Py.resize(MaxGEANTparticles);
  Pz.resize(MaxGEANTparticles);
  P.resize(MaxGEANTparticles);
  StartPointx.resize(MaxGEANTparticles);
  StartPointy.resize(MaxGEANTparticles);
  StartPointz.resize(MaxGEANTparticles);
  StartT.resize(MaxGEANTparticles); 
  EndT.resize(MaxGEANTparticles);    
  EndPointx.resize(MaxGEANTparticles);
  EndPointy.resize(MaxGEANTparticles);
  EndPointz.resize(MaxGEANTparticles);
  EndT.resize(MaxGEANTparticles);  
  theta.resize(MaxGEANTparticles);
  phi.resize(MaxGEANTparticles);
  theta_xz.resize(MaxGEANTparticles);
  theta_yz.resize(MaxGEANTparticles);
  pathlen.resize(MaxGEANTparticles);
  inTPCActive.resize(MaxGEANTparticles);
  StartPointx_tpcAV.resize(MaxGEANTparticles);
  StartPointy_tpcAV.resize(MaxGEANTparticles);
  StartPointz_tpcAV.resize(MaxGEANTparticles);
  EndPointx_tpcAV.resize(MaxGEANTparticles);
  EndPointy_tpcAV.resize(MaxGEANTparticles);
  EndPointz_tpcAV.resize(MaxGEANTparticles);    
  NumberDaughters.resize(MaxGEANTparticles);
  Mother.resize(MaxGEANTparticles);
  TrackId.resize(MaxGEANTparticles);
  process_primary.resize(MaxGEANTparticles);
  processname.resize(MaxGEANTparticles);
  MergedId.resize(MaxGEANTparticles);
  
  // auxiliary detector structure
  NAuxDets.resize(MaxGEANTparticles);
  AuxDetID.resize(MaxGEANTparticles);
  entryX.resize(MaxGEANTparticles);
  entryY.resize(MaxGEANTparticles);
  entryZ.resize(MaxGEANTparticles);
  entryT.resize(MaxGEANTparticles);
  exitX.resize(MaxGEANTparticles);
  exitY.resize(MaxGEANTparticles);
  exitZ.resize(MaxGEANTparticles);
  exitT.resize(MaxGEANTparticles);
  exitPx.resize(MaxGEANTparticles);
  exitPy.resize(MaxGEANTparticles);
  exitPz.resize(MaxGEANTparticles);
  CombinedEnergyDep.resize(MaxGEANTparticles);
  
} // microboone::AnalysisTreeDataStruct::ResizeGEANT()

void microboone::AnalysisTreeDataStruct::ResizeGenie(int nPrimaries) {
  
  // minimum size is 1, so that we always have an address
  MaxGeniePrimaries = (size_t) std::max(nPrimaries, 1);
  genie_primaries_pdg.resize(MaxGeniePrimaries);
  genie_Eng.resize(MaxGeniePrimaries);
  genie_Px.resize(MaxGeniePrimaries);
  genie_Py.resize(MaxGeniePrimaries);
  genie_Pz.resize(MaxGeniePrimaries);
  genie_P.resize(MaxGeniePrimaries);
  genie_status_code.resize(MaxGeniePrimaries);
  genie_mass.resize(MaxGeniePrimaries);
  genie_trackID.resize(MaxGeniePrimaries);
  genie_ND.resize(MaxGeniePrimaries);
  genie_mother.resize(MaxGeniePrimaries);

} // microboone::AnalysisTreeDataStruct::ResizeGenie()

void microboone::AnalysisTreeDataStruct::ResizeCry(int nPrimaries) {

  cry_primaries_pdg.resize(nPrimaries);
  cry_Eng.resize(nPrimaries);
  cry_Px.resize(nPrimaries);
  cry_Py.resize(nPrimaries);
  cry_Pz.resize(nPrimaries);
  cry_P.resize(nPrimaries);
  cry_StartPointx.resize(nPrimaries);
  cry_StartPointy.resize(nPrimaries);
  cry_StartPointz.resize(nPrimaries);  
  cry_status_code.resize(nPrimaries);
  cry_mass.resize(nPrimaries);
  cry_trackID.resize(nPrimaries);
  cry_ND.resize(nPrimaries);
  cry_mother.resize(nPrimaries);

} // microboone::AnalysisTreeDataStruct::ResizeCry()

void microboone::AnalysisTreeDataStruct::SetAddresses(
  TTree* pTree,
  const std::vector<std::string>& trackers
) {
  BranchCreator CreateBranch(pTree);

  CreateBranch("run",&run,"run/I");
  CreateBranch("subrun",&subrun,"subrun/I");
  CreateBranch("event",&event,"event/I");
  CreateBranch("evttime",&evttime,"evttime/D");
  CreateBranch("beamtime",&beamtime,"beamtime/D");
  CreateBranch("pot",&SubRunData.pot,"pot/D");
  CreateBranch("isdata",&isdata,"isdata/B");
  CreateBranch("taulife",&taulife,"taulife/D");

  if (hasHitInfo()){
    CreateBranch("no_hits",&no_hits,"no_hits/I");
    CreateBranch("hit_plane",hit_plane,"hit_plane[no_hits]/S");
    CreateBranch("hit_wire",hit_wire,"hit_wire[no_hits]/S");
    CreateBranch("hit_channel",hit_channel,"hit_channel[no_hits]/S");
    CreateBranch("hit_peakT",hit_peakT,"hit_peakT[no_hits]/F");
    CreateBranch("hit_charge",hit_charge,"hit_charge[no_hits]/F");
    CreateBranch("hit_ph",hit_ph,"hit_ph[no_hits]/F");
    CreateBranch("hit_startT",hit_startT,"hit_startT[no_hits]/F");
    CreateBranch("hit_endT",hit_endT,"hit_endT[no_hits]/F");
    CreateBranch("hit_trkid",hit_trkid,"hit_trkid[no_hits]/S");
    CreateBranch("hit_nelec",hit_nelec,"hit_nelec[no_hits]/F");
    CreateBranch("hit_energy",hit_energy,"hit_energy[no_hits]/F");
  }

  if (hasVertexInfo()){
    CreateBranch("nvtx",&nvtx,"nvtx/S");
    CreateBranch("vtx",vtx,"vtx[nvtx][3]/F");
  }  

  if (hasTrackInfo()){
    AutoResettingStringSteam sstr;
    sstr() << kMaxTrackHits;
    std::string MaxTrackHitsIndexStr("[" + sstr.str() + "]");

    kNTracker = trackers.size();
    CreateBranch("kNTracker",&kNTracker,"kNTracker/B");
    for(int i=0; i<kNTracker; i++){
      std::string TrackLabel = trackers[i];
      std::string BranchName;

      // note that if the tracker data has maximum number of tracks 0,
      // nothing is initialized (branches are not even created)
      TrackData[i].SetAddresses(pTree, TrackLabel);    
    } // for trackers
  } 

  if (hasGenieInfo()){
    CreateBranch("mcevts_truth",&mcevts_truth,"mcevts_truth/I");
    CreateBranch("nuPDG_truth",&nuPDG_truth,"nuPDG_truth/I");
    CreateBranch("ccnc_truth",&ccnc_truth,"ccnc_truth/I");
    CreateBranch("mode_truth",&mode_truth,"mode_truth/I");
    CreateBranch("enu_truth",&enu_truth,"enu_truth/F");
    CreateBranch("Q2_truth",&Q2_truth,"Q2_truth/F");
    CreateBranch("W_truth",&W_truth,"W_truth/F");
    CreateBranch("hitnuc_truth",&hitnuc_truth,"hitnuc_truth/I");
    CreateBranch("nuvtxx_truth",&nuvtxx_truth,"nuvtxx_truth/F");
    CreateBranch("nuvtxy_truth",&nuvtxy_truth,"nuvtxy_truth/F");
    CreateBranch("nuvtxz_truth",&nuvtxz_truth,"nuvtxz_truth/F");
    CreateBranch("nu_dcosx_truth",&nu_dcosx_truth,"nu_dcosx_truth/F");
    CreateBranch("nu_dcosy_truth",&nu_dcosy_truth,"nu_dcosy_truth/F");
    CreateBranch("nu_dcosz_truth",&nu_dcosz_truth,"nu_dcosz_truth/F");
    CreateBranch("lep_mom_truth",&lep_mom_truth,"lep_mom_truth/F");
    CreateBranch("lep_dcosx_truth",&lep_dcosx_truth,"lep_dcosx_truth/F");
    CreateBranch("lep_dcosy_truth",&lep_dcosy_truth,"lep_dcosy_truth/F");
    CreateBranch("lep_dcosz_truth",&lep_dcosz_truth,"lep_dcosz_truth/F");

    CreateBranch("tpx_flux",&tpx_flux,"tpx_flux/F");
    CreateBranch("tpy_flux",&tpy_flux,"tpy_flux/F");
    CreateBranch("tpz_flux",&tpz_flux,"tpz_flux/F");
    CreateBranch("tptype_flux",&tptype_flux,"tptype_flux/I");

    CreateBranch("genie_no_primaries",&genie_no_primaries,"genie_no_primaries/I");
    CreateBranch("genie_primaries_pdg",genie_primaries_pdg,"genie_primaries_pdg[genie_no_primaries]/I");
    CreateBranch("genie_Eng",genie_Eng,"genie_Eng[genie_no_primaries]/F");
    CreateBranch("genie_Px",genie_Px,"genie_Px[genie_no_primaries]/F");
    CreateBranch("genie_Py",genie_Py,"genie_Py[genie_no_primaries]/F");
    CreateBranch("genie_Pz",genie_Pz,"genie_Pz[genie_no_primaries]/F");
    CreateBranch("genie_P",genie_P,"genie_P[genie_no_primaries]/F");
    CreateBranch("genie_status_code",genie_status_code,"genie_status_code[genie_no_primaries]/I");
    CreateBranch("genie_mass",genie_mass,"genie_mass[genie_no_primaries]/F");
    CreateBranch("genie_trackID",genie_trackID,"genie_trackID[genie_no_primaries]/I");
    CreateBranch("genie_ND",genie_ND,"genie_ND[genie_no_primaries]/I");
    CreateBranch("genie_mother",genie_mother,"genie_mother[genie_no_primaries]/I");
  }

   if (hasCryInfo()){
    CreateBranch("mcevts_truthcry",&mcevts_truthcry,"mcevts_truthcry/I");  
    CreateBranch("cry_no_primaries",&cry_no_primaries,"cry_no_primaries/I");
    CreateBranch("cry_primaries_pdg",cry_primaries_pdg,"cry_primaries_pdg[cry_no_primaries]/I");
    CreateBranch("cry_Eng",cry_Eng,"cry_Eng[cry_no_primaries]/F");
    CreateBranch("cry_Px",cry_Px,"cry_Px[cry_no_primaries]/F");
    CreateBranch("cry_Py",cry_Py,"cry_Py[cry_no_primaries]/F");
    CreateBranch("cry_Pz",cry_Pz,"cry_Pz[cry_no_primaries]/F");
    CreateBranch("cry_P",cry_P,"cry_P[cry_no_primaries]/F");
    CreateBranch("cry_StartPointx",cry_StartPointx,"cry_StartPointx[cry_no_primaries]/F");
    CreateBranch("cry_StartPointy",cry_StartPointy,"cry_StartPointy[cry_no_primaries]/F");
    CreateBranch("cry_StartPointz",cry_StartPointz,"cry_StartPointz[cry_no_primaries]/F");   
    CreateBranch("cry_status_code",cry_status_code,"cry_status_code[cry_no_primaries]/I");
    CreateBranch("cry_mass",cry_mass,"cry_mass[cry_no_primaries]/F");
    CreateBranch("cry_trackID",cry_trackID,"cry_trackID[cry_no_primaries]/I");
    CreateBranch("cry_ND",cry_ND,"cry_ND[cry_no_primaries]/I");
    CreateBranch("cry_mother",cry_mother,"cry_mother[cry_no_primaries]/I");
  }  

  if (hasGeantInfo()){  
    CreateBranch("no_primaries",&no_primaries,"no_primaries/I");
    CreateBranch("geant_list_size",&geant_list_size,"geant_list_size/I");
    CreateBranch("geant_list_size_in_tpcAV",&geant_list_size_in_tpcAV,"geant_list_size_in_tpcAV/I");  
    CreateBranch("pdg",pdg,"pdg[geant_list_size]/I");
    CreateBranch("status",status,"status[geant_list_size]/I");
    CreateBranch("Mass",Mass,"Mass[geant_list_size]/F");
    CreateBranch("Eng",Eng,"Eng[geant_list_size]/F");
    CreateBranch("EndE",EndE,"EndE[geant_list_size]/F");
    CreateBranch("Px",Px,"Px[geant_list_size]/F");
    CreateBranch("Py",Py,"Py[geant_list_size]/F");
    CreateBranch("Pz",Pz,"Pz[geant_list_size]/F");
    CreateBranch("P",P,"P[geant_list_size]/F");
    CreateBranch("StartPointx",StartPointx,"StartPointx[geant_list_size]/F");
    CreateBranch("StartPointy",StartPointy,"StartPointy[geant_list_size]/F");
    CreateBranch("StartPointz",StartPointz,"StartPointz[geant_list_size]/F");
    CreateBranch("StartT",StartT,"StartT[geant_list_size]/F");
    CreateBranch("EndPointx",EndPointx,"EndPointx[geant_list_size]/F");
    CreateBranch("EndPointy",EndPointy,"EndPointy[geant_list_size]/F");
    CreateBranch("EndPointz",EndPointz,"EndPointz[geant_list_size]/F");
    CreateBranch("EndT",EndT,"EndT[geant_list_size]/F");
    CreateBranch("theta",theta,"theta[geant_list_size]/F");
    CreateBranch("phi",phi,"phi[geant_list_size]/F");
    CreateBranch("theta_xz",theta_xz,"theta_xz[geant_list_size]/F");
    CreateBranch("theta_yz",theta_yz,"theta_yz[geant_list_size]/F");
    CreateBranch("pathlen",pathlen,"pathlen[geant_list_size]/F");
    CreateBranch("inTPCActive",inTPCActive,"inTPCActive[geant_list_size]/I");  
    CreateBranch("StartPointx_tpcAV",StartPointx_tpcAV,"StartPointx_tpcAV[geant_list_size]/F");
    CreateBranch("StartPointy_tpcAV",StartPointy_tpcAV,"StartPointy_tpcAV[geant_list_size]/F");
    CreateBranch("StartPointz_tpcAV",StartPointz_tpcAV,"StartPointz_tpcAV[geant_list_size]/F");
    CreateBranch("EndPointx_tpcAV",EndPointx_tpcAV,"EndPointx_tpcAV[geant_list_size]/F");
    CreateBranch("EndPointy_tpcAV",EndPointy_tpcAV,"EndPointy_tpcAV[geant_list_size]/F");
    CreateBranch("EndPointz_tpcAV",EndPointz_tpcAV,"EndPointz_tpcAV[geant_list_size]/F");
    CreateBranch("NumberDaughters",NumberDaughters,"NumberDaughters[geant_list_size]/I");
    CreateBranch("Mother",Mother,"Mother[geant_list_size]/I");
    CreateBranch("TrackId",TrackId,"TrackId[geant_list_size]/I");
    CreateBranch("MergedId", MergedId, "MergedId[geant_list_size]/I");
    CreateBranch("process_primary",process_primary,"process_primary[geant_list_size]/I");
    CreateBranch("processname", processname);
  }

  if (hasAuxDetector()) {
    // Geant information is required to fill aux detector information.
    // if fSaveGeantInfo is not set to true, show an error message and quit!
    if (!hasGeantInfo()){
      throw art::Exception(art::errors::Configuration)
      << "Saving Auxiliary detector information requies saving GEANT information, "
      <<"please set fSaveGeantInfo flag to true in your fhicl file and rerun.\n"; 
    }    
    std::ostringstream sstr;
    sstr << "[" << kMaxAuxDets << "]";
    std::string MaxAuxDetIndexStr = sstr.str();
    CreateBranch("NAuxDets",     NAuxDets, "NAuxDets[geant_list_size]/s");
    CreateBranch("AuxDetID",     AuxDetID, "AuxDetID[geant_list_size]" + MaxAuxDetIndexStr + "/S");
    CreateBranch("AuxDetEntryX", entryX,   "AuxDetEntryX[geant_list_size]" + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetEntryY", entryY,   "AuxDetEntryY[geant_list_size]" + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetEntryZ", entryZ,   "AuxDetEntryZ[geant_list_size]" + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetEntryT", entryT,   "AuxDetEntryT[geant_list_size]" + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetExitX",  exitX,    "AuxDetExitX[geant_list_size]"  + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetExitY",  exitY,    "AuxDetExitY[geant_list_size]"  + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetExitZ",  exitZ,    "AuxDetExitZ[geant_list_size]"  + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetExitT",  exitT,    "AuxDetExitT[geant_list_size]"  + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetExitPx", exitPx,   "AuxDetExitPx[geant_list_size]" + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetExitPy", exitPy,   "AuxDetExitPy[geant_list_size]" + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetExitPz", exitPz,   "AuxDetExitPz[geant_list_size]" + MaxAuxDetIndexStr + "/F");
    CreateBranch("CombinedEnergyDep", CombinedEnergyDep,
      "CombinedEnergyDep[geant_list_size]" + MaxAuxDetIndexStr + "/F");
  } // if hasAuxDetector
  
} // microboone::AnalysisTreeDataStruct::SetAddresses()


//------------------------------------------------------------------------------
//---  AnalysisTree
//---

microboone::AnalysisTree::AnalysisTree(fhicl::ParameterSet const& pset) :
  EDAnalyzer(pset),
  fTree(nullptr), fData(nullptr),
  fDigitModuleLabel         (pset.get< std::string >("DigitModuleLabel")        ),
  fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel")         ),
  fLArG4ModuleLabel         (pset.get< std::string >("LArGeantModuleLabel")     ),
  fCalDataModuleLabel       (pset.get< std::string >("CalDataModuleLabel")      ),
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel")     ),
  fCryGenModuleLabel        (pset.get< std::string >("CryGenModuleLabel")       ), 
  fG4ModuleLabel            (pset.get< std::string >("G4ModuleLabel")           ),
  fVertexModuleLabel        (pset.get< std::string> ("VertexModuleLabel")       ),
  fTrackModuleLabel         (pset.get< std::vector<std::string> >("TrackModuleLabel")),
  fCalorimetryModuleLabel   (pset.get< std::vector<std::string> >("CalorimetryModuleLabel")),
  fParticleIDModuleLabel    (pset.get< std::vector<std::string> >("ParticleIDModuleLabel")   ),
  fPOTModuleLabel           (pset.get< std::string >("POTModuleLabel")          ),
  fUseBuffer                (pset.get< bool >("UseBuffers", false)),
  fSaveAuxDetInfo           (pset.get< bool >("SaveAuxDetInfo", false)),
  fSaveCryInfo              (pset.get< bool >("SaveCryInfo", false)),  
  fSaveGenieInfo	    (pset.get< bool >("SaveGenieInfo", false)), 
  fSaveGeantInfo	    (pset.get< bool >("SaveGeantInfo", false)), 
  fSaveHitInfo	            (pset.get< bool >("SaveHitInfo", false)), 
  fSaveTrackInfo	    (pset.get< bool >("SaveTrackInfo", false)), 
  fSaveVertexInfo	    (pset.get< bool >("SaveVertexInfo", false)),
  fCosmicTaggerAssocLabel  (pset.get<std::vector< std::string > >("CosmicTaggerAssocLabel") ),
  fFlashMatchAssocLabel (pset.get<std::vector< std::string > >("FlashMatchAssocLabel") ) 
{
  if (fSaveAuxDetInfo == true) fSaveGeantInfo = true;
  mf::LogInfo("AnalysisTree") << "Configuration:"
    << "\n  UseBuffers: " << std::boolalpha << fUseBuffer
    ;
  if (GetNTrackers() > kMaxTrackers) {
    throw art::Exception(art::errors::Configuration)
      << "AnalysisTree currently supports only up to " << kMaxTrackers
      << " tracking algorithms, but " << GetNTrackers() << " are specified."
      << "\nYou can increase kMaxTrackers and recompile.";
  } // if too many trackers
  if (fTrackModuleLabel.size() != fCalorimetryModuleLabel.size()){
    throw art::Exception(art::errors::Configuration)
      << "fTrackModuleLabel.size() = "<<fTrackModuleLabel.size()<<" does not match "
      << "fCalorimetryModuleLabel.size() = "<<fCalorimetryModuleLabel.size();
  }
  if (fTrackModuleLabel.size() != fParticleIDModuleLabel.size()){
    throw art::Exception(art::errors::Configuration)
      << "fTrackModuleLabel.size() = "<<fTrackModuleLabel.size()<<" does not match "
      << "fParticleIDModuleLabel.size() = "<<fParticleIDModuleLabel.size();
  }
} // microboone::AnalysisTree::AnalysisTree()

//-------------------------------------------------
microboone::AnalysisTree::~AnalysisTree()
{
  DestroyData();
}

void microboone::AnalysisTree::CreateTree(bool bClearData /* = false */) {
  if (!fTree) {
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("anatree","analysis tree");
  }
  CreateData(bClearData);
  SetAddresses();
} // microboone::AnalysisTree::CreateTree()


void microboone::AnalysisTree::beginSubRun(const art::SubRun& sr)
{

  art::Handle< sumdata::POTSummary > potListHandle;
  //sr.getByLabel(fPOTModuleLabel,potListHandle);

  if(sr.getByLabel(fPOTModuleLabel,potListHandle))
    SubRunData.pot=potListHandle->totpot;
  else
    SubRunData.pot=0.;

}

void microboone::AnalysisTree::analyze(const art::Event& evt)
{
  //services
  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<cheat::BackTracker> bt;
  art::ServiceHandle<util::DetectorProperties> detprop;
  art::ServiceHandle<util::LArProperties> LArProp;

  // collect the sizes which might me needed to resize the tree data structure:
  bool isMC = !evt.isRealData();
  
  // * hits
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);

  // * vertices
  art::Handle< std::vector<recob::Vertex> > vtxListHandle;
  std::vector<art::Ptr<recob::Vertex> > vtxlist;
  if (evt.getByLabel(fVertexModuleLabel,vtxListHandle))
    art::fill_ptr_vector(vtxlist, vtxListHandle);


  // * MC truth information
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
    art::fill_ptr_vector(mclist, mctruthListHandle);
    
  // *MC truth cosmic generator information
  art::Handle< std::vector<simb::MCTruth> > mctruthcryListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclistcry;
  if (fSaveCryInfo){
    if (evt.getByLabel(fCryGenModuleLabel,mctruthcryListHandle))
      art::fill_ptr_vector(mclistcry, mctruthcryListHandle);
  }       
    
  art::Ptr<simb::MCTruth> mctruthcry;
  int nCryPrimaries = 0;
   
  if (fSaveCryInfo){
    mctruthcry = mclistcry[0];      
    nCryPrimaries = mctruthcry->NParticles();  
  } 
  
  int nGeniePrimaries = 0, nGEANTparticles = 0;
  
  art::Ptr<simb::MCTruth> mctruth;
  int imc = 0;
  if (isMC) { //is MC
    // GENIE
    if (!mclist.empty()){//at least one mc record
      if (fSaveGenieInfo){
        //if (mclist[0]->NeutrinoSet()){//is neutrino
        //sometimes there can be multiple neutrino interactions in one spill
        //this is trying to find the primary interaction
        //by looking for the highest energy deposition
        std::map<art::Ptr<simb::MCTruth>,double> mctruthemap;
        for (size_t i = 0; i<hitlist.size(); i++){
          //if (hitlist[i]->View() == geo::kV){//collection view
          std::vector<sim::TrackIDE> eveIDs = bt->HitToEveID(hitlist[i]);
          for (size_t e = 0; e<eveIDs.size(); e++){
            art::Ptr<simb::MCTruth> ev_mctruth = bt->TrackIDToMCTruth(eveIDs[e].trackID);
            mctruthemap[ev_mctruth]+=eveIDs[e].energy;
          }
        //}
        }
        double maxenergy = -1;
        int imc0 = 0;
        for (std::map<art::Ptr<simb::MCTruth>,double>::iterator ii=mctruthemap.begin(); ii!=mctruthemap.end(); ++ii){
          if ((ii->second)>maxenergy){
            maxenergy = ii->second;
            mctruth = ii->first;
            imc = imc0;
          }
          imc0++;
        }

        imc = 0; //set imc to 0 to solve a confusion for BNB+cosmic files where there are two MCTruth
        mctruth = mclist[0];

        if (mctruth->NeutrinoSet()) nGeniePrimaries = mctruth->NParticles();
      } //end (fSaveGenieInfo)
      
      const sim::ParticleList& plist = bt->ParticleList();
      nGEANTparticles = plist.size();

      // to know the number of particles in AV would require
      // looking at all of them; so we waste some memory here
    } // if have MC truth
    LOG_DEBUG("AnalysisTree") << "Expected "
      << nGEANTparticles << " GEANT particles, "
      << nGeniePrimaries << " GENIE particles";
  } // if MC
  
  CreateData(); // tracker data is created with default constructor
  if (fSaveGenieInfo)
    fData->ResizeGenie(nGeniePrimaries);
  if (fSaveCryInfo)
    fData->ResizeCry(nCryPrimaries);
  if (fSaveGeantInfo)    
    fData->ResizeGEANT(nGEANTparticles);
  fData->ClearLocalData(); // don't bother clearing tracker data yet
  
//  const size_t Nplanes       = 3; // number of wire planes; pretty much constant...
  const size_t NTrackers = GetNTrackers(); // number of trackers passed into fTrackModuleLabel
  const size_t NHits     = hitlist.size(); // number of hits
  const size_t NVertices = vtxlist.size(); // number of vertices
  // make sure there is the data, the tree and everything;
  CreateTree();

  /// transfer the run and subrun data to the tree data object
//  fData->RunData = RunData;
  fData->SubRunData = SubRunData;

  fData->isdata = int(!isMC);
  
  std::vector< art::Handle< std::vector<recob::Track> > > trackListHandle(NTrackers);
  std::vector< std::vector<art::Ptr<recob::Track> > > tracklist(NTrackers);
  for (unsigned int it = 0; it < NTrackers; ++it){
    if (evt.getByLabel(fTrackModuleLabel[it],trackListHandle[it]))
      art::fill_ptr_vector(tracklist[it], trackListHandle[it]);
  }

  art::Handle< std::vector<simb::MCFlux> > mcfluxListHandle;
  std::vector<art::Ptr<simb::MCFlux> > fluxlist;
  if (evt.getByLabel(fGenieGenModuleLabel,mcfluxListHandle))
    art::fill_ptr_vector(fluxlist, mcfluxListHandle);

  std::vector<const sim::AuxDetSimChannel*> fAuxDetSimChannels;
  if (fSaveAuxDetInfo){
    evt.getView(fLArG4ModuleLabel, fAuxDetSimChannels);
  }

  std::vector<const sim::SimChannel*> fSimChannels;
  evt.getView(fLArG4ModuleLabel, fSimChannels);

  fData->run = evt.run();
  fData->subrun = evt.subRun();
  fData->event = evt.id().event();

  art::Timestamp ts = evt.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  fData->evttime = tts.AsDouble();

  //copied from MergeDataPaddles.cxx
  art::Handle< raw::BeamInfo > beam;
  if (evt.getByLabel("beam",beam)){
    fData->beamtime = (double)beam->get_t_ms();
    fData->beamtime/=1000.; //in second
  }

//  std::cout<<detprop->NumberTimeSamples()<<" "<<detprop->ReadOutWindowSize()<<std::endl;
//  std::cout<<geom->DetHalfHeight()*2<<" "<<geom->DetHalfWidth()*2<<" "<<geom->DetLength()<<std::endl;
//  std::cout<<geom->Nwires(0)<<" "<<geom->Nwires(1)<<" "<<geom->Nwires(2)<<std::endl;

  //hit information
  if (fSaveHitInfo){
    fData->no_hits = (int) NHits;
    if (NHits > kMaxHits) {
      // got this error? consider increasing kMaxHits
      // (or ask for a redesign using vectors)
      mf::LogError("AnalysisTree:limits") << "event has " << NHits
        << " hits, only kMaxHits=" << kMaxHits << " stored in tree";
    }
    for (size_t i = 0; i < NHits && i < kMaxHits ; ++i){//loop over hits
      fData->hit_channel[i] = hitlist[i]->Channel();
      fData->hit_plane[i]   = hitlist[i]->WireID().Plane;
      fData->hit_wire[i]    = hitlist[i]->WireID().Wire;
      fData->hit_peakT[i]   = hitlist[i]->PeakTime();
      fData->hit_charge[i]  = hitlist[i]->Integral();
      fData->hit_ph[i]  = hitlist[i]->PeakAmplitude();
      fData->hit_startT[i] = hitlist[i]->PeakTimeMinusRMS();
      fData->hit_endT[i] = hitlist[i]->PeakTimePlusRMS();
      /*
      for (unsigned int it=0; it<fTrackModuleLabel.size();++it){
        art::FindManyP<recob::Track> fmtk(hitListHandle,evt,fTrackModuleLabel[it]);
        if (fmtk.at(i).size()!=0){
          hit_trkid[it][i] = fmtk.at(i)[0]->ID();
        }
        else
          hit_trkid[it][i] = 0;
      }
      */

      if (!evt.isRealData()){
         fData -> hit_nelec[i] = 0;
         fData -> hit_energy[i] = 0;
         const sim::SimChannel* chan = 0;
         for(size_t sc = 0; sc < fSimChannels.size(); ++sc){
           if(fSimChannels[sc]->Channel() == hitlist[i]->Channel()) chan = fSimChannels[sc];
         }
         if (chan){
           const std::map<unsigned short, std::vector<sim::IDE> >& tdcidemap = chan->TDCIDEMap();
           for(auto mapitr = tdcidemap.begin(); mapitr != tdcidemap.end(); mapitr++){
             // loop over the vector of IDE objects.
             const std::vector<sim::IDE> idevec = (*mapitr).second;
             for(size_t iv = 0; iv < idevec.size(); ++iv){
                fData -> hit_nelec[i] += idevec[iv].numElectrons;
                fData -> hit_energy[i] += idevec[iv].energy;
             }
           }
         }
       }
    }

    if (evt.getByLabel(fHitsModuleLabel,hitListHandle)){
      //Find tracks associated with hits
      art::FindManyP<recob::Track> fmtk(hitListHandle,evt,fTrackModuleLabel[0]);
      for (size_t i = 0; i < NHits && i < kMaxHits ; ++i){//loop over hits
        if (fmtk.isValid()){
	  if (fmtk.at(i).size()!=0){
	    fData->hit_trkid[i] = fmtk.at(i)[0]->ID();
	  }
	  else
	    fData->hit_trkid[i] = -1;
        }
      }
    }
  }// end (fSaveHitInfo) 

  //vertex information
  if (fSaveVertexInfo){
    fData->nvtx = NVertices;
    if (NVertices > kMaxVertices){
      // got this error? consider increasing kMaxVerticestra
      // (or ask for a redesign using vectors)
      mf::LogError("AnalysisTree:limits") << "event has " << NVertices
        << " vertices, only kMaxVertices=" << kMaxVertices << " stored in tree";
    }
    for (size_t i = 0; i < NVertices && i < kMaxVertices ; ++i){//loop over hits
      Double_t xyz[3] = {};
      vtxlist[i]->XYZ(xyz);
      for (size_t j = 0; j<3; ++j) fData->vtx[i][j] = xyz[j];
    }
  }// end (fSaveVertexInfo)
    
  //track information for multiple trackers
  if (fSaveTrackInfo){
    for (unsigned int iTracker=0; iTracker < NTrackers; ++iTracker){
      AnalysisTreeDataStruct::TrackDataStruct& TrackerData = fData->GetTrackerData(iTracker);
    
      size_t NTracks = tracklist[iTracker].size();
      // allocate enough space for this number of tracks (but at least for one of them!)
      TrackerData.SetMaxTracks(std::max(NTracks, (size_t) 1));
      TrackerData.Clear(); // clear all the data
    
      TrackerData.ntracks = (int) NTracks;
    
      // now set the tree addresses to the newly allocated memory;
      // this creates the tree branches in case they are not there yet
      SetTrackerAddresses(iTracker);
      if (NTracks > TrackerData.GetMaxTracks()) {
        // got this error? it might be a bug,
        // since we are supposed to have allocated enough space to fit all tracks
        mf::LogError("AnalysisTree:limits") << "event has " << NTracks
          << " " << fTrackModuleLabel[iTracker] << " tracks, only "
          << TrackerData.GetMaxTracks() << " stored in tree";
      }
    
      //call the track momentum algorithm that gives you momentum based on track range
      trkf::TrackMomentumCalculator trkm;
       
      for(size_t iTrk=0; iTrk < NTracks; ++iTrk){//loop over tracks
      
        //Cosmic Tagger information
        art::FindManyP<anab::CosmicTag> fmct(trackListHandle[iTracker],evt,fCosmicTaggerAssocLabel[iTracker]);
        if (fmct.isValid()){          
          TrackerData.trkncosmictags_tagger[iTrk]     = fmct.at(iTrk).size();
          if (fmct.at(iTrk).size()>0){
            if(fmct.at(iTrk).size()>1)
              std::cerr << "\n Warning : more than one cosmic tag per track in module! assigning the first tag to the track" << fCosmicTaggerAssocLabel[iTracker];
            TrackerData.trkcosmicscore_tagger[iTrk] = fmct.at(iTrk).at(0)->CosmicScore();
            TrackerData.trkcosmictype_tagger[iTrk] = fmct.at(iTrk).at(0)->CosmicType();
          }
        }

        //Flash match compatibility information
        //Unlike CosmicTagger, Flash match doesn't assign a cosmic tag for every track. For those tracks, AnalysisTree initializes them with -9999 or -99999
        art::FindManyP<anab::CosmicTag> fmbfm(trackListHandle[iTracker],evt,fFlashMatchAssocLabel[iTracker]);
        if (fmbfm.isValid()){  
          TrackerData.trkncosmictags_flashmatch[iTrk] = fmbfm.at(iTrk).size();
          if (fmbfm.at(iTrk).size()>0){
            if(fmbfm.at(iTrk).size()>1) 
              std::cerr << "\n Warning : more than one cosmic tag per track in module! assigning the first tag to the track" << fFlashMatchAssocLabel[iTracker];
  	      TrackerData.trkcosmicscore_flashmatch[iTrk] = fmbfm.at(iTrk).at(0)->CosmicScore();
              TrackerData.trkcosmictype_flashmatch[iTrk] = fmbfm.at(iTrk).at(0)->CosmicType();
	    //std::cout<<"\n"<<evt.event()<<"\t"<<iTrk<<"\t"<<fmbfm.at(iTrk).at(0)->CosmicScore()<<"\t"<<fmbfm.at(iTrk).at(0)->CosmicType();
          }
        }
     			 	   
        art::Ptr<recob::Track> ptrack(trackListHandle[iTracker], iTrk);
        const recob::Track& track = *ptrack;
      
        TVector3 pos, dir_start, dir_end, end;        

        double tlen = 0., mom = 0.;
        int TrackID = -1;
      
        int ntraj = 0;
        //we need to use Bezier methods for Bezier tracks
        if (fTrackModuleLabel[iTracker].find("beziertracker")!=std::string::npos) {
          trkf::BezierTrack btrack(*ptrack);
          ntraj = btrack.NSegments();
          if(ntraj > 0) {
            double xyz[3];
            btrack.GetTrackPoint(0,xyz);
            pos.SetXYZ(xyz[0],xyz[1],xyz[2]);
            btrack.GetTrackDirection(0,xyz);
            dir_start.SetXYZ(xyz[0],xyz[1],xyz[2]);
            btrack.GetTrackDirection(1,xyz);
            dir_end.SetXYZ(xyz[0],xyz[1],xyz[2]);
            btrack.GetTrackPoint(1,xyz);
            end.SetXYZ(xyz[0],xyz[1],xyz[2]);

            tlen = btrack.GetLength();
            if (btrack.NumberFitMomentum() > 0)
              mom = btrack.VertexMomentum();
            // fill bezier track reco branches
            TrackID = iTrk;  //bezier has some screwed up track IDs
          }
        }
        else {   //use the normal methods for other kinds of tracks
          ntraj = track.NumberTrajectoryPoints();
          if (ntraj > 0) {
            pos       = track.Vertex();
            dir_start = track.VertexDirection();
            dir_end   = track.EndDirection();
            end       = track.End();

            tlen        = length(track);
            if(track.NumberFitMomentum() > 0)
              mom = track.VertexMomentum();
            // fill non-bezier-track reco branches
            TrackID = track.ID();
          }
        }
      
        if (ntraj > 0) {
          double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
          double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());
          double dpos = bdist(pos);
          double dend = bdist(end);
        
          TrackerData.trkId[iTrk]                 = TrackID;
          TrackerData.trkstartx[iTrk]             = pos.X();
          TrackerData.trkstarty[iTrk]             = pos.Y();
          TrackerData.trkstartz[iTrk]             = pos.Z();
          TrackerData.trkstartd[iTrk]		  = dpos;
          TrackerData.trkendx[iTrk]		  = end.X();
          TrackerData.trkendy[iTrk]		  = end.Y();
          TrackerData.trkendz[iTrk]		  = end.Z();
          TrackerData.trkendd[iTrk]		  = dend;
          TrackerData.trktheta[iTrk]		  = dir_start.Theta();
          TrackerData.trkphi[iTrk]		  = dir_start.Phi();
          TrackerData.trkstartdcosx[iTrk]	  = dir_start.X();
          TrackerData.trkstartdcosy[iTrk]	  = dir_start.Y();
          TrackerData.trkstartdcosz[iTrk]	  = dir_start.Z();
          TrackerData.trkenddcosx[iTrk] 	  = dir_end.X();
          TrackerData.trkenddcosy[iTrk] 	  = dir_end.Y();
          TrackerData.trkenddcosz[iTrk] 	  = dir_end.Z();
          TrackerData.trkthetaxz[iTrk]  	  = theta_xz;
          TrackerData.trkthetayz[iTrk]  	  = theta_yz;
          TrackerData.trkmom[iTrk]		  = mom;
          TrackerData.trklen[iTrk]		  = tlen;
          TrackerData.trkmomrange[iTrk] 	  = trkm.GetTrackMomentum(tlen,13);
          TrackerData.trkmommschi2[iTrk]	  = trkm.GetMomentumMultiScatterChi2(ptrack);
          TrackerData.trkmommsllhd[iTrk]	  = trkm.GetMomentumMultiScatterLLHD(ptrack);
        } // if we have trajectory

      // find vertices associated with this track
      /*
      art::FindMany<recob::Vertex> fmvtx(trackListHandle[iTracker], evt, fVertexModuleLabel[iTracker]);
      if(fmvtx.isValid()) {
        std::vector<const recob::Vertex*> verts = fmvtx.at(iTrk);
        // should have two at most
        for(size_t ivx = 0; ivx < verts.size(); ++ivx) {
          verts[ivx]->XYZ(xyz);
          // find the vertex in TrackerData to get the index
          short theVtx = -1;
          for(short jvx = 0; jvx < TrackerData.nvtx; ++jvx) {
            if(TrackerData.vtx[jvx][2] == xyz[2]) {
              theVtx = jvx;
              break;
            }
          } // jvx
          // decide if it should be assigned to the track Start or End.
          // A simple dz test should suffice
          if(fabs(xyz[2] - TrackerData.trkstartz[iTrk]) < 
             fabs(xyz[2] - TrackerData.trkendz[iTrk])) {
            TrackerData.trksvtxid[iTrk] = theVtx;
          } else {
            TrackerData.trkevtxid[iTrk] = theVtx;
          }
        } // vertices
      } // fmvtx.isValid()
      */
      Float_t minsdist = 10000;
      Float_t minedist = 10000;
      for (int ivx = 0; ivx < fData->nvtx && ivx < kMaxVertices; ++ivx){
        Float_t sdist = sqrt(pow(TrackerData.trkstartx[iTrk]-fData->vtx[ivx][0],2)+
                             pow(TrackerData.trkstarty[iTrk]-fData->vtx[ivx][1],2)+
                             pow(TrackerData.trkstartz[iTrk]-fData->vtx[ivx][2],2));
        Float_t edist = sqrt(pow(TrackerData.trkendx[iTrk]-fData->vtx[ivx][0],2)+
                             pow(TrackerData.trkendy[iTrk]-fData->vtx[ivx][1],2)+
                             pow(TrackerData.trkendz[iTrk]-fData->vtx[ivx][2],2));
        if (sdist<minsdist){
          minsdist = sdist;
          if (minsdist<10) TrackerData.trksvtxid[iTrk] = ivx;
        }
        if (edist<minedist){
          minedist = edist;
          if (minedist<10) TrackerData.trkevtxid[iTrk] = ivx;
        }
      }
      // find particle ID info
      art::FindMany<anab::ParticleID> fmpid(trackListHandle[iTracker], evt, fParticleIDModuleLabel[iTracker]);
      if(fmpid.isValid()) {
        std::vector<const anab::ParticleID*> pids = fmpid.at(iTrk);
        //if(pids.size() > 1) {
          //mf::LogError("AnalysisTree:limits")
            //<< "the " << fTrackModuleLabel[iTracker] << " track #" << iTrk
            //<< " has " << pids.size() 
            //<< " set of ParticleID variables. Only one stored in the tree";
        //}
        for (size_t ipid = 0; ipid < pids.size(); ++ipid){
	  if (!pids[ipid]->PlaneID().isValid) continue;
	  int planenum = pids[ipid]->PlaneID().Plane;
	  if (planenum<0||planenum>2) continue;
          TrackerData.trkpidpdg[iTrk][planenum] = pids[ipid]->Pdg();
          TrackerData.trkpidchi[iTrk][planenum] = pids[ipid]->MinChi2();
          TrackerData.trkpidchipr[iTrk][planenum] = pids[ipid]->Chi2Proton();
          TrackerData.trkpidchika[iTrk][planenum] = pids[ipid]->Chi2Kaon();
          TrackerData.trkpidchipi[iTrk][planenum] = pids[ipid]->Chi2Pion();
          TrackerData.trkpidchimu[iTrk][planenum] = pids[ipid]->Chi2Muon();
          TrackerData.trkpidpida[iTrk][planenum] = pids[ipid]->PIDA();
        }
      } // fmpid.isValid()
      
      

      art::FindMany<anab::Calorimetry> fmcal(trackListHandle[iTracker], evt, fCalorimetryModuleLabel[iTracker]);
      if (fmcal.isValid()){
        std::vector<const anab::Calorimetry*> calos = fmcal.at(iTrk);
        if (calos.size() > TrackerData.GetMaxPlanesPerTrack(iTrk)) {
          // if you get this message, there is probably a bug somewhere since
          // the calorimetry planes should be 3.
          mf::LogError("AnalysisTree:limits")
            << "the " << fTrackModuleLabel[iTracker] << " track #" << iTrk
            << " has " << calos.size() << " planes for calorimetry , only "
            << TrackerData.GetMaxPlanesPerTrack(iTrk) << " stored in tree";
        }
        for (size_t ical = 0; ical<calos.size(); ++ical){
	  if (!calos[ical]) continue;
	  if (!calos[ical]->PlaneID().isValid) continue;
	  int planenum = calos[ical]->PlaneID().Plane;
	  if (planenum<0||planenum>2) continue;
          TrackerData.trkke[iTrk][planenum]    = calos[ical]->KineticEnergy();
          TrackerData.trkrange[iTrk][planenum] = calos[ical]->Range();
          //For now make the second argument as 13 for muons. 
          TrackerData.trkpitchc[iTrk][planenum]= calos[ical] -> TrkPitchC();
          const size_t NHits = calos[ical] -> dEdx().size();
          TrackerData.ntrkhits[iTrk][planenum] = (int) NHits;
          if (NHits > TrackerData.GetMaxHitsPerTrack(iTrk, planenum)) {
            // if you get this error, you'll have to increase kMaxTrackHits
            mf::LogError("AnalysisTree:limits")
              << "the " << fTrackModuleLabel[iTracker] << " track #" << iTrk
              << " has " << NHits << " hits on calorimetry plane #" << planenum
              <<", only "
              << TrackerData.GetMaxHitsPerTrack(iTrk, planenum) << " stored in tree";
          }
          for(size_t iTrkHit = 0; iTrkHit < NHits && iTrkHit < TrackerData.GetMaxHitsPerTrack(iTrk, planenum); ++iTrkHit) {
            TrackerData.trkdedx[iTrk][planenum][iTrkHit]  = (calos[ical] -> dEdx())[iTrkHit];
            TrackerData.trkdqdx[iTrk][planenum][iTrkHit]  = (calos[ical] -> dQdx())[iTrkHit];
            TrackerData.trkresrg[iTrk][planenum][iTrkHit] = (calos[ical] -> ResidualRange())[iTrkHit];
            const auto& TrkPos = (calos[ical] -> XYZ())[iTrkHit];
            auto& TrkXYZ = TrackerData.trkxyz[iTrk][planenum][iTrkHit];
            TrkXYZ[0] = TrkPos.X();
            TrkXYZ[1] = TrkPos.Y();
            TrkXYZ[2] = TrkPos.Z();
          } // for track hits
        } // for calorimetry info
        if(TrackerData.ntrkhits[iTrk][0] > TrackerData.ntrkhits[iTrk][1] && TrackerData.ntrkhits[iTrk][0] > TrackerData.ntrkhits[iTrk][2]) TrackerData.trkpidbestplane[iTrk] = 0;
        else if(TrackerData.ntrkhits[iTrk][1] > TrackerData.ntrkhits[iTrk][0] && TrackerData.ntrkhits[iTrk][1] > TrackerData.ntrkhits[iTrk][2]) TrackerData.trkpidbestplane[iTrk] = 1;
        else if(TrackerData.ntrkhits[iTrk][2] > TrackerData.ntrkhits[iTrk][0] && TrackerData.ntrkhits[iTrk][2] > TrackerData.ntrkhits[iTrk][1]) TrackerData.trkpidbestplane[iTrk] = 2;
        else if(TrackerData.ntrkhits[iTrk][2] == TrackerData.ntrkhits[iTrk][0] && TrackerData.ntrkhits[iTrk][2] > TrackerData.ntrkhits[iTrk][1]) TrackerData.trkpidbestplane[iTrk] = 2;
        else if(TrackerData.ntrkhits[iTrk][2] == TrackerData.ntrkhits[iTrk][1] && TrackerData.ntrkhits[iTrk][2] > TrackerData.ntrkhits[iTrk][0]) TrackerData.trkpidbestplane[iTrk] = 2;
        else if(TrackerData.ntrkhits[iTrk][1] == TrackerData.ntrkhits[iTrk][0] && TrackerData.ntrkhits[iTrk][1] > TrackerData.ntrkhits[iTrk][2]) TrackerData.trkpidbestplane[iTrk] = 0;
        else if(TrackerData.ntrkhits[iTrk][1] == TrackerData.ntrkhits[iTrk][0] && TrackerData.ntrkhits[iTrk][1] == TrackerData.ntrkhits[iTrk][2]) TrackerData.trkpidbestplane[iTrk] = 2;
      } // if has calorimetry info

      //track truth information
      if (isMC){
        //get the hits on each plane
        art::FindManyP<recob::Hit>      fmht(trackListHandle[iTracker], evt, fTrackModuleLabel[iTracker]);
        std::vector< art::Ptr<recob::Hit> > allHits = fmht.at(iTrk);
        std::vector< art::Ptr<recob::Hit> > hits[kNplanes];

        for(size_t ah = 0; ah < allHits.size(); ++ah){
          if (/* allHits[ah]->WireID().Plane >= 0 && */ // always true
            allHits[ah]->WireID().Plane <  3){
            hits[allHits[ah]->WireID().Plane].push_back(allHits[ah]);
          }
        }
        
        for (size_t ipl = 0; ipl < 3; ++ipl){
          double maxe = 0;
          HitsPurity(hits[ipl],TrackerData.trkidtruth[iTrk][ipl],TrackerData.trkpurtruth[iTrk][ipl],maxe);
        //std::cout<<"\n"<<iTracker<<"\t"<<iTrk<<"\t"<<ipl<<"\t"<<trkidtruth[iTracker][iTrk][ipl]<<"\t"<<trkpurtruth[iTracker][iTrk][ipl]<<"\t"<<maxe;
          if (TrackerData.trkidtruth[iTrk][ipl]>0){
            const art::Ptr<simb::MCTruth> mc = bt->TrackIDToMCTruth(TrackerData.trkidtruth[iTrk][ipl]);
            TrackerData.trkorigin[iTrk][ipl] = mc->Origin();
            const simb::MCParticle *particle = bt->TrackIDToParticle(TrackerData.trkidtruth[iTrk][ipl]);
            double tote = 0;
            std::vector<sim::IDE> vide(bt->TrackIDToSimIDE(TrackerData.trkidtruth[iTrk][ipl]));
            for (const sim::IDE& ide: vide) {
               tote += ide.energy;
               TrackerData.trksimIDEenergytruth[iTrk][ipl] = ide.energy;
               TrackerData.trksimIDExtruth[iTrk][ipl] = ide.x;
               TrackerData.trksimIDEytruth[iTrk][ipl] = ide.y;
               TrackerData.trksimIDEztruth[iTrk][ipl] = ide.z;
            }
            TrackerData.trkpdgtruth[iTrk][ipl] = particle->PdgCode();
            TrackerData.trkefftruth[iTrk][ipl] = maxe/(tote/kNplanes); //tote include both induction and collection energies
          //std::cout<<"\n"<<trkpdgtruth[iTracker][iTrk][ipl]<<"\t"<<trkefftruth[iTracker][iTrk][ipl];
          }
        }
      }//end if (isMC)
    }//end loop over track
  }//end loop over track module labels
 }// end (fSaveTrackInfo) 
  
  /*trkf::TrackMomentumCalculator trkm;  
  std::cout<<"\t"<<trkm.GetTrackMomentum(200,2212)<<"\t"<<trkm.GetTrackMomentum(-10, 13)<<"\t"<<trkm.GetTrackMomentum(300,-19)<<"\n";
*/
  //mc truth information
  if (isMC){
    if (fSaveCryInfo){ 
      //store cry (cosmic generator information) 
      fData->mcevts_truthcry = mclistcry.size();
      fData->cry_no_primaries = nCryPrimaries;
      //fData->cry_no_primaries;
      for(Int_t iPartc = 0; iPartc < mctruthcry->NParticles(); ++iPartc){
        const simb::MCParticle& partc(mctruthcry->GetParticle(iPartc));
        fData->cry_primaries_pdg[iPartc]=partc.PdgCode();
        fData->cry_Eng[iPartc]=partc.E();
        fData->cry_Px[iPartc]=partc.Px();
        fData->cry_Py[iPartc]=partc.Py();
        fData->cry_Pz[iPartc]=partc.Pz();
        fData->cry_P[iPartc]=partc.P();
	fData->cry_StartPointx[iPartc] = partc.Vx();
	fData->cry_StartPointy[iPartc] = partc.Vy();
	fData->cry_StartPointz[iPartc] = partc.Vz();	
        fData->cry_status_code[iPartc]=partc.StatusCode();
        fData->cry_mass[iPartc]=partc.Mass();
        fData->cry_trackID[iPartc]=partc.TrackId();
        fData->cry_ND[iPartc]=partc.NumberDaughters();
        fData->cry_mother[iPartc]=partc.Mother();
      } // for cry particles  
    }// end fSaveCryInfo   
    //save neutrino interaction information
    fData->mcevts_truth = mclist.size();
    if (fData->mcevts_truth > 0){//at least one mc record
    if (fSaveGenieInfo){
      if (mctruth->NeutrinoSet()){
        fData->nuPDG_truth = mctruth->GetNeutrino().Nu().PdgCode();
        fData->ccnc_truth = mctruth->GetNeutrino().CCNC();
        fData->mode_truth = mctruth->GetNeutrino().Mode();
        fData->Q2_truth = mctruth->GetNeutrino().QSqr();
        fData->W_truth = mctruth->GetNeutrino().W();
        fData->hitnuc_truth = mctruth->GetNeutrino().HitNuc();
        fData->enu_truth = mctruth->GetNeutrino().Nu().E();
        fData->nuvtxx_truth = mctruth->GetNeutrino().Nu().Vx();
        fData->nuvtxy_truth = mctruth->GetNeutrino().Nu().Vy();
        fData->nuvtxz_truth = mctruth->GetNeutrino().Nu().Vz();
        if (mctruth->GetNeutrino().Nu().P()){
          fData->nu_dcosx_truth = mctruth->GetNeutrino().Nu().Px()/mctruth->GetNeutrino().Nu().P();
          fData->nu_dcosy_truth = mctruth->GetNeutrino().Nu().Py()/mctruth->GetNeutrino().Nu().P();
          fData->nu_dcosz_truth = mctruth->GetNeutrino().Nu().Pz()/mctruth->GetNeutrino().Nu().P();
        }
        fData->lep_mom_truth = mctruth->GetNeutrino().Lepton().P();
        if (mctruth->GetNeutrino().Lepton().P()){
          fData->lep_dcosx_truth = mctruth->GetNeutrino().Lepton().Px()/mctruth->GetNeutrino().Lepton().P();
          fData->lep_dcosy_truth = mctruth->GetNeutrino().Lepton().Py()/mctruth->GetNeutrino().Lepton().P();
          fData->lep_dcosz_truth = mctruth->GetNeutrino().Lepton().Pz()/mctruth->GetNeutrino().Lepton().P();
        }
        //flux information
        art::Ptr<simb::MCFlux>  mcflux = fluxlist[imc];
        fData->tpx_flux = mcflux->ftpx;
        fData->tpy_flux = mcflux->ftpy;
        fData->tpz_flux = mcflux->ftpz;
        fData->tptype_flux = mcflux->ftptype;

        //genie particles information
        fData->genie_no_primaries = mctruth->NParticles();

        size_t StoreParticles = std::min((size_t) fData->genie_no_primaries, fData->GetMaxGeniePrimaries());
        if (fData->genie_no_primaries > (int) StoreParticles) {
          // got this error? it might be a bug,
          // since the structure should have enough room for everything
          mf::LogError("AnalysisTree:limits") << "event has "
            << fData->genie_no_primaries << " MC particles, only "
            << StoreParticles << " stored in tree";
        }
        for(size_t iPart = 0; iPart < StoreParticles; ++iPart){
          const simb::MCParticle& part(mctruth->GetParticle(iPart));
          fData->genie_primaries_pdg[iPart]=part.PdgCode();
          fData->genie_Eng[iPart]=part.E();
          fData->genie_Px[iPart]=part.Px();
          fData->genie_Py[iPart]=part.Py();
          fData->genie_Pz[iPart]=part.Pz();
          fData->genie_P[iPart]=part.P();
          fData->genie_status_code[iPart]=part.StatusCode();
          fData->genie_mass[iPart]=part.Mass();
          fData->genie_trackID[iPart]=part.TrackId();
          fData->genie_ND[iPart]=part.NumberDaughters();
          fData->genie_mother[iPart]=part.Mother();
        } // for particle
      } //if neutrino set
    }// end (fSaveGenieInfo)  

      //GEANT particles information
      if (fSaveGeantInfo){ 
        const sim::ParticleList& plist = bt->ParticleList();
        
        std::string pri("primary");
        int primary=0;
        int active = 0;
        int geant_particle=0;
        sim::ParticleList::const_iterator itPart = plist.begin(),
          pend = plist.end(); // iterator to pairs (track id, particle)
        	  
        for(size_t iPart = 0; (iPart < plist.size()) && (itPart != pend); ++iPart)
        {
          const simb::MCParticle* pPart = (itPart++)->second;
          if (!pPart) {
            throw art::Exception(art::errors::LogicError)
              << "GEANT particle #" << iPart << " returned a null pointer";
          }
          
          ++geant_particle;
          bool isPrimary = pPart->Process() == pri;
          if (isPrimary) ++primary;
          
          int TrackID = pPart->TrackId();

          TVector3 mcstart, mcend;
          double plen = length(*pPart, mcstart, mcend);

          bool isActive = plen != 0;
          if (plen) active++;

          if (iPart < fData->GetMaxGEANTparticles()) {
	   if (pPart->E()>0.01){
            fData->process_primary[iPart] = int(isPrimary);
            fData->processname[iPart]= pPart->Process();
            fData->Mother[iPart]=pPart->Mother();
            fData->TrackId[iPart]=TrackID;
            fData->pdg[iPart]=pPart->PdgCode();
            fData->status[iPart] = pPart->StatusCode();
            fData->Eng[iPart]=pPart->E();
	    fData->EndE[iPart]=pPart->EndE();
            fData->Mass[iPart]=pPart->Mass();
            fData->Px[iPart]=pPart->Px();
            fData->Py[iPart]=pPart->Py();
            fData->Pz[iPart]=pPart->Pz();
            fData->P[iPart]=pPart->Momentum().Vect().Mag();
            fData->StartPointx[iPart]=pPart->Vx();
            fData->StartPointy[iPart]=pPart->Vy();
            fData->StartPointz[iPart]=pPart->Vz();
            fData->StartT[iPart] = pPart->T();
            fData->EndPointx[iPart]=pPart->EndPosition()[0];
            fData->EndPointy[iPart]=pPart->EndPosition()[1];
            fData->EndPointz[iPart]=pPart->EndPosition()[2];
            fData->EndT[iPart] = pPart->EndT();
            fData->theta[iPart] = pPart->Momentum().Theta();
            fData->phi[iPart] = pPart->Momentum().Phi();
            fData->theta_xz[iPart] = std::atan2(pPart->Px(), pPart->Pz());
            fData->theta_yz[iPart] = std::atan2(pPart->Py(), pPart->Pz());
            fData->pathlen[iPart]  = plen;
            fData->NumberDaughters[iPart]=pPart->NumberDaughters();
            fData->inTPCActive[iPart] = int(isActive);
            if (isActive){	  
              fData->StartPointx_tpcAV[iPart] = mcstart.X();
              fData->StartPointy_tpcAV[iPart] = mcstart.Y();
              fData->StartPointz_tpcAV[iPart] = mcstart.Z();
              fData->EndPointx_tpcAV[iPart] = mcend.X();
              fData->EndPointy_tpcAV[iPart] = mcend.Y();
              fData->EndPointz_tpcAV[iPart] = mcend.Z();
            }		       
           } 
            //access auxiliary detector parameters
            if (fSaveAuxDetInfo) {
              unsigned short nAD = 0; // number of cells that particle hit
              
              // find deposit of this particle in each of the detector cells
              for (const sim::AuxDetSimChannel* c: fAuxDetSimChannels) {
        	
        	// find if this cell has a contribution (IDE) from this particle,
        	// and which one
        	const std::vector<sim::AuxDetIDE>& setOfIDEs = c->AuxDetIDEs();
        	// using a C++ "lambda" function here; this one:
        	// - sees only TrackID from the current scope
        	// - takes one parameter: the AuxDetIDE to be tested
        	// - returns if that IDE belongs to the track we are looking for
        	std::vector<sim::AuxDetIDE>::const_iterator iIDE
        	  = std::find_if(
        	    setOfIDEs.begin(), setOfIDEs.end(),
        	    [TrackID](const sim::AuxDetIDE& IDE){ return IDE.trackID == TrackID; }
        	  );
        	if (iIDE == setOfIDEs.end()) continue;
        	
        	// now iIDE points to the energy released by the track #i (TrackID)
        	
              // look for IDE with matching trackID
              // find trackIDs stored in setOfIDEs with the same trackID, but negative,
              // this is an untracked particle who's energy should be added as deposited by this original trackID
              float totalE = 0.; // total energy deposited around by the GEANT particle in this cell
              for(const auto& adtracks: setOfIDEs) {
                 if( fabs(adtracks.trackID) == TrackID )
                   totalE += adtracks.energyDeposited;
              } // for
              
              // fill the structure
              if (nAD < kMaxAuxDets) {
                fData->AuxDetID[iPart][nAD] = c->AuxDetID();
                fData->entryX[iPart][nAD]   = iIDE->entryX;
                fData->entryY[iPart][nAD]   = iIDE->entryY;
                fData->entryZ[iPart][nAD]   = iIDE->entryZ;
                fData->entryT[iPart][nAD]   = iIDE->entryT;
                fData->exitX[iPart][nAD]    = iIDE->exitX;
                fData->exitY[iPart][nAD]    = iIDE->exitY;
                fData->exitZ[iPart][nAD]    = iIDE->exitZ;
                fData->exitT[iPart][nAD]    = iIDE->exitT;
                fData->exitPx[iPart][nAD]   = iIDE->exitMomentumX;
                fData->exitPy[iPart][nAD]   = iIDE->exitMomentumY;
                fData->exitPz[iPart][nAD]   = iIDE->exitMomentumZ;
                fData->CombinedEnergyDep[iPart][nAD] = totalE;
              }
              ++nAD;
            } // for aux det sim channels
            fData->NAuxDets[iPart] = nAD; 
            
            if (nAD > kMaxAuxDets) {
              // got this error? consider increasing kMaxAuxDets
              mf::LogError("AnalysisTree:limits") << "particle #" << iPart
                << " touches " << nAD << " auxiliary detector cells, only "
                << kMaxAuxDets << " of them are saved in the tree";
            } // if too many detector cells
          } // if (fSaveAuxDetInfo) 
        }
        else if (iPart == fData->GetMaxGEANTparticles()) {
          // got this error? it might be a bug,
          // since the structure should have enough room for everything
          mf::LogError("AnalysisTree:limits") << "event has "
            << plist.size() << " MC particles, only "
            << fData->GetMaxGEANTparticles() << " will be stored in tree";
        }     
      } // for particles
            
      fData->geant_list_size_in_tpcAV = active;
      fData->no_primaries = primary;
      fData->geant_list_size = geant_particle;
      
      LOG_DEBUG("AnalysisTree") << "Counted "
        << fData->geant_list_size << " GEANT particles ("
        << fData->geant_list_size_in_tpcAV << " in AV), "
        << fData->no_primaries << " primaries, "
        << fData->genie_no_primaries << " GENIE particles";
      
      FillWith(fData->MergedId, 0);

      // helper map track ID => index
      std::map<int, size_t> TrackIDtoIndex;
      const size_t nTrackIDs = fData->TrackId.size();
      for (size_t index = 0; index < nTrackIDs; ++index)
        TrackIDtoIndex.emplace(fData->TrackId[index], index);
      
      // for each particle, consider all the direct ancestors with the same
      // PDG ID, and mark them as belonging to the same "group"
      // (having the same MergedId)
      int currentMergedId = 1;
      for(size_t iPart = fData->geant_list_size; iPart-- > 0; ) {
        // if the particle already belongs to a group, don't bother
        if (fData->MergedId[iPart]) continue;
        // the particle starts its own group
        fData->MergedId[iPart] = currentMergedId;
        
        // look in the ancestry, one by one
        int currentMotherTrackId = fData->Mother[iPart];
        while(currentMotherTrackId > 0) {
          // find the mother (we have its track ID in currentMotherTrackId)
          std::map<int, size_t>::const_iterator iMother
            = TrackIDtoIndex.find(currentMotherTrackId);
          if (iMother == TrackIDtoIndex.end()) break; // no mother found
          size_t currentMotherIndex = iMother->second;
          // if the mother particle is of a different type,
          // don't bother with iPart ancestry any further
          if (fData->pdg[iPart] != fData->pdg[currentMotherIndex]) break;
          
          // group the "current mother" (actually, ancestor) with iPart
          fData->MergedId[currentMotherIndex] = currentMergedId;
          // then consider the grandmother (one generation earlier)
          currentMotherTrackId = fData->Mother[currentMotherIndex];
        } // while ancestry exists
        ++currentMergedId;
      } // for merging check
     } // if (fSaveGeantInfo) 
            
    }//if (mcevts_truth)
  }//if (isMC){
  fData->taulife = LArProp->ElectronLifetime();
  fTree->Fill();
  
  if (mf::isDebugEnabled()) {
    // use mf::LogDebug instead of LOG_DEBUG because we reuse it in many lines;
    // thus, we protect this part of the code with the line above
    mf::LogDebug logStream("AnalysisTreeStructure");
    logStream
      << "Tree data structure contains:"
      << "\n - " << fData->no_hits << " hits (" << fData->GetMaxHits() << ")"
      << "\n - " << fData->genie_no_primaries << " genie primaries (" << fData->GetMaxGeniePrimaries() << ")"
      << "\n - " << fData->geant_list_size << " GEANT particles (" << fData->GetMaxGEANTparticles() << "), "
        << fData->no_primaries << " primaries"
      << "\n - " << fData->geant_list_size_in_tpcAV << " GEANT particles in AV "
      << "\n - " << ((int) fData->kNTracker) << " trackers:"
      ;
    
    size_t iTracker = 0;
    for (auto tracker = fData->TrackData.cbegin();
      tracker != fData->TrackData.cend(); ++tracker, ++iTracker
    ) {
      logStream
         << "\n -> " << tracker->ntracks << " " << fTrackModuleLabel[iTracker]
           << " tracks (" << tracker->GetMaxTracks() << ")"
         ;
      for (int iTrk = 0; iTrk < tracker->ntracks; ++iTrk) {
        logStream << "\n    [" << iTrk << "] "<< tracker->ntrkhits[iTrk][0];
        for (size_t ipl = 1; ipl < tracker->GetMaxPlanesPerTrack(iTrk); ++ipl)
          logStream << " + " << tracker->ntrkhits[iTrk][ipl];
        logStream << " hits (" << tracker->GetMaxHitsPerTrack(iTrk, 0);
        for (size_t ipl = 1; ipl < tracker->GetMaxPlanesPerTrack(iTrk); ++ipl)
          logStream << " + " << tracker->GetMaxHitsPerTrack(iTrk, ipl);
        logStream << ")";
      } // for tracks
    } // for trackers
  } // if logging enabled
  
  // if we don't use a permanent buffer (which can be huge),
  // delete the current buffer, and we'll create a new one on the next event
  if (!fUseBuffer) {
    LOG_DEBUG("AnalysisTreeStructure") << "Freeing the tree data structure";
    DestroyData();
  }
} // microboone::AnalysisTree::analyze()

void microboone::AnalysisTree::HitsPurity(std::vector< art::Ptr<recob::Hit> > const& hits, Int_t& trackid, Float_t& purity, double& maxe){

  trackid = -1;
  purity = -1;

  art::ServiceHandle<cheat::BackTracker> bt;

  std::map<int,double> trkide;

  for(size_t h = 0; h < hits.size(); ++h){

    art::Ptr<recob::Hit> hit = hits[h];
    std::vector<sim::IDE> ides;
    //bt->HitToSimIDEs(hit,ides);
    std::vector<sim::TrackIDE> eveIDs = bt->HitToEveID(hit);

    for(size_t e = 0; e < eveIDs.size(); ++e){
      //std::cout<<h<<" "<<e<<" "<<eveIDs[e].trackID<<" "<<eveIDs[e].energy<<" "<<eveIDs[e].energyFrac<<std::endl;
      trkide[eveIDs[e].trackID] += eveIDs[e].energy;
    }
  }

  maxe = -1;
  double tote = 0;
  for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
    tote += ii->second;
    if ((ii->second)>maxe){
      maxe = ii->second;
      trackid = ii->first;
    }
  }

  //std::cout << "the total energy of this reco track is: " << tote << std::endl;

  if (tote>0){
    purity = maxe/tote;
  }
}

// Calculate distance to boundary.
double microboone::AnalysisTree::bdist(const TVector3& pos)
{
  // Get geometry.
  art::ServiceHandle<geo::Geometry> geom;

  double d1 = pos.X();                             // Distance to right side (wires).
  double d2 = 2.*geom->DetHalfWidth() - pos.X();   // Distance to left side (cathode).
  double d3 = pos.Y() + geom->DetHalfHeight();     // Distance to bottom.
  double d4 = geom->DetHalfHeight() - pos.Y();     // Distance to top.
  double d5 = pos.Z();                             // Distance to front.
  double d6 = geom->DetLength() - pos.Z();           // Distance to back.

  double result = std::min(std::min(std::min(std::min(std::min(d1, d2), d3), d4), d5), d6);
  return result;
}

// Length of reconstructed track, trajectory by trajectory.
double microboone::AnalysisTree::length(const recob::Track& track)
{
  double result = 0.;
  TVector3 disp = track.LocationAtPoint(0);
  int n = track.NumberTrajectoryPoints();

  for(int i = 1; i < n; ++i) {
    const TVector3& pos = track.LocationAtPoint(i);
    disp -= pos;
    result += disp.Mag();
    disp = pos;
  }
  return result;
}

// Length of MC particle, trajectory by trajectory.
double microboone::AnalysisTree::length(const simb::MCParticle& part, TVector3& start, TVector3& end)
{
  // Get geometry.
  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::DetectorProperties> detprop;

  // Get active volume boundary.
  double xmin = 0.;
  double xmax = 2.*geom->DetHalfWidth();
  double ymin = -geom->DetHalfHeight();
  double ymax = geom->DetHalfHeight();
  double zmin = 0.;
  double zmax = geom->DetLength();
  double vDrift = 160*pow(10,-6);

  double result = 0.;
  TVector3 disp;
  int n = part.NumberTrajectoryPoints();
  bool first = true;

  for(int i = 0; i < n; ++i) {
    // check if the particle is inside a TPC
   double mypos[3] = {part.Vx(i), part.Vy(i), part.Vz(i)};
   if (mypos[0] >= xmin && mypos[0] <= xmax && mypos[1] >= ymin && mypos[1] <= ymax && mypos[2] >= zmin && mypos[2] <= zmax){
     double xGen   = part.Vx(i);
     double tGen   = part.T(i);
     // Doing some manual shifting to account for
     // an interaction not occuring with the beam dump
     // we will reconstruct an x distance different from
     // where the particle actually passed to to the time
     // being different from in-spill interactions
     double newX = xGen+(tGen*vDrift);
     if (newX < -xmax || newX > (2*xmax)) continue;
     
     TVector3 pos(newX,part.Vy(i),part.Vz(i));
     if(first){
      start = pos;
     }
     else {
      disp -= pos;
      result += disp.Mag();
     }
     first = false;
     disp = pos;
     end = pos;
   }
  }
  return result;
}


namespace microboone{

  DEFINE_ART_MODULE(AnalysisTree)

}

#endif

