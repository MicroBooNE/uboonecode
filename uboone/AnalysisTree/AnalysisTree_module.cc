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
// [ ] turn the truth/GEANT information into vectors
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
// algorithm, and each has a fixed size pool for hits (do ROOT trees support
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

#include <cstring> // std::memcpy()
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional> // std::mem_fun_ref
#include <initializer_list>

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"
#include "TString.h"
#include "TTimeStamp.h"

constexpr int kNplanes       = 3;     //number of wire planes
constexpr int kMaxTrack      = 1000;  //maximum number of tracks
constexpr int kMaxClusters   = 1000;  //maximum number of clusters
constexpr int kMaxHits       = 20000; //maximum number of hits;
constexpr int kMaxPrimaries  = 20000;  //maximum number of primary particles
constexpr int kMaxTrackHits  = 1000;  //maximum number of hits on a track
constexpr int kMaxTrackers   = 10;    //number of trackers passed into fTrackModuleLabel


/// total_extent\<T\>::value has the total number of elements of an array
template <typename T>
struct total_extent {
  using value_type = size_t;
  static constexpr value_type value
    = sizeof(T) / sizeof(typename std::remove_all_extents<T>::type);
}; // total_extent<>


namespace microboone {

  class AnalysisTreeDataStruct {
      public:
    
    /// A wrapper to a C array
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
      
      /// Return the total data size of this object in bytes
      static constexpr size_t data_size() { return sizeof(array); }
      
      /// Return the total size of this object in bytes
      static constexpr size_t memory_size() { return sizeof(This_t); }
      
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
    
    class TrackDataStruct {
        public:
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
      PlaneData_t<Int_t>      trkpdgtruth; //true pdg code
      PlaneData_t<Float_t>    trkefftruth; //completeness
      PlaneData_t<Float_t>    trkpurtruth; //purity of track
      PlaneData_t<Float_t>    trkpitchc;
      PlaneData_t<Short_t>    ntrkhits;
      HitData_t<Float_t>      trkdedx;
      HitData_t<Float_t>      trkdqdx;
      HitData_t<Float_t>      trkresrg;
      HitCoordData_t<Float_t> trkxyz;

      // more track info
      TrackData_t<Short_t> trkId;
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
      
      /// Return the total size of data from dynamic vectors in bytes
      size_t dynamic_data_size() const
        {
          return trkke.size() * sizeof(typename decltype(trkke)::value_type)
            + trkrange.size() * sizeof(typename decltype(trkrange)::value_type)
            + trkidtruth.size() * sizeof(typename decltype(trkidtruth)::value_type)
            + trkpdgtruth.size() * sizeof(typename decltype(trkpdgtruth)::value_type)
            + trkefftruth.size() * sizeof(typename decltype(trkefftruth)::value_type)
            + trkpurtruth.size() * sizeof(typename decltype(trkpurtruth)::value_type)
            + trkpitchc.size() * sizeof(typename decltype(trkpitchc)::value_type)
            + ntrkhits.size() * sizeof(typename decltype(ntrkhits)::value_type)
            + trkdedx.size() * sizeof(typename decltype(trkdedx)::value_type)
            + trkdqdx.size() * sizeof(typename decltype(trkdqdx)::value_type)
            + trkresrg.size() * sizeof(typename decltype(trkresrg)::value_type)
            + trkxyz.size() * sizeof(typename decltype(trkxyz)::value_type)
            + trkId.size() * sizeof(typename decltype(trkId)::value_type)
            + trkstartx.size() * sizeof(typename decltype(trkstartx)::value_type)
            + trkstarty.size() * sizeof(typename decltype(trkstarty)::value_type)
            + trkstartz.size() * sizeof(typename decltype(trkstartz)::value_type)
            + trkstartd.size() * sizeof(typename decltype(trkstartd)::value_type)
            + trkendx.size() * sizeof(typename decltype(trkendx)::value_type)
            + trkendy.size() * sizeof(typename decltype(trkendy)::value_type)
            + trkendz.size() * sizeof(typename decltype(trkendz)::value_type)
            + trkendd.size() * sizeof(typename decltype(trkendd)::value_type)
            + trktheta.size() * sizeof(typename decltype(trktheta)::value_type)
            + trkphi.size() * sizeof(typename decltype(trkphi)::value_type)
            + trkstartdcosx.size() * sizeof(typename decltype(trkstartdcosx)::value_type)
            + trkstartdcosy.size() * sizeof(typename decltype(trkstartdcosy)::value_type)
            + trkstartdcosz.size() * sizeof(typename decltype(trkstartdcosz)::value_type)
            + trkenddcosx.size() * sizeof(typename decltype(trkenddcosx)::value_type)
            + trkenddcosy.size() * sizeof(typename decltype(trkenddcosy)::value_type)
            + trkenddcosz.size() * sizeof(typename decltype(trkenddcosz)::value_type)
            + trkthetaxz.size() * sizeof(typename decltype(trkthetaxz)::value_type)
            + trkthetayz.size() * sizeof(typename decltype(trkthetayz)::value_type)
            + trkmom.size() * sizeof(typename decltype(trkmom)::value_type)
            + trklen.size() * sizeof(typename decltype(trklen)::value_type)
           ;
        } // dynamic_data_size()
      
      /// Return the total size allocated for data from dynamic vectors in bytes
      size_t dynamic_data_capacity() const
        {
          return trkke.capacity() * sizeof(typename decltype(trkke)::value_type)
            + trkrange.capacity() * sizeof(typename decltype(trkrange)::value_type)
            + trkidtruth.capacity() * sizeof(typename decltype(trkidtruth)::value_type)
            + trkpdgtruth.capacity() * sizeof(typename decltype(trkpdgtruth)::value_type)
            + trkefftruth.capacity() * sizeof(typename decltype(trkefftruth)::value_type)
            + trkpurtruth.capacity() * sizeof(typename decltype(trkpurtruth)::value_type)
            + trkpitchc.capacity() * sizeof(typename decltype(trkpitchc)::value_type)
            + ntrkhits.capacity() * sizeof(typename decltype(ntrkhits)::value_type)
            + trkdedx.capacity() * sizeof(typename decltype(trkdedx)::value_type)
            + trkdqdx.capacity() * sizeof(typename decltype(trkdqdx)::value_type)
            + trkresrg.capacity() * sizeof(typename decltype(trkresrg)::value_type)
            + trkxyz.capacity() * sizeof(typename decltype(trkxyz)::value_type)
            + trkId.capacity() * sizeof(typename decltype(trkId)::value_type)
            + trkstartx.capacity() * sizeof(typename decltype(trkstartx)::value_type)
            + trkstarty.capacity() * sizeof(typename decltype(trkstarty)::value_type)
            + trkstartz.capacity() * sizeof(typename decltype(trkstartz)::value_type)
            + trkstartd.capacity() * sizeof(typename decltype(trkstartd)::value_type)
            + trkendx.capacity() * sizeof(typename decltype(trkendx)::value_type)
            + trkendy.capacity() * sizeof(typename decltype(trkendy)::value_type)
            + trkendz.capacity() * sizeof(typename decltype(trkendz)::value_type)
            + trkendd.capacity() * sizeof(typename decltype(trkendd)::value_type)
            + trktheta.capacity() * sizeof(typename decltype(trktheta)::value_type)
            + trkphi.capacity() * sizeof(typename decltype(trkphi)::value_type)
            + trkstartdcosx.capacity() * sizeof(typename decltype(trkstartdcosx)::value_type)
            + trkstartdcosy.capacity() * sizeof(typename decltype(trkstartdcosy)::value_type)
            + trkstartdcosz.capacity() * sizeof(typename decltype(trkstartdcosz)::value_type)
            + trkenddcosx.capacity() * sizeof(typename decltype(trkenddcosx)::value_type)
            + trkenddcosy.capacity() * sizeof(typename decltype(trkenddcosy)::value_type)
            + trkenddcosz.capacity() * sizeof(typename decltype(trkenddcosz)::value_type)
            + trkthetaxz.capacity() * sizeof(typename decltype(trkthetaxz)::value_type)
            + trkthetayz.capacity() * sizeof(typename decltype(trkthetayz)::value_type)
            + trkmom.capacity() * sizeof(typename decltype(trkmom)::value_type)
            + trklen.capacity() * sizeof(typename decltype(trklen)::value_type)
            ;
        } // dynamic_data_capacity()
      
      /// Return the total size of this object in bytes
      size_t data_size() const
        { return sizeof(MaxTracks) + sizeof(ntracks) + dynamic_data_size(); }
      
      /// Return the total size of this object in bytes
      size_t memory_size() const
        { return sizeof(*this) + dynamic_data_capacity(); }
      
    }; // class TrackDataStruct
    
    
    /// information from the run
/*    struct RunData_t {
        public:
      RunData_t() { Clear(); }
      void Clear() {}
    };
*/
    /// information from the subrun
    struct SubRunData_t {
      SubRunData_t() { Clear(); }
      void Clear() { pot = -99999.; }
      Double_t pot; //protons on target
    };

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

    //hit information
    Int_t    no_hits;                  //number of hits
    Char_t   hit_plane[kMaxHits];      //plane number
    Short_t  hit_wire[kMaxHits];       //wire number
    Short_t  hit_channel[kMaxHits];    //channel ID
    Double_t hit_peakT[kMaxHits];      //peak time
    Float_t  hit_charge[kMaxHits];     //charge (area)
    Float_t  hit_ph[kMaxHits];         //amplitude
    Short_t  hit_trkid[kMaxTrackers][kMaxHits];      //is this hit associated with a reco track?

    //track information
    Char_t   kNTracker;
    std::vector<TrackDataStruct> TrackData;
    
    //mctruth information
    Int_t     mcevts_truth;    //number of neutrino Int_teractions in the spill
    Int_t     nuPDG_truth;     //neutrino PDG code
    Int_t     ccnc_truth;      //0=CC 1=NC
    Int_t     mode_truth;      //0=QE/El, 1=RES, 2=DIS, 3=Coherent production
    Double_t  enu_truth;       //true neutrino energy
    Double_t  Q2_truth;        //Momentum transfer squared
    Double_t  W_truth;         //hadronic invariant mass
    Int_t     hitnuc_truth;    //hit nucleon
    Double_t  nuvtxx_truth;    //neutrino vertex x
    Double_t  nuvtxy_truth;    //neutrino vertex y
    Double_t  nuvtxz_truth;    //neutrino vertex z
    Double_t  nu_dcosx_truth;  //neutrino dcos x
    Double_t  nu_dcosy_truth;  //neutrino dcos y
    Double_t  nu_dcosz_truth;  //neutrino dcos z
    Double_t  lep_mom_truth;   //lepton momentum
    Double_t  lep_dcosx_truth; //lepton dcos x
    Double_t  lep_dcosy_truth; //lepton dcos y
    Double_t  lep_dcosz_truth; //lepton dcos z

    //flux information
    Double_t  tpx_flux;        //Px of parent particle leaving BNB target
    Double_t  tpy_flux;        //Py of parent particle leaving BNB target
    Double_t  tpz_flux;        //Pz of parent particle leaving BNB target
    Int_t     tptype_flux;     //Type of parent particle leaving BNB target

    //genie information
    Int_t     genie_no_primaries;
    Int_t     genie_primaries_pdg[kMaxPrimaries];
    Double_t  genie_Eng[kMaxPrimaries];
    Double_t  genie_Px[kMaxPrimaries];
    Double_t  genie_Py[kMaxPrimaries];
    Double_t  genie_Pz[kMaxPrimaries];
    Double_t  genie_P[kMaxPrimaries];
    Int_t     genie_status_code[kMaxPrimaries];
    Double_t  genie_mass[kMaxPrimaries];
    Int_t     genie_trackID[kMaxPrimaries];
    Int_t     genie_ND[kMaxPrimaries];
    Int_t     genie_mother[kMaxPrimaries];

    //geant information
    Int_t     no_primaries;      //number of primary geant particles
    Int_t     geant_list_size;  //number of all geant particles
    Int_t     geant_list_size_in_tpcFV;
    Int_t     pdg[kMaxPrimaries];
    Double_t  Eng[kMaxPrimaries];
    Double_t  Px[kMaxPrimaries];
    Double_t  Py[kMaxPrimaries];
    Double_t  Pz[kMaxPrimaries];
    Double_t  StartPointx[kMaxPrimaries];
    Double_t  StartPointy[kMaxPrimaries];
    Double_t  StartPointz[kMaxPrimaries];
    Double_t  EndPointx[kMaxPrimaries];
    Double_t  EndPointy[kMaxPrimaries];
    Double_t  EndPointz[kMaxPrimaries];
    Int_t     NumberDaughters[kMaxPrimaries];
    Int_t     TrackId[kMaxPrimaries];
    Int_t     Mother[kMaxPrimaries];
    Int_t     process_primary[kMaxPrimaries];
    Int_t     MergedId[kMaxPrimaries]; //geant track segments, which belong to the same particle, get the same

    // more geant information
    Int_t   geant_tpcFV_status[kMaxPrimaries];
    Int_t   geant_tpcFV_trackId[kMaxPrimaries];
    Int_t   geant_tpcFV_pdg[kMaxPrimaries];

    Double_t  geant_tpcFV_orig_E[kMaxPrimaries];
    Double_t  geant_tpcFV_orig_px[kMaxPrimaries];
    Double_t  geant_tpcFV_orig_py[kMaxPrimaries];
    Double_t  geant_tpcFV_orig_pz[kMaxPrimaries];
    Double_t  geant_tpcFV_orig_startx[kMaxPrimaries];
    Double_t  geant_tpcFV_orig_starty[kMaxPrimaries];
    Double_t  geant_tpcFV_orig_startz[kMaxPrimaries];
    Double_t  geant_tpcFV_orig_startt[kMaxPrimaries];
    Double_t  geant_tpcFV_orig_endx[kMaxPrimaries];
    Double_t  geant_tpcFV_orig_endy[kMaxPrimaries];
    Double_t  geant_tpcFV_orig_endz[kMaxPrimaries];
    Double_t  geant_tpcFV_orig_endt[kMaxPrimaries];

    Double_t  geant_tpcFV_startx[kMaxPrimaries];          // starting x position.
    Double_t  geant_tpcFV_starty[kMaxPrimaries];          // starting y position.
    Double_t  geant_tpcFV_startz[kMaxPrimaries];          // starting z position.
    Double_t  geant_tpcFV_startd[kMaxPrimaries];          // starting distance to boundary.
    Double_t  geant_tpcFV_endx[kMaxPrimaries];          // ending x position.
    Double_t  geant_tpcFV_endy[kMaxPrimaries];          // ending y position.
    Double_t  geant_tpcFV_endz[kMaxPrimaries];          // ending z position.
    Double_t  geant_tpcFV_endd[kMaxPrimaries];          // ending distance to boundary.
    Double_t  geant_tpcFV_theta[kMaxPrimaries];          // theta.
    Double_t  geant_tpcFV_phi[kMaxPrimaries];          // phi.
    Double_t  geant_tpcFV_theta_xz[kMaxPrimaries];    // theta_xz.
    Double_t  geant_tpcFV_theta_yz[kMaxPrimaries];    // theta_yz.
    Double_t  geant_tpcFV_mom[kMaxPrimaries];         // momentum.
    Double_t  geant_tpcFV_len[kMaxPrimaries];         // length.

    /// Constructor; clears all fields
    AnalysisTreeDataStruct(size_t nTrackers = 0)
      { SetTrackers(nTrackers); Clear(); }

    TrackDataStruct& GetTrackerData(size_t iTracker)
      { return TrackData.at(iTracker); }
    const TrackDataStruct& GetTrackerData(size_t iTracker) const
      { return TrackData.at(iTracker); }
    
    /// Clear all fields
    void Clear();
    
    /// Allocates data structures for the given number of trackers (no Clear())
    void SetTrackers(size_t nTrackers) { TrackData.resize(nTrackers); }

    /// Returns the number of trackers for which data structures are allocated
    size_t GetNTrackers() const { return TrackData.size(); }
    
    /// Returns the number of hits for which memory is allocated
    size_t GetMaxHits() const { return kMaxHits; }
    
    /// Returns the number of trackers for which memory is allocated
    size_t GetMaxTrackers() const { return TrackData.capacity(); }
    
    /// Returns the number of GEANT primaries for which memory is allocated
    size_t GetMaxGEANTePrimaries() const { return kMaxPrimaries; }
    
    /// Returns the number of GENIE primaries for which memory is allocated
    size_t GetMaxGeniePrimaries() const { return kMaxPrimaries; }
    
    /// Connect this object with a tree
    void SetAddresses(TTree* pTree, const std::vector<std::string>& trackers);

      /// Return the total size of data from dynamic vectors in bytes
      size_t dynamic_data_size() const
        {
          size_t total = 0;
          for(const auto& TrackerData: TrackData)
            total += TrackerData.dynamic_data_size();
          return total;
        } // dynamic_data_size()
      
      /// Return the total size allocated for data from dynamic vectors in bytes
      size_t dynamic_data_capacity() const
        {
          size_t total = 0;
          for(const auto& TrackerData: TrackData)
            total += TrackerData.dynamic_data_capacity();
          return total;
        } // dynamic_data_capacity()
      
      /// Return the total size of this object in bytes
      size_t data_size() const
        { return sizeof(*this) + dynamic_data_size(); }
      
      /// Return the total size of this object in bytes
      size_t memory_size() const
        { return sizeof(*this) + dynamic_data_capacity(); }
      
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
    std::string fG4ModuleLabel;
    std::string fClusterModuleLabel;
    std::string fKingaModuleLabel;
    std::string fLineMergerModuleLabel;
    std::string fDbscanModuleLabel;
    std::string fFuzzyModuleLabel;
    std::vector<std::string> fTrackModuleLabel;
    std::string fEndPoint2DModuleLabel;
    std::string fVertexModuleLabel;
    std::string fPOTModuleLabel;
    std::vector<std::string> fCalorimetryModuleLabel;
    std::string fParticleIDModuleLabel;
    bool fUseBuffer; ///< whether to use a permanent buffer (faster, huge memory)

    size_t GetNTrackers() const { return fTrackModuleLabel.size(); }
    
    // make sure the data structure exists
    void CreateData(bool bClearData = false)
      {
        if (!fData) fData = new AnalysisTreeDataStruct(GetNTrackers());
        else {
          fData->SetTrackers(GetNTrackers());
          if (bClearData) fData->Clear();
        }
      } // CreateData()
    void SetAddresses()
      {
        CheckData("SetAddress()"); CheckTree("SetAddress()");
        fData->SetAddresses(fTree, fTrackModuleLabel);
      } // SetAddresses()
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
    void UpdateAddresses() { CreateData(); SetAddresses(); }
    
    /// Create the output tree and the data structures, if needed
    void CreateTree(bool bClearData = false);
    
    /// Destroy the local buffers (existing branches will point to invalid address!)
    void DestroyData() { if (fData) { delete fData; fData = nullptr; } }
    
    void CheckData(std::string caller) const
      {
        if (fData) return;
        throw art::Exception(art::errors::LogicError)
          << "AnalysisTree::" << caller << ": no data";
      } // CheckData()
    void CheckTree(std::string caller) const
      {
        if (fTree) return;
        throw art::Exception(art::errors::LogicError)
          << "AnalysisTree::" << caller << ": no tree";
      } // CheckData()
  }; // class microboone::AnalysisTree
} // namespace microboone


namespace {
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
    { FillWith(data.begin(), data.end(), value); }

} // local namespace


//------------------------------------------------------------------------------
//---  AnalysisTreeDataStruct::TrackDataStruct
//---

void microboone::AnalysisTreeDataStruct::TrackDataStruct::Resize(size_t nTracks)
{
  MaxTracks = nTracks;
  
  trkId.resize(MaxTracks);
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
  trklen.resize(MaxTracks);
  
  trkke.resize(MaxTracks);
  trkrange.resize(MaxTracks);
  trkidtruth.resize(MaxTracks);
  trkpdgtruth.resize(MaxTracks);
  trkefftruth.resize(MaxTracks);
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
  FillWith(trklen       , -99999.);
  
  for (size_t iTrk = 0; iTrk < MaxTracks; ++iTrk){
    
    // the following are BoxedArray's;
    // their iterators traverse all the array dimensions
    FillWith(trkke[iTrk]      , -99999.);
    FillWith(trkrange[iTrk]   , -99999.);
    FillWith(trkidtruth[iTrk] , -99999 );
    FillWith(trkpdgtruth[iTrk], -99999 );
    FillWith(trkefftruth[iTrk], -99999.);
    FillWith(trkpurtruth[iTrk], -99999.);
    FillWith(trkpitchc[iTrk]  , -99999.);
    FillWith(ntrkhits[iTrk]   ,  -9999 );
    
    FillWith(trkdedx[iTrk], 0.);
    FillWith(trkdqdx[iTrk], 0.);
    FillWith(trkresrg[iTrk], 0.);
    
    FillWith(trkxyz[iTrk], 0.);
    
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
  
  BranchName = "trkke_" + TrackLabel;
  CreateBranch(BranchName, trkke, BranchName + NTracksIndexStr + "[3]/F");
  
  BranchName = "trkrange_" + TrackLabel;
  CreateBranch(BranchName, trkrange, BranchName + NTracksIndexStr + "[3]/F");
  
  BranchName = "trkidtruth_" + TrackLabel;
  CreateBranch(BranchName, trkidtruth, BranchName + NTracksIndexStr + "[3]/I");
  
  BranchName = "trkpdgtruth_" + TrackLabel;
  CreateBranch(BranchName, trkpdgtruth, BranchName + NTracksIndexStr + "[3]/I");
  
  BranchName = "trkefftruth_" + TrackLabel;
  CreateBranch(BranchName, trkefftruth, BranchName + NTracksIndexStr + "[3]/F");
  
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
  
  BranchName = "trklen_" + TrackLabel;
  CreateBranch(BranchName, trklen, BranchName + NTracksIndexStr + "/F");
} // microboone::AnalysisTreeDataStruct::TrackDataStruct::SetAddresses()

//------------------------------------------------------------------------------
//---  AnalysisTreeDataStruct
//---

void microboone::AnalysisTreeDataStruct::Clear() {

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
  
  std::fill(hit_plane, hit_plane + sizeof(hit_plane)/sizeof(hit_plane[0]), -99);
  std::fill(hit_wire, hit_wire + sizeof(hit_wire)/sizeof(hit_wire[0]), -9999);
  std::fill(hit_channel, hit_channel + sizeof(hit_channel)/sizeof(hit_channel[0]), -9999);
  std::fill(hit_peakT, hit_peakT + sizeof(hit_peakT)/sizeof(hit_peakT[0]), -99999.);
  std::fill(hit_charge, hit_charge + sizeof(hit_charge)/sizeof(hit_charge[0]), -99999.);
  std::fill(hit_ph, hit_ph + sizeof(hit_ph)/sizeof(hit_ph[0]), -99999.);

  for (size_t iTrk = 0; iTrk < kMaxTrackers; ++iTrk) {
    std::fill(hit_trkid[iTrk], hit_trkid[iTrk] + kMaxHits, -9999);
  }
  
  std::for_each
    (TrackData.begin(), TrackData.end(), std::mem_fun_ref(&TrackDataStruct::Clear));

  mcevts_truth = -99999;
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
  no_primaries = 0;
  geant_list_size=0;
  geant_list_size_in_tpcFV = 0;
  for (int i = 0; i<kMaxPrimaries; ++i){
    pdg[i] = -99999;
    Eng[i] = -99999;
    Px[i] = -99999;
    Py[i] = -99999;
    Pz[i] = -99999;
    StartPointx[i] = -99999;
    StartPointy[i] = -99999;
    StartPointz[i] = -99999;
    EndPointx[i] = -99999;
    EndPointy[i] = -99999;
    EndPointz[i] = -99999;
    NumberDaughters[i] = -99999;
    Mother[i] = -99999;
    TrackId[i] = -99999;
    process_primary[i] = -99999;
    genie_primaries_pdg[i] = -99999;
    genie_Eng[i] = -99999;
    genie_Px[i] = -99999;
    genie_Py[i] = -99999;
    genie_Pz[i] = -99999;
    genie_P[i] = -99999;
    genie_status_code[i] = -99999;
    genie_mass[i] = -99999;
    genie_trackID[i] = -99999;
    genie_ND[i] = -99999;
    genie_mother[i] = -99999;

    geant_tpcFV_status[i] = -99999;
    geant_tpcFV_trackId[i] = -99999;
    geant_tpcFV_pdg[i] = -99999;

    geant_tpcFV_orig_E[i] = -99999;
    geant_tpcFV_orig_px[i] = -99999;
    geant_tpcFV_orig_py[i] = -99999;
    geant_tpcFV_orig_pz[i] = -99999;
    geant_tpcFV_orig_startx[i] = -99999;
    geant_tpcFV_orig_starty[i] = -99999;
    geant_tpcFV_orig_startz[i] = -99999;
    geant_tpcFV_orig_startt[i] = -99999;
    geant_tpcFV_orig_endx[i] = -99999;
    geant_tpcFV_orig_endy[i] = -99999;
    geant_tpcFV_orig_endz[i] = -99999;
    geant_tpcFV_orig_endt[i] = -99999;

    geant_tpcFV_startx[i] = -99999;
    geant_tpcFV_starty[i] = -99999;
    geant_tpcFV_startz[i] = -99999;
    geant_tpcFV_startd[i] = -99999;
    geant_tpcFV_endx[i] = -99999;
    geant_tpcFV_endy[i] = -99999;
    geant_tpcFV_endz[i] = -99999;
    geant_tpcFV_endd[i] = -99999;
    geant_tpcFV_theta[i] = -99999;
    geant_tpcFV_phi[i] = -99999;
    geant_tpcFV_theta_xz[i] = -99999;
    geant_tpcFV_theta_yz[i] = -99999;
    geant_tpcFV_mom[i] = -99999;
    geant_tpcFV_len[i] = -99999;
  }

} // microboone::AnalysisTreeDataStruct::Clear()


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

  CreateBranch("no_hits",&no_hits,"no_hits/I");
  CreateBranch("hit_plane",hit_plane,"hit_plane[no_hits]/B");
  CreateBranch("hit_wire",hit_wire,"hit_wire[no_hits]/S");
  CreateBranch("hit_channel",hit_channel,"hit_channel[no_hits]/S");
  CreateBranch("hit_peakT",hit_peakT,"hit_peakT[no_hits]/D");
  CreateBranch("hit_charge",hit_charge,"hit_charge[no_hits]/F");
  CreateBranch("hit_ph",hit_ph,"hit_ph[no_hits]/F");

  AutoResettingStringSteam sstr;
  sstr() << kMaxTrackHits;
  std::string MaxTrackHitsIndexStr("[" + sstr.str() + "]");

  kNTracker = trackers.size();
  CreateBranch("kNTracker",&kNTracker,"kNTracker/B");
  for(int i=0; i<kNTracker; i++){
    std::string TrackLabel = trackers[i];
    std::string BranchName;

    BranchName = "hit_trkid_" + TrackLabel;
    CreateBranch(BranchName, hit_trkid[i], BranchName + "[no_hits]/S");

    // note that if the tracker data has maximum number of tracks 0,
    // nothing is initialized (branches are not even created)
    TrackData[i].SetAddresses(pTree, TrackLabel);
    
  } // for trackers

  CreateBranch("mcevts_truth",&mcevts_truth,"mcevts_truth/I");
  CreateBranch("nuPDG_truth",&nuPDG_truth,"nuPDG_truth/I");
  CreateBranch("ccnc_truth",&ccnc_truth,"ccnc_truth/I");
  CreateBranch("mode_truth",&mode_truth,"mode_truth/I");
  CreateBranch("enu_truth",&enu_truth,"enu_truth/D");
  CreateBranch("Q2_truth",&Q2_truth,"Q2_truth/D");
  CreateBranch("W_truth",&W_truth,"W_truth/D");
  CreateBranch("hitnuc_truth",&hitnuc_truth,"hitnuc_truth/I");
  CreateBranch("nuvtxx_truth",&nuvtxx_truth,"nuvtxx_truth/D");
  CreateBranch("nuvtxy_truth",&nuvtxy_truth,"nuvtxy_truth/D");
  CreateBranch("nuvtxz_truth",&nuvtxz_truth,"nuvtxz_truth/D");
  CreateBranch("nu_dcosx_truth",&nu_dcosx_truth,"nu_dcosx_truth/D");
  CreateBranch("nu_dcosy_truth",&nu_dcosy_truth,"nu_dcosy_truth/D");
  CreateBranch("nu_dcosz_truth",&nu_dcosz_truth,"nu_dcosz_truth/D");
  CreateBranch("lep_mom_truth",&lep_mom_truth,"lep_mom_truth/D");
  CreateBranch("lep_dcosx_truth",&lep_dcosx_truth,"lep_dcosx_truth/D");
  CreateBranch("lep_dcosy_truth",&lep_dcosy_truth,"lep_dcosy_truth/D");
  CreateBranch("lep_dcosz_truth",&lep_dcosz_truth,"lep_dcosz_truth/D");

  CreateBranch("tpx_flux",&tpx_flux,"tpx_flux/D");
  CreateBranch("tpy_flux",&tpy_flux,"tpy_flux/D");
  CreateBranch("tpz_flux",&tpz_flux,"tpz_flux/D");
  CreateBranch("tptype_flux",&tptype_flux,"tptype_flux/I");

  CreateBranch("genie_no_primaries",&genie_no_primaries,"genie_no_primaries/I");
  CreateBranch("genie_primaries_pdg",genie_primaries_pdg,"genie_primaries_pdg[genie_no_primaries]/I");
  CreateBranch("genie_Eng",genie_Eng,"genie_Eng[genie_no_primaries]/D");
  CreateBranch("genie_Px",genie_Px,"genie_Px[genie_no_primaries]/D");
  CreateBranch("genie_Py",genie_Py,"genie_Py[genie_no_primaries]/D");
  CreateBranch("genie_Pz",genie_Pz,"genie_Pz[genie_no_primaries]/D");
  CreateBranch("genie_P",genie_P,"genie_P[genie_no_primaries]/D");
  CreateBranch("genie_status_code",genie_status_code,"genie_status_code[genie_no_primaries]/I");
  CreateBranch("genie_mass",genie_mass,"genie_mass[genie_no_primaries]/D");
  CreateBranch("genie_trackID",genie_trackID,"genie_trackID[genie_no_primaries]/I");
  CreateBranch("genie_ND",genie_ND,"genie_ND[genie_no_primaries]/I");
  CreateBranch("genie_mother",genie_mother,"genie_mother[genie_no_primaries]/I");

  CreateBranch("no_primaries",&no_primaries,"no_primaries/I");
  CreateBranch("geant_list_size",&geant_list_size,"geant_list_size/I");

  CreateBranch("pdg",pdg,"pdg[geant_list_size]/I");
  CreateBranch("Eng",Eng,"Eng[geant_list_size]/D");
  CreateBranch("Px",Px,"Px[geant_list_size]/D");
  CreateBranch("Py",Py,"Py[geant_list_size]/D");
  CreateBranch("Pz",Pz,"Pz[geant_list_size]/D");
  CreateBranch("StartPointx",StartPointx,"StartPointx[geant_list_size]/D");
  CreateBranch("StartPointy",StartPointy,"StartPointy[geant_list_size]/D");
  CreateBranch("StartPointz",StartPointz,"StartPointz[geant_list_size]/D");
  CreateBranch("EndPointx",EndPointx,"EndPointx[geant_list_size]/D");
  CreateBranch("EndPointy",EndPointy,"EndPointy[geant_list_size]/D");
  CreateBranch("EndPointz",EndPointz,"EndPointz[geant_list_size]/D");
  CreateBranch("NumberDaughters",NumberDaughters,"NumberDaughters[geant_list_size]/I");
  CreateBranch("Mother",Mother,"Mother[geant_list_size]/I");
  CreateBranch("TrackId",TrackId,"TrackId[geant_list_size]/I");
  CreateBranch("MergedId", MergedId, "MergedId[geant_list_size]/I");
  CreateBranch("process_primary",process_primary,"process_primary[geant_list_size]/I");

  CreateBranch("geant_list_size_in_tpcFV",&geant_list_size_in_tpcFV,"geant_list_size_in_tpcFV/I");
  CreateBranch("geant_tpcFV_pdg", geant_tpcFV_pdg, "geant_tpcFV_pdg[geant_list_size_in_tpcFV]/I");
  CreateBranch("geant_tpcFV_status", geant_tpcFV_status, "geant_tpcFV_status[geant_list_size_in_tpcFV]/I");
  CreateBranch("geant_tpcFV_trackId", geant_tpcFV_trackId, "geant_tpcFV_trackId[geant_list_size_in_tpcFV]/I");
  CreateBranch("geant_tpcFV_orig_E", geant_tpcFV_orig_E, "geant_tpcFV_orig_E[geant_list_size_in_tpcFV]/F");
  CreateBranch("geant_tpcFV_orig_px", geant_tpcFV_orig_px, "geant_tpcFV_orig_px[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_orig_py", geant_tpcFV_orig_py, "geant_tpcFV_orig_py[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_orig_pz", geant_tpcFV_orig_pz, "geant_tpcFV_orig_pz[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_orig_startx", geant_tpcFV_orig_startx, "geant_tpcFV_orig_startx[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_orig_starty", geant_tpcFV_orig_starty, "geant_tpcFV_orig_starty[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_orig_startz", geant_tpcFV_orig_startz, "geant_tpcFV_orig_startz[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_orig_startt", geant_tpcFV_orig_startt, "geant_tpcFV_orig_startt[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_orig_endx", geant_tpcFV_orig_endx, "geant_tpcFV_orig_endx[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_orig_endy", geant_tpcFV_orig_endy, "geant_tpcFV_orig_endy[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_orig_endz", geant_tpcFV_orig_endz, "geant_tpcFV_orig_endz[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_orig_endt", geant_tpcFV_orig_endt, "geant_tpcFV_orig_endt[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_startx", geant_tpcFV_startx, "geant_tpcFV_startx[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_starty", geant_tpcFV_starty, "geant_tpcFV_starty[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_startz", geant_tpcFV_startz, "geant_tpcFV_startz[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_startd", geant_tpcFV_startd, "geant_tpcFV_startd[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_endx", geant_tpcFV_endx, "geant_tpcFV_endx[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_endy", geant_tpcFV_endy, "geant_tpcFV_endy[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_endz", geant_tpcFV_endz, "geant_tpcFV_endz[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_endd", geant_tpcFV_endd, "geant_tpcFV_endd[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_theta", geant_tpcFV_theta, "geant_tpcFV_theta[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_phi", geant_tpcFV_phi, "geant_tpcFV_phi[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_theta_xz", geant_tpcFV_theta_xz, "geant_tpcFV_theta_xz[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_theta_yz", geant_tpcFV_theta_yz, "geant_tpcFV_theta_yz[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_mom", geant_tpcFV_mom, "geant_tpcFV_mom[geant_list_size_in_tpcFV]/D");
  CreateBranch("geant_tpcFV_len", geant_tpcFV_len, "geant_tpcFV_len[geant_list_size_in_tpcFV]/D");
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
  fG4ModuleLabel            (pset.get< std::string >("G4ModuleLabel")           ),
  fClusterModuleLabel       (pset.get< std::string >("ClusterModuleLabel")      ),
  fKingaModuleLabel         (pset.get< std::string >("KingaModuleLabel")        ),
  fLineMergerModuleLabel    (pset.get< std::string >("LineMergerModuleLabel")   ),
  fDbscanModuleLabel        (pset.get< std::string >("DbscanModuleLabel")       ),
  fFuzzyModuleLabel         (pset.get< std::string >("FuzzyModuleLabel")       ),
  fTrackModuleLabel         (pset.get< std::vector<std::string> >("TrackModuleLabel")),
  fEndPoint2DModuleLabel    (pset.get< std::string >("EndPoint2DModuleLabel")   ),
  fVertexModuleLabel        (pset.get< std::string >("VertexModuleLabel")       ),
  fPOTModuleLabel           (pset.get< std::string >("POTModuleLabel")          ),
  fCalorimetryModuleLabel   (pset.get< std::vector<std::string> >("CalorimetryModuleLabel")),
  fParticleIDModuleLabel    (pset.get< std::string >("ParticleIDModuleLabel")   ),
  fUseBuffer                (pset.get< bool >("UseBuffers", false)              )
{
  mf::LogInfo("AnalysisTree") << "Configuration:"
    << "\n  UseBuffers: " << std::boolalpha << fUseBuffer
    ;
}

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
  // collect the sizes which might me needed to resize the tree data structure:
  
  // hits
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);

//  const size_t Nplanes       = 3; // number of wire planes; pretty much constant...
  const size_t NTrackers = GetNTrackers(); // number of trackers passed into fTrackModuleLabel
  const size_t NHits     = hitlist.size(); // number of hits
  
  // make sure there is the data, the tree and everything
  CreateTree(true);

  mf::LogDebug("AnalysisTreeStructure") << "After initialization, tree data structure has "
    << fData->data_size() << " bytes in data, " << fData->memory_size()
    << " allocated";
  
  /// transfer the run and subrun data to the tree data object
//  fData->RunData = RunData;
  fData->SubRunData = SubRunData;

  std::vector< art::Handle< std::vector<recob::Track> > > trackListHandle(NTrackers);
  std::vector< std::vector<art::Ptr<recob::Track> > > tracklist(NTrackers);
  for (unsigned int it = 0; it < NTrackers; ++it){
    if (evt.getByLabel(fTrackModuleLabel[it],trackListHandle[it]))
      art::fill_ptr_vector(tracklist[it], trackListHandle[it]);
  }

  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
    art::fill_ptr_vector(mclist, mctruthListHandle);

  art::Handle< std::vector<simb::MCFlux> > mcfluxListHandle;
  std::vector<art::Ptr<simb::MCFlux> > fluxlist;
  if (evt.getByLabel(fGenieGenModuleLabel,mcfluxListHandle))
    art::fill_ptr_vector(fluxlist, mcfluxListHandle);

  //services
  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<cheat::BackTracker> bt;
  art::ServiceHandle<util::DetectorProperties> detprop;
  art::ServiceHandle<util::LArProperties> LArProp;

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
  fData->isdata = evt.isRealData()? 1: 0;

  //hit information
  fData->no_hits = (int) NHits;
  if (NHits > kMaxHits) {
    // got this error? consider increasing kMaxHits
    // (or ask for a redesign using vectors)
    mf::LogError("AnalysisTree:limits") << "event has " << NHits
      << " hits, only kMaxHits=" << kMaxHits << " stored in tree";
  }
  for (size_t i = 0; i< NHits; ++i){//loop over hits
    fData->hit_channel[i] = hitlist[i]->Channel();
    fData->hit_plane[i]   = hitlist[i]->WireID().Plane;
    fData->hit_wire[i]    = hitlist[i]->WireID().Wire;
    fData->hit_peakT[i]   = hitlist[i]->PeakTime();
    fData->hit_charge[i]  = hitlist[i]->Charge();
    fData->hit_charge[i]  = hitlist[i]->Charge(true);
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
  }

  //track information for multiple trackers
  for (unsigned int iTracker=0; iTracker < NTrackers; ++iTracker){
    AnalysisTreeDataStruct::TrackDataStruct& TrackerData = fData->GetTrackerData(iTracker);
    
    size_t NTracks = tracklist[iTracker].size();
    TrackerData.ntracks = (int) NTracks;
    // allocate enough space for this number of tracks (but at least for one of them!)
    TrackerData.SetMaxTracks(std::max(NTracks, (size_t) 1));
    // now set the tree addresses to the newly allocated memory;
    // this creates the tree branches in case they are not there yet
    SetTrackerAddresses(iTracker);
    if (NTracks > TrackerData.GetMaxTracks()) {
      // got this error? consider increasing kMaxTrack
      // (or ask for a redesign using vectors)
      mf::LogError("AnalysisTree:limits") << "event has " << NTracks
        << " " << fTrackModuleLabel[iTracker] << " tracks, only "
        << TrackerData.GetMaxTracks() << " stored in tree";
    }
    for(size_t iTrk=0; iTrk < NTracks; ++iTrk){//loop over tracks

      art::Ptr<recob::Track> ptrack(trackListHandle[iTracker], iTrk);
      const recob::Track& track = *ptrack;
      
      TVector3 pos, dir_start, dir_end, end;
      double tlen = 0., mom = 0.;
      
      //we need to use Bezier methods for Bezier tracks
      if(fTrackModuleLabel[iTracker].compare("beziertracker")==0){
          trkf::BezierTrack btrack(*ptrack);
          int ntraj = btrack.NSegments();
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

            tlen        = btrack.GetLength();
            if (btrack.NumberFitMomentum() > 0)
              mom = btrack.VertexMomentum();
            // fill bezier track reco branches
            TrackerData.trkId[iTrk]         = iTrk;  //bezier has some screwed up track IDs
        }
      }
      else {   //use the normal methods for other kinds of tracks
        int ntraj = track.NumberTrajectoryPoints();
        if (ntraj > 0) {
            pos       = track.Vertex();
            dir_start = track.VertexDirection();
            dir_end   = track.EndDirection();
            end       = track.End();

            tlen        = length(track);
            if(track.NumberFitMomentum() > 0)
              mom = track.VertexMomentum();
            // fill non-bezier-track reco branches
            TrackerData.trkId[iTrk]                    = track.ID();
        }
      }
      double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
      double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());
      double dpos = bdist(pos);
      double dend = bdist(end);
      
      TrackerData.trkstartx[iTrk]             = pos.X();
      TrackerData.trkstarty[iTrk]             = pos.Y();
      TrackerData.trkstartz[iTrk]             = pos.Z();
      TrackerData.trkstartd[iTrk]             = dpos;
      TrackerData.trkendx[iTrk]               = end.X();
      TrackerData.trkendy[iTrk]               = end.Y();
      TrackerData.trkendz[iTrk]               = end.Z();
      TrackerData.trkendd[iTrk]               = dend;
      TrackerData.trktheta[iTrk]              = dir_start.Theta();
      TrackerData.trkphi[iTrk]                = dir_start.Phi();
      TrackerData.trkstartdcosx[iTrk]         = dir_start.X();
      TrackerData.trkstartdcosy[iTrk]         = dir_start.Y();
      TrackerData.trkstartdcosz[iTrk]         = dir_start.Z();
      TrackerData.trkenddcosx[iTrk]           = dir_end.X();
      TrackerData.trkenddcosy[iTrk]           = dir_end.Y();
      TrackerData.trkenddcosz[iTrk]           = dir_end.Z();
      TrackerData.trkthetaxz[iTrk]            = theta_xz;
      TrackerData.trkthetayz[iTrk]            = theta_yz;
      TrackerData.trkmom[iTrk]                = mom;
      TrackerData.trklen[iTrk]                = tlen;

      art::FindMany<anab::Calorimetry> fmcal(trackListHandle[iTracker], evt, fCalorimetryModuleLabel[iTracker]);
      if (fmcal.isValid()){
        std::vector<const anab::Calorimetry*> calos = fmcal.at(iTrk);
        //std::cout<<"calo size "<<calos.size()<<std::endl;
        if (calos.size() > TrackerData.GetMaxPlanesPerTrack(iTrk)) {
          // if you get this message, there is probably a bug somewhere since
          // the calorimetry planes should be 3.
          mf::LogError("AnalysisTree:limits")
            << "the " << fTrackModuleLabel[iTracker] << " track #" << iTrk
            << " has " << calos.size() << " planes for calorimetry , only "
            << TrackerData.GetMaxPlanesPerTrack(iTrk) << " stored in tree";
        }
        for (size_t ipl = 0; ipl<calos.size(); ++ipl){
          TrackerData.trkke[iTrk][ipl]    = calos[ipl]->KineticEnergy();
          TrackerData.trkrange[iTrk][ipl] = calos[ipl]->Range();
          TrackerData.trkpitchc[iTrk][ipl]= calos[ipl] -> TrkPitchC();
          const size_t NHits = calos[ipl] -> dEdx().size();
          TrackerData.ntrkhits[iTrk][ipl] = (int) NHits;
          if (NHits > TrackerData.GetMaxHitsPerTrack(iTrk)) {
            // if you get this error, you'll have to increase kMaxTrackHits
            mf::LogError("AnalysisTree:limits")
              << "the " << fTrackModuleLabel[iTracker] << " track #" << iTrk
              << " has " << NHits << " hits on calorimetry plane #" << ipl
              <<", only "
              << TrackerData.GetMaxHitsPerTrack(iTrk, ipl) << " stored in tree";
          }
          for(size_t iTrkHit = 0; iTrkHit < NHits; ++iTrkHit) {
            TrackerData.trkdedx[iTrk][ipl][iTrkHit]  = (calos[ipl] -> dEdx())[iTrkHit];
            TrackerData.trkdqdx[iTrk][ipl][iTrkHit]  = (calos[ipl] -> dQdx())[iTrkHit];
            TrackerData.trkresrg[iTrk][ipl][iTrkHit] = (calos[ipl] -> ResidualRange())[iTrkHit];
            const auto& TrkPos = (calos[ipl] -> XYZ())[iTrkHit];
            auto& TrkXYZ = TrackerData.trkxyz[iTrk][ipl][iTrkHit];
            TrkXYZ[0] = TrkPos.X();
            TrkXYZ[1] = TrkPos.Y();
            TrkXYZ[2] = TrkPos.Z();
          } // for track hits
        } // for calorimetry info
      } // if has calorimetry info

      //track truth information
      if (!fData->isdata){
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
            const simb::MCParticle *particle = bt->TrackIDToParticle(TrackerData.trkidtruth[iTrk][ipl]);
            const std::vector<sim::IDE> vide = bt->TrackIDToSimIDE(TrackerData.trkidtruth[iTrk][ipl]);
            double tote = 0;
            for (size_t iide = 0; iide<vide.size(); ++iide){
              tote += vide[iide].energy;
            }
            TrackerData.trkpdgtruth[iTrk][ipl] = particle->PdgCode();
            TrackerData.trkefftruth[iTrk][ipl] = maxe/(tote/kNplanes); //tote include both induction and collection energies
          //std::cout<<"\n"<<trkpdgtruth[iTracker][iTrk][ipl]<<"\t"<<trkefftruth[iTracker][iTrk][ipl];
          }
        }
      }//end if (!isdata)
    }//end loop over track
  }//end loop over track module labels

  //mc truth information
  if (!fData->isdata){//is MC
    const sim::ParticleList& plist = bt->ParticleList();
    //save neutrino interaction information
    fData->mcevts_truth = mclist.size();
    if (fData->mcevts_truth){//at least one mc record
      //if (mclist[0]->NeutrinoSet()){//is neutrino
      //sometimes there can be multiple neutrino interactions in one spill
      //this is trying to find the primary interaction
      //by looking for the highest energy deposition
      std::map<art::Ptr<simb::MCTruth>,double> mctruthemap;
      for (size_t i = 0; i<hitlist.size(); i++){
        //if (hitlist[i]->View() == geo::kV){//collection view
        std::vector<cheat::TrackIDE> eveIDs = bt->HitToEveID(hitlist[i]);
        for (size_t e = 0; e<eveIDs.size(); e++){
          art::Ptr<simb::MCTruth> mctruth = bt->TrackIDToMCTruth(eveIDs[e].trackID);
          mctruthemap[mctruth]+=eveIDs[e].energy;
        }
        //}
      }
      art::Ptr<simb::MCTruth> mctruth = mclist[0];
      double maxenergy = -1;
      int imc = 0;
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
        fData->genie_no_primaries=mctruth->NParticles();

        for(int j = 0; j < mctruth->NParticles(); ++j){
          simb::MCParticle part(mctruth->GetParticle(j));
          
          fData->genie_primaries_pdg[j]=part.PdgCode();
          fData->genie_Eng[j]=part.E();
          fData->genie_Px[j]=part.Px();
          fData->genie_Py[j]=part.Py();
          fData->genie_Pz[j]=part.Pz();
          fData->genie_P[j]=part.Px();
          fData->genie_status_code[j]=part.StatusCode();
          fData->genie_mass[j]=part.Mass();
          fData->genie_trackID[j]=part.TrackId();
          fData->genie_ND[j]=part.NumberDaughters();
          fData->genie_mother[j]=part.Mother();
        }
      }

      //GEANT particles information
      std::vector<const simb::MCParticle* > geant_part;
      int i = 0;
      fData->geant_list_size_in_tpcFV = 0;
      for(size_t p = 0; p < plist.size(); ++p)
      {
        geant_part.push_back(plist.Particle(p));
        assert(plist.Particle(p) != 0);
        TVector3 mcstart, mcend;
        double plen = length(*plist.Particle(p), mcstart, mcend);
        if ( (plen==0) || plist.Particle(p)->PdgCode() > 10000) continue;
        else{
          fData->geant_tpcFV_pdg[i] = plist.Particle(p)->PdgCode();
          double mctheta_xz = std::atan2(plist.Particle(p)->Px(), plist.Particle(p)->Pz());
          double mctheta_yz = std::atan2(plist.Particle(p)->Py(), plist.Particle(p)->Pz());

          fData->geant_tpcFV_trackId[i] = plist.Particle(p)->TrackId();
          fData->geant_tpcFV_status[i]  = plist.Particle(p)->StatusCode();
          //
          fData->geant_tpcFV_orig_E[i]             = plist.Particle(p)->E();
          fData->geant_tpcFV_orig_px[i]     = plist.Particle(p)->Px();
          fData->geant_tpcFV_orig_py[i]     = plist.Particle(p)->Py();
          fData->geant_tpcFV_orig_pz[i]     = plist.Particle(p)->Pz();
          fData->geant_tpcFV_orig_startx[i] = plist.Particle(p)->Vx();
          fData->geant_tpcFV_orig_starty[i] = plist.Particle(p)->Vy();
          fData->geant_tpcFV_orig_startz[i] = plist.Particle(p)->Vz();
          fData->geant_tpcFV_orig_startt[i] = plist.Particle(p)->T();
          fData->geant_tpcFV_orig_endx[i]   = plist.Particle(p)->EndX();
          fData->geant_tpcFV_orig_endy[i]   = plist.Particle(p)->EndY();
          fData->geant_tpcFV_orig_endz[i]   = plist.Particle(p)->EndZ();
          fData->geant_tpcFV_orig_endt[i]   = plist.Particle(p)->EndT();
          //
          fData->geant_tpcFV_startx[i]   = mcstart.X();
          fData->geant_tpcFV_starty[i]   = mcstart.Y();
          fData->geant_tpcFV_startz[i]   = mcstart.Z();
          fData->geant_tpcFV_endx[i]          = mcend.X();
          fData->geant_tpcFV_endy[i]          = mcend.Y();
          fData->geant_tpcFV_endz[i]          = mcend.Z();
          fData->geant_tpcFV_theta[i]          = plist.Particle(p)->Momentum().Theta();
          fData->geant_tpcFV_phi[i]          = plist.Particle(p)->Momentum().Phi();
          fData->geant_tpcFV_theta_xz[i] = mctheta_xz;
          fData->geant_tpcFV_theta_yz[i] = mctheta_yz;
          fData->geant_tpcFV_mom[i]          = plist.Particle(p)->Momentum().Vect().Mag();
          fData->geant_tpcFV_len[i]          = plen;
          i++;
        }
      }
      fData->geant_list_size_in_tpcFV = i;

      std::string pri("primary");
      int primary=0;
      int geant_particle=0;
      //determine the number of primary particles from geant:

      for( unsigned int i = 0; i < geant_part.size(); ++i ){
        geant_particle++;
        if(geant_part[i]->Process()==pri){
          primary++;
        }
      }

      fData->no_primaries=primary;
      fData->geant_list_size=geant_particle;
      //std::cout<<"\n"<<geant_list_size<<"\n";

      for( unsigned int i = 0; i < geant_part.size(); ++i ){

        if(geant_part[i]->Process()==pri){
          fData->process_primary[i]=1;
        }
        else{
          fData->process_primary[i]=0;
        }

        fData->Mother[i]=geant_part[i]->Mother();
        fData->TrackId[i]=geant_part[i]->TrackId();
        fData->pdg[i]=geant_part[i]->PdgCode();
        fData->Eng[i]=geant_part[i]->E();
        fData->Px[i]=geant_part[i]->Px();
        fData->Py[i]=geant_part[i]->Py();
        fData->Pz[i]=geant_part[i]->Pz();
        fData->StartPointx[i]=geant_part[i]->Vx();
        fData->StartPointy[i]=geant_part[i]->Vy();
        fData->StartPointz[i]=geant_part[i]->Vz();
        fData->EndPointx[i]=geant_part[i]->EndPosition()[0];
        fData->EndPointy[i]=geant_part[i]->EndPosition()[1];
        fData->EndPointz[i]=geant_part[i]->EndPosition()[2];
        fData->NumberDaughters[i]=geant_part[i]->NumberDaughters();
      }

    int currentMergedId = 1;
    int currentMotherId = 0;
    int currentMotherTrackId = 0;
    for(int i = 0; i < fData->geant_list_size; i++) {
       fData->MergedId[i] = 0;
    }

    for(int i = fData->geant_list_size - 1; i >= 0; i--) {
       if(fData->MergedId[i] == 0) {
          fData->MergedId[i] = currentMergedId;
          currentMotherId = fData->Mother[i];
          currentMotherTrackId = -1;
          while(currentMotherId > 0) {
             for(int j = 0; j < fData->geant_list_size; j++) {
                if(fData->TrackId[j] == currentMotherId) currentMotherTrackId = j;
             }
             if(fData->pdg[i] == fData->pdg[currentMotherTrackId]) {
                fData->MergedId[currentMotherTrackId] = currentMergedId;
                currentMotherId = fData->Mother[currentMotherTrackId];
             }
             else currentMotherId = 0;
          }
          currentMergedId++;
       }
    }
   }//if (mcevts_truth){//at least one mc record
  }//if (!isdata){
  fData->taulife = LArProp->ElectronLifetime();
  fTree->Fill();
  
  if (mf::MessageDrop::instance()->debugEnabled) {
    mf::LogDebug logStream("AnalysisTreeStructure");
    logStream
      << "Tree data structure has "
      << fData->data_size() << " bytes in data (" << fData->memory_size()
      << " allocated):"
      << "\n - " << fData->no_hits << " hits (" << fData->GetMaxHits() << ")"
      << "\n - " << fData->genie_no_primaries << " genie primaries (" << fData->GetMaxGeniePrimaries() << ")"
      << "\n - " << fData->no_primaries << " GEANT primaries (" << fData->GetMaxGEANTePrimaries() << ")"
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
    mf::LogDebug("AnalysisTreeStructure") << "Freeing the tree data structure";
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
    std::vector<cheat::TrackIDE> eveIDs = bt->HitToEveID(hit);

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
    //double momentum = track.MomentumAtPoint(i);
    //std::cout<<"\n"<<i<<"\t"<<momentum<<"\n";
    disp -= pos;
    result += disp.Mag();
    disp = pos;
  }
  return result;
}

// Length of MC particle, tracjectory by tracjectory.

double microboone::AnalysisTree::length(const simb::MCParticle& part, TVector3& start, TVector3& end)
{
  // Get geometry.
  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::DetectorProperties> detprop;

  // Get fiducial volume boundary.
  //double xmin = 0.;
  //double xmax = 2.*geom->DetHalfWidth();
  double ymin = -geom->DetHalfHeight();
  double ymax = geom->DetHalfHeight();
  double zmin = 0.;
  double zmax = geom->DetLength();

  const double fSamplingRate = 500;
  const double fReadOutWindowSize = 3200;
  int whatSpill=1;
  if (detprop->NumberTimeSamples()==3200)
           whatSpill = 0;
  else
         whatSpill = 1;

  double result = 0.;
  TVector3 disp;
  int n = part.NumberTrajectoryPoints();
  bool first = true;

  for(int i = 0; i < n; ++i) {
    try{
      // check if the particle is inside a TPC
      double mypos[3] = {part.Vx(i), part.Vy(i), part.Vz(i)};
      unsigned int tpc   = 0;
      unsigned int cstat = 0;
      geom->PositionToTPC(mypos, tpc, cstat);
    }
    catch(cet::exception &e){
      continue;
    }
    if(part.Vx(i) < (2.0*geom->DetHalfWidth()/fReadOutWindowSize)*(whatSpill*fReadOutWindowSize - part.T(i)*1./fSamplingRate ) ) continue;
    if(part.Vx(i) > (2.0*geom->DetHalfWidth()/fReadOutWindowSize)*((whatSpill+1) *fReadOutWindowSize - part.T(i)*1./fSamplingRate ) )
      continue;
    if(part.Vy(i) < ymin || part.Vy(i) > ymax) continue;
    if(part.Vz(i) < zmin || part.Vz(i) > zmax) continue;
    // Doing some manual shifting to account for
    // an interaction not occuring with the beam dump
    // we will reconstruct an x distance different from
    // where the particle actually passed to to the time
    // being different from in-spill interactions
    double newX = -(2.0*geom->DetHalfWidth()/fReadOutWindowSize)*(whatSpill*fReadOutWindowSize - part.T(i)*1./fSamplingRate ) + part.Vx(i);

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
  return result;
}

namespace microboone{

  DEFINE_ART_MODULE(AnalysisTree)

}

#endif

