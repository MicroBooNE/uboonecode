////////////////////////////////////////////////////////////////////////
// $Id: DBSCANfinderAna.cxx,v 1.36 2010/09/15  bpage Exp $
//
// module to create a TTree for analysis
//
// \author tjyang@fnal.gov, sowjanyag@phys.ksu.edu
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

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional>

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

namespace microboone {

  class AnalysisTreeDataStruct {
      public:
    AnalysisTreeDataStruct() { Clear(); }

    void Clear();

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
    Short_t  ntracks[kMaxTrackers];             //number of reconstructed tracks
    Float_t  trkke[kMaxTrackers][kMaxTrack][kNplanes];
    Float_t  trkrange[kMaxTrackers][kMaxTrack][kNplanes];
    Int_t    trkidtruth[kMaxTrackers][kMaxTrack][kNplanes]; //true geant trackid
    Int_t    trkpdgtruth[kMaxTrackers][kMaxTrack][kNplanes]; //true pdg code
    Float_t  trkefftruth[kMaxTrackers][kMaxTrack][kNplanes]; //completeness
    Float_t  trkpurtruth[kMaxTrackers][kMaxTrack][kNplanes]; //purity of track
    Float_t  trkpitchc[kMaxTrackers][kMaxTrack][kNplanes];
    Short_t  ntrkhits[kMaxTrackers][kMaxTrack][kNplanes];
    Float_t  trkdedx[kMaxTrackers][kMaxTrack][kNplanes][kMaxTrackHits];
    Float_t  trkdqdx[kMaxTrackers][kMaxTrack][kNplanes][kMaxTrackHits];
    Float_t  trkresrg[kMaxTrackers][kMaxTrack][kNplanes][kMaxTrackHits];
    Float_t  trkxyz[kMaxTrackers][kMaxTrack][kNplanes][kMaxTrackHits][3];

    // more track info
    Short_t   trkId[kMaxTrackers][kMaxTrack];
    Double_t  trkstartx[kMaxTrackers][kMaxTrack];      // starting x position.
    Double_t  trkstarty[kMaxTrackers][kMaxTrack];      // starting y position.
    Double_t  trkstartz[kMaxTrackers][kMaxTrack];      // starting z position.
    Double_t  trkstartd[kMaxTrackers][kMaxTrack];      // starting distance to boundary.
    Double_t  trkendx[kMaxTrackers][kMaxTrack];        // ending x position.
    Double_t  trkendy[kMaxTrackers][kMaxTrack];        // ending y position.
    Double_t  trkendz[kMaxTrackers][kMaxTrack];        // ending z position.
    Double_t  trkendd[kMaxTrackers][kMaxTrack];        // ending distance to boundary.
    Float_t   trktheta[kMaxTrackers][kMaxTrack];       // theta.
    Float_t   trkphi[kMaxTrackers][kMaxTrack];         // phi.
    Float_t   trkstartdcosx[kMaxTrackers][kMaxTrack];
    Float_t   trkstartdcosy[kMaxTrackers][kMaxTrack];
    Float_t   trkstartdcosz[kMaxTrackers][kMaxTrack];
    Float_t   trkenddcosx[kMaxTrackers][kMaxTrack];
    Float_t   trkenddcosy[kMaxTrackers][kMaxTrack];
    Float_t   trkenddcosz[kMaxTrackers][kMaxTrack];
    Float_t   trkthetaxz[kMaxTrackers][kMaxTrack];    // theta_xz.
    Float_t   trkthetayz[kMaxTrackers][kMaxTrack];    // theta_yz.
    Float_t   trkmom[kMaxTrackers][kMaxTrack];              // momentum.
    Float_t   trklen[kMaxTrackers][kMaxTrack];              // length.

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

    /// Connect this object with a tree
    void SetAddresses(TTree* pTree, const std::vector<std::string>& trackers);

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
          if (!pBranch)
            pTree->Branch(name.c_str(), address, leaflist.c_str() /*, bufsize */);
          else if (pBranch->GetAddress() != address) pBranch->SetAddress(address);
        } // operator()
      void operator()
        (std::string name, void* address, const std::stringstream& leaflist /*, int bufsize = 32000 */)
        { return this->operator() (name, address, leaflist.str()); }
      //@}
    }; // class BranchCreator

  }; // class AnalysisTreeDataStruct

  class AnalysisTree : public art::EDAnalyzer {

  public:

    explicit AnalysisTree(fhicl::ParameterSet const& pset);
    virtual ~AnalysisTree();

    /// read access to event
    void analyze(const art::Event& evt);
    void beginJob();
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

    // make sure the data structure exists
    void CreateData(bool bClearData = false)
      {
        if (!fData) fData = new AnalysisTreeDataStruct; // constructor Clear()s
        else if (bClearData) fData->Clear();
      } // CreateData()
    void SetAddresses()
      {
        if (!fData) {
          throw art::Exception(art::errors::LogicError)
            << "AnalysisTree::SetAddress(): no data";
        }
        if (!fTree) {
          throw art::Exception(art::errors::LogicError)
            << "AnalysisTree::SetAddress(): no tree";
        }
        fData->SetAddresses(fTree, fTrackModuleLabel);
      } // SetAddresses()
    void UpdateAddresses() { CreateData(); SetAddresses(); }
    
    /// Create the output tree and the data structures, if needed
    void CreateTree(bool bClearData = false);
  
  }; // class microboone::AnalysisTree
} // namespace microboone

namespace {
  class AutoResettingStringSteam: public std::ostringstream {
      public:
    AutoResettingStringSteam& operator() () { str(""); return *this; }
  }; // class AutoResettingStringSteam
} // local namespace


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
  for (int i = 0; i<kMaxHits; ++i){
    hit_plane[i] = -99;
    hit_wire[i] = -9999;
    hit_channel[i] = -9999;
    hit_peakT[i] = -99999;
    hit_charge[i] = -99999;
    hit_ph[i] = -99999;
  }

  for (int i = 0; i < kMaxTrackers; ++i){
    ntracks[i] = 0;
    std::fill(hit_trkid[i], hit_trkid[i] + kMaxHits, -9999);
    for (int j = 0; j < kMaxTrack; ++j){
      trkId[i][j]     = -9999;
      trkstartx[i][j] = -99999;
      trkstarty[i][j] = -99999;
      trkstartz[i][j] = -99999;
      trkstartd[i][j] = -99999;
      trkendx[i][j] = -99999;
      trkendy[i][j] = -99999;
      trkendz[i][j] = -99999;
      trkendd[i][j] = -99999;
      trktheta[i][j] = -99999;
      trkphi[i][j] = -99999;
      trkstartdcosx[i][j] = -99999;
      trkstartdcosy[i][j] = -99999;
      trkstartdcosz[i][j] = -99999;
      trkenddcosx[i][j] = -99999;
      trkenddcosy[i][j] = -99999;
      trkenddcosz[i][j] = -99999;
      trkthetaxz[i][j] = -99999;
      trkthetayz[i][j] = -99999;
      trkmom[i][j] = -99999;
      trklen[i][j] = -99999;
      for (int k = 0; k < kNplanes; ++k){
        trkke[i][j][k]       = -99999;
        trkrange[i][j][k]    = -99999;
        trkidtruth[i][j][k]  = -99999;
        trkpdgtruth[i][j][k] = -99999;
        trkefftruth[i][j][k] = -99999;
        trkpurtruth[i][j][k] = -99999;
        trkpitchc[i][j][k]   = -99999;
        ntrkhits[i][j][k]    = -9999;
        for(int l = 0; l < kMaxTrackHits; ++l) {
         trkdedx[i][j][k][l]  = 0;
         trkdqdx[i][j][k][l]  = 0;
         trkresrg[i][j][k][l] = 0;
         for (int m = 0; m<3; ++m)
           trkxyz[i][j][k][l][m] = 0;
        }
      }
    }
  }

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
    CreateBranch(BranchName, hit_trkid[i], BranchName + "[nohits]/S");

    BranchName = "ntracks_" + TrackLabel;
    CreateBranch(BranchName, &ntracks[i], BranchName + "/S");
    std::string NTracksIndexStr = "[" + BranchName + "]";
    
    BranchName = "trkId_" + TrackLabel;
    CreateBranch(BranchName, trkId[i], BranchName + NTracksIndexStr + "/S");
    
    BranchName = "trkke_" + TrackLabel;
    CreateBranch(BranchName, trkke[i], BranchName + NTracksIndexStr + "[3]/F");
    
    BranchName = "trkrange_" + TrackLabel;
    CreateBranch(BranchName, trkrange[i], BranchName + NTracksIndexStr + "[3]/F");
    
    BranchName = "trkidtruth_" + TrackLabel;
    CreateBranch(BranchName, trkidtruth[i], BranchName + NTracksIndexStr + "[3]/I");
    
    BranchName = "trkpdgtruth_" + TrackLabel;
    CreateBranch(BranchName, trkpdgtruth[i], BranchName + NTracksIndexStr + "[3]/I");
    
    BranchName = "trkefftruth_" + TrackLabel;
    CreateBranch(BranchName, trkefftruth[i], BranchName + NTracksIndexStr + "[3]/F");
    
    BranchName = "trkpurtruth_" + TrackLabel;
    CreateBranch(BranchName, trkpurtruth[i], BranchName + NTracksIndexStr + "[3]/F");
    
    BranchName = "trkpitchc_" + TrackLabel;
    CreateBranch(BranchName, trkpitchc[i], BranchName + NTracksIndexStr + "[3]/F");
    
    BranchName = "ntrkhits_" + TrackLabel;
    CreateBranch(BranchName, ntrkhits[i], BranchName + NTracksIndexStr + "[3]/S");
    
    // trkdedx_<TrackerName>[ntracks_<TrackerName>][3][<kMaxTrackHits>]/F"
    BranchName = "trkdedx_" + TrackLabel;
    CreateBranch(BranchName, trkdedx[i], BranchName + NTracksIndexStr + "[3]" + MaxTrackHitsIndexStr + "/F");
    
    // trkdqdx_<TrackerName>[ntracks_<TrackerName>][3][<kMaxTrackHits>]/F"
    BranchName = "trkdqdx_" + TrackLabel;
    CreateBranch(BranchName, trkdqdx[i], BranchName + NTracksIndexStr + "[3]" + MaxTrackHitsIndexStr + "/F");
    
    // trkresrg_<TrackerName>[ntracks_<TrackerName>][3][<kMaxTrackHits>]/F"
    BranchName = "trkresrg_" + TrackLabel;
    CreateBranch(BranchName, trkresrg[i], BranchName + NTracksIndexStr + "[3]" + MaxTrackHitsIndexStr + "/F");
    
    // trkresrg_<TrackerName>[ntracks_<TrackerName>][3][<kMaxTrackHits>]/F"
    BranchName = "trkxyz_" + TrackLabel;
    CreateBranch(BranchName, trkxyz[i], BranchName + NTracksIndexStr + "[3]" + MaxTrackHitsIndexStr + "/F");
    
    BranchName = "trkstartx_" + TrackLabel;
    CreateBranch(BranchName, trkstartx[i], BranchName + NTracksIndexStr + "/D");
    
    BranchName = "trkstarty_" + TrackLabel;
    CreateBranch(BranchName, trkstarty[i], BranchName + NTracksIndexStr + "/D");
    
    BranchName = "trkstartz_" + TrackLabel;
    CreateBranch(BranchName, trkstartz[i], BranchName + NTracksIndexStr + "/D");
    
    BranchName = "trkstartd_" + TrackLabel;
    CreateBranch(BranchName, trkstartd[i], BranchName + NTracksIndexStr + "/D");
    
    BranchName = "trkendx_" + TrackLabel;
    CreateBranch(BranchName, trkendx[i], BranchName + NTracksIndexStr + "/D");
    
    BranchName = "trkendy_" + TrackLabel;
    CreateBranch(BranchName, trkendy[i], BranchName + NTracksIndexStr + "/D");
    
    BranchName = "trkendz_" + TrackLabel;
    CreateBranch(BranchName, trkendz[i], BranchName + NTracksIndexStr + "/D");
    
    BranchName = "trkendd_" + TrackLabel;
    CreateBranch(BranchName, trkendd[i], BranchName + NTracksIndexStr + "/D");
    
    BranchName = "trktheta_" + TrackLabel;
    CreateBranch(BranchName, trktheta[i], BranchName + NTracksIndexStr + "/F");
    
    BranchName = "trkphi_" + TrackLabel;
    CreateBranch(BranchName, trkphi[i], BranchName + NTracksIndexStr + "/F");
    
    BranchName = "trkstartdcosx_" + TrackLabel;
    CreateBranch(BranchName, trkstartdcosx[i], BranchName + NTracksIndexStr + "/F");
    
    BranchName = "trkstartdcosy_" + TrackLabel;
    CreateBranch(BranchName, trkstartdcosy[i], BranchName + NTracksIndexStr + "/F");
    
    BranchName = "trkstartdcosz_" + TrackLabel;
    CreateBranch(BranchName, trkstartdcosz[i], BranchName + NTracksIndexStr + "/F");
    
    BranchName = "trkenddcosx_" + TrackLabel;
    CreateBranch(BranchName, trkenddcosx[i], BranchName + NTracksIndexStr + "/F");
    
    BranchName = "trkenddcosy_" + TrackLabel;
    CreateBranch(BranchName, trkenddcosy[i], BranchName + NTracksIndexStr + "/F");
    
    BranchName = "trkenddcosz_" + TrackLabel;
    CreateBranch(BranchName, trkenddcosz[i], BranchName + NTracksIndexStr + "/F");
    
    BranchName = "trkthetaxz_" + TrackLabel;
    CreateBranch(BranchName, trkthetaxz[i], BranchName + NTracksIndexStr + "/F");
    
    BranchName = "trkthetayz_" + TrackLabel;
    CreateBranch(BranchName, trkthetayz[i], BranchName + NTracksIndexStr + "/F");
    
    BranchName = "trkmom_" + TrackLabel;
    CreateBranch(BranchName, trkmom[i], BranchName + NTracksIndexStr + "/F");
    
    BranchName = "trklen_" + TrackLabel;
    CreateBranch(BranchName, trklen[i], BranchName + NTracksIndexStr + "/F");
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
  fParticleIDModuleLabel    (pset.get< std::string >("ParticleIDModuleLabel")   )
{
}

//-------------------------------------------------
microboone::AnalysisTree::~AnalysisTree()
{
  delete fData;
  fData = nullptr;
}

void microboone::AnalysisTree::CreateTree(bool bClearData /* = false */) {
  if (!fTree) {
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("anatree","analysis tree");
  }
  CreateData(bClearData);
  SetAddresses();
} // microboone::AnalysisTree::CreateTree()


void microboone::AnalysisTree::beginJob(){
  CreateTree(); // TODO move this to the event (optionally)
}

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
  // make sure there is the data, the tree and everything
  CreateTree(true);

  /// transfer the run and subrun data to the tree data object
//  fData->RunData = RunData;
  fData->SubRunData = SubRunData;

  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);

  std::vector< art::Handle< std::vector<recob::Track> > > trackListHandle(fData->kNTracker);
  std::vector< std::vector<art::Ptr<recob::Track> > > tracklist(fData->kNTracker);
  for (unsigned int it = 0; it < fTrackModuleLabel.size(); ++it){
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
  fData->no_hits=hitlist.size();
  for (size_t i = 0; i<hitlist.size(); ++i){//loop over hits
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
  for (unsigned int it1=0; it1<fTrackModuleLabel.size(); ++it1){
    fData->ntracks[it1]=tracklist[it1].size();
    for(unsigned int i=0; i<tracklist[it1].size();++i){//loop over tracks

      art::Ptr<recob::Track> ptrack(trackListHandle[it1], i);
      const recob::Track& track = *ptrack;
      
      TVector3 pos, dir_start, dir_end, end;
      double tlen = 0., mom = 0.;
      
      //we need to use Bezier methods for Bezier tracks
      if(fTrackModuleLabel[it1].compare("beziertracker")==0){
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
            fData->trkId[it1][i]         = i;  //bezier has some screwed up track IDs
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
            fData->trkId[it1][i]                    = track.ID();
        }
      }
      double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
      double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());
      double dpos = bdist(pos);
      double dend = bdist(end);
      
      fData->trkstartx[it1][i]             = pos.X();
      fData->trkstarty[it1][i]             = pos.Y();
      fData->trkstartz[it1][i]             = pos.Z();
      fData->trkstartd[it1][i]             = dpos;
      fData->trkendx[it1][i]               = end.X();
      fData->trkendy[it1][i]               = end.Y();
      fData->trkendz[it1][i]               = end.Z();
      fData->trkendd[it1][i]               = dend;
      fData->trktheta[it1][i]              = dir_start.Theta();
      fData->trkphi[it1][i]                = dir_start.Phi();
      fData->trkstartdcosx[it1][i]         = dir_start.X();
      fData->trkstartdcosy[it1][i]         = dir_start.Y();
      fData->trkstartdcosz[it1][i]         = dir_start.Z();
      fData->trkenddcosx[it1][i]           = dir_end.X();
      fData->trkenddcosy[it1][i]           = dir_end.Y();
      fData->trkenddcosz[it1][i]           = dir_end.Z();
      fData->trkthetaxz[it1][i]            = theta_xz;
      fData->trkthetayz[it1][i]            = theta_yz;
      fData->trkmom[it1][i]                = mom;
      fData->trklen[it1][i]                = tlen;

      art::FindMany<anab::Calorimetry> fmcal(trackListHandle[it1], evt, fCalorimetryModuleLabel[it1]);
      if (fmcal.isValid()){
        std::vector<const anab::Calorimetry*> calos = fmcal.at(i);
        //std::cout<<"calo size "<<calos.size()<<std::endl;
        for (size_t j = 0; j<calos.size(); ++j){
          fData->trkke[it1][i][j]    = calos[j]->KineticEnergy();
          fData->trkrange[it1][i][j] = calos[j]->Range();
          fData->trkpitchc[it1][i][j]= calos[j] -> TrkPitchC();
          fData->ntrkhits[it1][i][j] = calos[j] -> dEdx().size();
          for(int k = 0; k < fData->ntrkhits[it1][i][j]; ++k) {
            fData->trkdedx[it1][i][j][k]  = (calos[j] -> dEdx())[k];
            fData->trkdqdx[it1][i][j][k]  = (calos[j] -> dQdx())[k];
            fData->trkresrg[it1][i][j][k] = (calos[j] -> ResidualRange())[k];
            fData->trkxyz[it1][i][j][k][0] = (calos[j] -> XYZ())[k].X();
            fData->trkxyz[it1][i][j][k][1] = (calos[j] -> XYZ())[k].Y();
            fData->trkxyz[it1][i][j][k][2] = (calos[j] -> XYZ())[k].Z();
          } // for track hits
        } // for calorimetry info
      } // if has calorimetry info

      //track truth information
      if (!fData->isdata){
        //get the hits on each plane
        art::FindManyP<recob::Hit>      fmht(trackListHandle[it1], evt, fTrackModuleLabel[it1]);
        std::vector< art::Ptr<recob::Hit> > allHits = fmht.at(i);
        std::vector< art::Ptr<recob::Hit> > hits[kNplanes];

        for(size_t ah = 0; ah < allHits.size(); ++ah){
          if (/* allHits[ah]->WireID().Plane >= 0 && */ // always true
            allHits[ah]->WireID().Plane <  3){
            hits[allHits[ah]->WireID().Plane].push_back(allHits[ah]);
          }
        }
        
        for (size_t ipl = 0; ipl < 3; ++ipl){
          double maxe = 0;
          HitsPurity(hits[ipl],fData->trkidtruth[it1][i][ipl],fData->trkpurtruth[it1][i][ipl],maxe);
        //std::cout<<"\n"<<it1<<"\t"<<i<<"\t"<<ipl<<"\t"<<trkidtruth[it1][i][ipl]<<"\t"<<trkpurtruth[it1][i][ipl]<<"\t"<<maxe;
          if (fData->trkidtruth[it1][i][ipl]>0){
            const simb::MCParticle *particle = bt->TrackIDToParticle(fData->trkidtruth[it1][i][ipl]);
            const std::vector<sim::IDE> vide = bt->TrackIDToSimIDE(fData->trkidtruth[it1][i][ipl]);
            double tote = 0;
            for (size_t iide = 0; iide<vide.size(); ++iide){
              tote += vide[iide].energy;
            }
            fData->trkpdgtruth[it1][i][ipl] = particle->PdgCode();
            fData->trkefftruth[it1][i][ipl] = maxe/(tote/kNplanes); //tote include both induction and collection energies
          //std::cout<<"\n"<<trkpdgtruth[it1][i][ipl]<<"\t"<<trkefftruth[it1][i][ipl];
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

