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

const int kNplanes       = 3;     //number of wire planes
const int kMaxTrack      = 1000;  //maximum number of tracks
const int kMaxClusters   = 1000;  //maximum number of clusters
const int kMaxHits       = 20000; //maximum number of hits;
const int kMaxPrimaries  = 20000;  //maximum number of primary particles
const int kMaxTrackHits  = 1000;  //maximum number of hits on a track
const int kMaxTrackers   = 10;    //number of trackers passed into fTrackModuleLabel 

namespace microboone {
   
  class AnalysisTree : public art::EDAnalyzer {

  public:
          
    explicit AnalysisTree(fhicl::ParameterSet const& pset); 
    virtual ~AnalysisTree();
 
    /// read access to event
    void analyze(const art::Event& evt);
    void beginJob();
    void beginSubRun(const art::SubRun& sr);

  private:
    
    void   HitsPurity(std::vector< art::Ptr<recob::Hit> > const& hits, int& trackid, double& purity, double& maxe);
    void   ResetVars();
    double length(const recob::Track& track);
    double length(const simb::MCParticle& part, TVector3& start, TVector3& end);
    double bdist(const TVector3& pos);

    TTree* fTree;
    //run information
    int    run;                  //run number
    int    subrun;               //subrun number
    int    event;                //event number
    double evttime;              //event time in sec
    double beamtime;             //beam time
    double pot;                  //protons on target
    double taulife;              //electron lifetime
    int    isdata;               //flag, 0=MC 1=data

    //hit information
    int no_hits;                     //number of hits
    int hit_plane[kMaxHits];         //plane number
    int hit_wire[kMaxHits];          //wire number
    int hit_channel[kMaxHits];       //channel ID
    double hit_peakT[kMaxHits];      //peak time
    double hit_charge[kMaxHits];     //charge (area)
    double hit_ph[kMaxHits];         //amplitude
    int    hit_trkid[kMaxTrackers][kMaxHits];      //is this hit associated with a reco track?

    //track information
    int    kNTracker;
    int    ntracks[kMaxTrackers];             //number of reconstructed tracks
    double trkke[kMaxTrackers][kMaxTrack][kNplanes];
    double trkrange[kMaxTrackers][kMaxTrack][kNplanes];
    int    trkidtruth[kMaxTrackers][kMaxTrack][kNplanes]; //true geant trackid
    int    trkpdgtruth[kMaxTrackers][kMaxTrack][kNplanes]; //true pdg code
    double trkefftruth[kMaxTrackers][kMaxTrack][kNplanes]; //completeness
    double trkpurtruth[kMaxTrackers][kMaxTrack][kNplanes]; //purity of track    
    double trkpitchc[kMaxTrackers][kMaxTrack][kNplanes];
    int    ntrkhits[kMaxTrackers][kMaxTrack][kNplanes];
    double trkdedx[kMaxTrackers][kMaxTrack][kNplanes][kMaxTrackHits];
    double trkdqdx[kMaxTrackers][kMaxTrack][kNplanes][kMaxTrackHits];
    double trkresrg[kMaxTrackers][kMaxTrack][kNplanes][kMaxTrackHits];
    double trkxyz[kMaxTrackers][kMaxTrack][kNplanes][kMaxTrackHits][3];
    
    // more track info
    char   S_temp[50],S1_temp[50], S2_temp[50], Str_temp[50];
    int    trkId[kMaxTrackers][kMaxTrack];
    double trkstartx[kMaxTrackers][kMaxTrack];      // starting x position.
    double trkstarty[kMaxTrackers][kMaxTrack];      // starting y position.
    double trkstartz[kMaxTrackers][kMaxTrack];      // starting z position.
    double trkstartd[kMaxTrackers][kMaxTrack];      // starting distance to boundary.
    double trkendx[kMaxTrackers][kMaxTrack];        // ending x position.
    double trkendy[kMaxTrackers][kMaxTrack];        // ending y position.
    double trkendz[kMaxTrackers][kMaxTrack];        // ending z position.
    double trkendd[kMaxTrackers][kMaxTrack];        // ending distance to boundary.
    double trktheta[kMaxTrackers][kMaxTrack];       // theta.
    double trkphi[kMaxTrackers][kMaxTrack];	      // phi.
    double trkstartdcosx[kMaxTrackers][kMaxTrack]; 
    double trkstartdcosy[kMaxTrackers][kMaxTrack]; 
    double trkstartdcosz[kMaxTrackers][kMaxTrack]; 
    double trkenddcosx[kMaxTrackers][kMaxTrack];
    double trkenddcosy[kMaxTrackers][kMaxTrack];
    double trkenddcosz[kMaxTrackers][kMaxTrack];
    double trkthetaxz[kMaxTrackers][kMaxTrack];    // theta_xz.
    double trkthetayz[kMaxTrackers][kMaxTrack];    // theta_yz.
    double trkmom[kMaxTrackers][kMaxTrack];	      // momentum.
    double trklen[kMaxTrackers][kMaxTrack];	      // length.
       
    //mctruth information
    int    mcevts_truth;    //number of neutrino interactions in the spill
    int    nuPDG_truth;     //neutrino PDG code
    int    ccnc_truth;      //0=CC 1=NC
    int    mode_truth;      //0=QE/El, 1=RES, 2=DIS, 3=Coherent production
    double enu_truth;       //true neutrino energy
    double Q2_truth;        //Momentum transfer squared
    double W_truth;         //hadronic invariant mass
    int    hitnuc_truth;    //hit nucleon
    double nuvtxx_truth;    //neutrino vertex x
    double nuvtxy_truth;    //neutrino vertex y
    double nuvtxz_truth;    //neutrino vertex z
    double nu_dcosx_truth;  //neutrino dcos x
    double nu_dcosy_truth;  //neutrino dcos y
    double nu_dcosz_truth;  //neutrino dcos z
    double lep_mom_truth;   //lepton momentum
    double lep_dcosx_truth; //lepton dcos x
    double lep_dcosy_truth; //lepton dcos y
    double lep_dcosz_truth; //lepton dcos z

    //flux information
    double tpx_flux;        //Px of parent particle leaving BNB target
    double tpy_flux;        //Py of parent particle leaving BNB target
    double tpz_flux;        //Pz of parent particle leaving BNB target
    int    tptype_flux;     //Type of parent particle leaving BNB target

    //genie information
    int    genie_no_primaries;
    int    genie_primaries_pdg[kMaxPrimaries];
    double genie_Eng[kMaxPrimaries];
    double genie_Px[kMaxPrimaries];
    double genie_Py[kMaxPrimaries];
    double genie_Pz[kMaxPrimaries];
    double genie_P[kMaxPrimaries];
    int    genie_status_code[kMaxPrimaries];
    double genie_mass[kMaxPrimaries];
    int    genie_trackID[kMaxPrimaries];
    int    genie_ND[kMaxPrimaries];
    int    genie_mother[kMaxPrimaries];
 
    //geant information
    int    no_primaries;      //number of primary geant particles
    int    geant_list_size;  //number of all geant particles
    int    geant_list_size_in_tpcFV;
    int    pdg[kMaxPrimaries];
    double Eng[kMaxPrimaries];
    double Px[kMaxPrimaries];
    double Py[kMaxPrimaries];
    double Pz[kMaxPrimaries];
    double StartPointx[kMaxPrimaries];
    double StartPointy[kMaxPrimaries];
    double StartPointz[kMaxPrimaries];
    double EndPointx[kMaxPrimaries];
    double EndPointy[kMaxPrimaries];
    double EndPointz[kMaxPrimaries];
    int    NumberDaughters[kMaxPrimaries];
    int    TrackId[kMaxPrimaries];
    int    Mother[kMaxPrimaries];
    int    process_primary[kMaxPrimaries];
    int    MergedId[kMaxPrimaries]; //geant track segments, which belong to the same particle, get the same
    
    // more geant information
    int   geant_tpcFV_status[kMaxPrimaries];
    int   geant_tpcFV_trackId[kMaxPrimaries];
    int   geant_tpcFV_pdg[kMaxPrimaries];
    
    double geant_tpcFV_orig_E[kMaxPrimaries];
    double geant_tpcFV_orig_px[kMaxPrimaries];
    double geant_tpcFV_orig_py[kMaxPrimaries];
    double geant_tpcFV_orig_pz[kMaxPrimaries];
    double geant_tpcFV_orig_startx[kMaxPrimaries];
    double geant_tpcFV_orig_starty[kMaxPrimaries];
    double geant_tpcFV_orig_startz[kMaxPrimaries];
    double geant_tpcFV_orig_startt[kMaxPrimaries];
    double geant_tpcFV_orig_endx[kMaxPrimaries];
    double geant_tpcFV_orig_endy[kMaxPrimaries];
    double geant_tpcFV_orig_endz[kMaxPrimaries];
    double geant_tpcFV_orig_endt[kMaxPrimaries];
    
    double geant_tpcFV_startx[kMaxPrimaries];	  // starting x position.
    double geant_tpcFV_starty[kMaxPrimaries];	  // starting y position.
    double geant_tpcFV_startz[kMaxPrimaries];	  // starting z position.
    double geant_tpcFV_startd[kMaxPrimaries];	  // starting distance to boundary.
    double geant_tpcFV_endx[kMaxPrimaries];	  // ending x position.
    double geant_tpcFV_endy[kMaxPrimaries];	  // ending y position.
    double geant_tpcFV_endz[kMaxPrimaries];	  // ending z position.
    double geant_tpcFV_endd[kMaxPrimaries];	  // ending distance to boundary.
    double geant_tpcFV_theta[kMaxPrimaries];	  // theta.
    double geant_tpcFV_phi[kMaxPrimaries];	  // phi.
    double geant_tpcFV_theta_xz[kMaxPrimaries];    // theta_xz.
    double geant_tpcFV_theta_yz[kMaxPrimaries];    // theta_yz.
    double geant_tpcFV_mom[kMaxPrimaries];         // momentum.
    double geant_tpcFV_len[kMaxPrimaries];         // length.


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

  };
}

//-------------------------------------------------

microboone::AnalysisTree::AnalysisTree(fhicl::ParameterSet const& pset) : 
  EDAnalyzer(pset),
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
}

void microboone::AnalysisTree::beginJob(){

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("anatree","analysis tree");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("evttime",&evttime,"evttime/D");
  fTree->Branch("beamtime",&beamtime,"beamtime/D");
  fTree->Branch("pot",&pot,"pot/D");
  fTree->Branch("isdata",&isdata,"isdata/I");
  fTree->Branch("taulife",&taulife,"taulife/D");

  fTree->Branch("no_hits",&no_hits,"no_hits/I");
  fTree->Branch("hit_plane",hit_plane,"hit_plane[no_hits]/I");
  fTree->Branch("hit_wire",hit_wire,"hit_wire[no_hits]/I");
  fTree->Branch("hit_channel",hit_channel,"hit_channel[no_hits]/I");
  fTree->Branch("hit_peakT",hit_peakT,"hit_peakT[no_hits]/D");
  fTree->Branch("hit_charge",hit_charge,"hit_charge[no_hits]/D");
  fTree->Branch("hit_ph",hit_ph,"hit_ph[no_hits]/D");

  kNTracker = fTrackModuleLabel.size();
  fTree->Branch("kNTracker",&kNTracker,"kNTracker/I");
  for(int i=0; i<kNTracker; i++){
    sprintf(Str_temp,"hit_trkid_%s",fTrackModuleLabel[i].c_str());
    sprintf(S1_temp,"hit_trkid_%s[no_hits]/I",fTrackModuleLabel[i].c_str()); 	
    fTree->Branch(Str_temp,hit_trkid[i],S1_temp);
    sprintf(S_temp,"ntracks_%s",fTrackModuleLabel[i].c_str());
    sprintf(S1_temp,"ntracks_%s/I",fTrackModuleLabel[i].c_str());
    fTree->Branch(S_temp, &ntracks[i], S1_temp);
    sprintf(S1_temp,"trkId_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkId_%s[%s]/I",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkId[i], S2_temp);  
    sprintf(S1_temp,"trkke_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkke_%s[%s][3]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkke[i], S2_temp);
    sprintf(S1_temp,"trkrange_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkrange_%s[%s][3]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkrange[i], S2_temp);
    sprintf(S1_temp,"trkidtruth_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkidtruth_%s[%s][3]/I",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkidtruth[i], S2_temp);
    sprintf(S1_temp,"trkpdgtruth_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkpdgtruth_%s[%s][3]/I",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkpdgtruth[i], S2_temp);
    sprintf(S1_temp,"trkefftruth_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkefftruth_%s[%s][3]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkefftruth[i], S2_temp);
    sprintf(S1_temp,"trkpurtruth_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkpurtruth_%s[%s][3]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkpurtruth[i], S2_temp);
    sprintf(S1_temp,"trkpitchc_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkpitchc_%s[%s][3]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp,trkpitchc[i],S2_temp);
    sprintf(S1_temp,"ntrkhits_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"ntrkhits_%s[%s][3]/I",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, ntrkhits[i], S2_temp);    
    sprintf(S1_temp,"trkdedx_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkdedx_%s[%s][3][1000]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp,trkdedx[i],S2_temp);    
    sprintf(S1_temp,"trkdqdx_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkdqdx_%s[%s][3][1000]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp,trkdqdx[i],S2_temp);    
    sprintf(S1_temp,"trkresrg_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkresrg_%s[%s][3][1000]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp,trkresrg[i],S2_temp);
    sprintf(S1_temp,"trkxyz_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkxyz_%s[%s][3][1000][3]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp,trkxyz[i],S2_temp);
    sprintf(S1_temp,"trkstartx_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkstartx_%s[%s]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkstartx[i], S2_temp);
    sprintf(S1_temp,"trkstarty_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkstarty_%s[%s]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkstarty[i], S2_temp);
    sprintf(S1_temp,"trkstartz_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkstartz_%s[%s]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkstartz[i], S2_temp);
    sprintf(S1_temp,"trkstartd_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkstartd_%s[%s]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkstartd[i], S2_temp);
    sprintf(S1_temp,"trkendx_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkendx_%s[%s]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkendx[i], S2_temp);
    sprintf(S1_temp,"trkendy_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkendy_%s[%s]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkendy[i], S2_temp);
    sprintf(S1_temp,"trkendz_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkendz_%s[%s]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkendz[i], S2_temp);
    sprintf(S1_temp,"trkendd_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkendd_%s[%s]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkendd[i], S2_temp);
    sprintf(S1_temp,"trktheta_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trktheta_%s[%s]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trktheta[i], S2_temp);
    sprintf(S1_temp,"trkphi_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkphi_%s[%s]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkphi[i], S2_temp);
    sprintf(S1_temp,"trkstartdcosx_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkstartdcosx_%s[%s]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkstartdcosx[i], S2_temp);
    sprintf(S1_temp,"trkstartdcosy_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkstartdcosy_%s[%s]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkstartdcosy[i], S2_temp);
    sprintf(S1_temp,"trkstartdcosz_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkstartdcosz_%s[%s]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkstartdcosz[i], S2_temp);
    sprintf(S1_temp,"trkenddcosx_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trksenddcosx_%s[%s]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkenddcosx[i], S2_temp);
    sprintf(S1_temp,"trkenddcosy_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trksenddcosy_%s[%s]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkenddcosy[i], S2_temp);
    sprintf(S1_temp,"trkenddcosz_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkenddcosz_%s[%s]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkenddcosz[i], S2_temp);
    sprintf(S1_temp,"trkthetaxz_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkthetaxz_%s[%s]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkthetaxz[i], S2_temp);
    sprintf(S1_temp,"trkthetayz_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkthetayz_%s[%s]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkthetayz[i], S2_temp);
    sprintf(S1_temp,"trkmom_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trkmom_%s[%s]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trkmom[i], S2_temp);
    sprintf(S1_temp,"trklen_%s",fTrackModuleLabel[i].c_str());
    sprintf(S2_temp,"trklen_%s[%s]/D",fTrackModuleLabel[i].c_str(),S_temp);
    fTree->Branch(S1_temp, trklen[i], S2_temp);
  }
  
  fTree->Branch("mcevts_truth",&mcevts_truth,"mcevts_truth/I");
  fTree->Branch("nuPDG_truth",&nuPDG_truth,"nuPDG_truth/I");
  fTree->Branch("ccnc_truth",&ccnc_truth,"ccnc_truth/I");
  fTree->Branch("mode_truth",&mode_truth,"mode_truth/I");
  fTree->Branch("enu_truth",&enu_truth,"enu_truth/D");
  fTree->Branch("Q2_truth",&Q2_truth,"Q2_truth/D");
  fTree->Branch("W_truth",&W_truth,"W_truth/D");
  fTree->Branch("hitnuc_truth",&hitnuc_truth,"hitnuc_truth/I");
  fTree->Branch("nuvtxx_truth",&nuvtxx_truth,"nuvtxx_truth/D");
  fTree->Branch("nuvtxy_truth",&nuvtxy_truth,"nuvtxy_truth/D");
  fTree->Branch("nuvtxz_truth",&nuvtxz_truth,"nuvtxz_truth/D");
  fTree->Branch("nu_dcosx_truth",&nu_dcosx_truth,"nu_dcosx_truth/D");
  fTree->Branch("nu_dcosy_truth",&nu_dcosy_truth,"nu_dcosy_truth/D");
  fTree->Branch("nu_dcosz_truth",&nu_dcosz_truth,"nu_dcosz_truth/D");
  fTree->Branch("lep_mom_truth",&lep_mom_truth,"lep_mom_truth/D");
  fTree->Branch("lep_dcosx_truth",&lep_dcosx_truth,"lep_dcosx_truth/D");
  fTree->Branch("lep_dcosy_truth",&lep_dcosy_truth,"lep_dcosy_truth/D");
  fTree->Branch("lep_dcosz_truth",&lep_dcosz_truth,"lep_dcosz_truth/D");

  fTree->Branch("tpx_flux",&tpx_flux,"tpx_flux/D");
  fTree->Branch("tpy_flux",&tpy_flux,"tpy_flux/D");
  fTree->Branch("tpz_flux",&tpz_flux,"tpz_flux/D");
  fTree->Branch("tptype_flux",&tptype_flux,"tptype_flux/I");

  fTree->Branch("genie_no_primaries",&genie_no_primaries,"genie_no_primaries/I");
  fTree->Branch("genie_primaries_pdg",genie_primaries_pdg,"genie_primaries_pdg[genie_no_primaries]/I");
  fTree->Branch("genie_Eng",genie_Eng,"genie_Eng[genie_no_primaries]/D");
  fTree->Branch("genie_Px",genie_Px,"genie_Px[genie_no_primaries]/D");
  fTree->Branch("genie_Py",genie_Py,"genie_Py[genie_no_primaries]/D");
  fTree->Branch("genie_Pz",genie_Pz,"genie_Pz[genie_no_primaries]/D");
  fTree->Branch("genie_P",genie_P,"genie_P[genie_no_primaries]/D");
  fTree->Branch("genie_status_code",genie_status_code,"genie_status_code[genie_no_primaries]/I");
  fTree->Branch("genie_mass",genie_mass,"genie_mass[genie_no_primaries]/D");
  fTree->Branch("genie_trackID",genie_trackID,"genie_trackID[genie_no_primaries]/I");
  fTree->Branch("genie_ND",genie_ND,"genie_ND[genie_no_primaries]/I");
  fTree->Branch("genie_mother",genie_mother,"genie_mother[genie_no_primaries]/I");

  fTree->Branch("no_primaries",&no_primaries,"no_primaries/I");
  fTree->Branch("geant_list_size",&geant_list_size,"geant_list_size/I");
  
  fTree->Branch("pdg",pdg,"pdg[geant_list_size]/I");
  fTree->Branch("Eng",Eng,"Eng[geant_list_size]/D");
  fTree->Branch("Px",Px,"Px[geant_list_size]/D");
  fTree->Branch("Py",Py,"Py[geant_list_size]/D");
  fTree->Branch("Pz",Pz,"Pz[geant_list_size]/D");
  fTree->Branch("StartPointx",StartPointx,"StartPointx[geant_list_size]/D");
  fTree->Branch("StartPointy",StartPointy,"StartPointy[geant_list_size]/D");
  fTree->Branch("StartPointz",StartPointz,"StartPointz[geant_list_size]/D");
  fTree->Branch("EndPointx",EndPointx,"EndPointx[geant_list_size]/D");
  fTree->Branch("EndPointy",EndPointy,"EndPointy[geant_list_size]/D");
  fTree->Branch("EndPointz",EndPointz,"EndPointz[geant_list_size]/D");
  fTree->Branch("NumberDaughters",NumberDaughters,"NumberDaughters[geant_list_size]/I");
  fTree->Branch("Mother",Mother,"Mother[geant_list_size]/I");
  fTree->Branch("TrackId",TrackId,"TrackId[geant_list_size]/I");
  fTree->Branch("MergedId", MergedId, "MergedId[geant_list_size]/I");
  fTree->Branch("process_primary",process_primary,"process_primary[geant_list_size]/I");  
  
  fTree->Branch("geant_list_size_in_tpcFV",&geant_list_size_in_tpcFV,"geant_list_size_in_tpcFV/I");  
  fTree->Branch("geant_tpcFV_pdg", geant_tpcFV_pdg, "geant_tpcFV_pdg[geant_list_size_in_tpcFV]/I");
  fTree->Branch("geant_tpcFV_status", geant_tpcFV_status, "geant_tpcFV_status[geant_list_size_in_tpcFV]/I");
  fTree->Branch("geant_tpcFV_trackId", geant_tpcFV_trackId, "geant_tpcFV_trackId[geant_list_size_in_tpcFV]/I");
  fTree->Branch("geant_tpcFV_orig_E", geant_tpcFV_orig_E, "geant_tpcFV_orig_E[geant_list_size_in_tpcFV]/F");
  fTree->Branch("geant_tpcFV_orig_px", geant_tpcFV_orig_px, "geant_tpcFV_orig_px[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_orig_py", geant_tpcFV_orig_py, "geant_tpcFV_orig_py[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_orig_pz", geant_tpcFV_orig_pz, "geant_tpcFV_orig_pz[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_orig_startx", geant_tpcFV_orig_startx, "geant_tpcFV_orig_startx[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_orig_starty", geant_tpcFV_orig_starty, "geant_tpcFV_orig_starty[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_orig_startz", geant_tpcFV_orig_startz, "geant_tpcFV_orig_startz[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_orig_startt", geant_tpcFV_orig_startt, "geant_tpcFV_orig_startt[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_orig_endx", geant_tpcFV_orig_endx, "geant_tpcFV_orig_endx[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_orig_endy", geant_tpcFV_orig_endy, "geant_tpcFV_orig_endy[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_orig_endz", geant_tpcFV_orig_endz, "geant_tpcFV_orig_endz[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_orig_endt", geant_tpcFV_orig_endt, "geant_tpcFV_orig_endt[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_startx", geant_tpcFV_startx, "geant_tpcFV_startx[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_starty", geant_tpcFV_starty, "geant_tpcFV_starty[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_startz", geant_tpcFV_startz, "geant_tpcFV_startz[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_startd", geant_tpcFV_startd, "geant_tpcFV_startd[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_endx", geant_tpcFV_endx, "geant_tpcFV_endx[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_endy", geant_tpcFV_endy, "geant_tpcFV_endy[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_endz", geant_tpcFV_endz, "geant_tpcFV_endz[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_endd", geant_tpcFV_endd, "geant_tpcFV_endd[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_theta", geant_tpcFV_theta, "geant_tpcFV_theta[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_phi", geant_tpcFV_phi, "geant_tpcFV_phi[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_theta_xz", geant_tpcFV_theta_xz, "geant_tpcFV_theta_xz[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_theta_yz", geant_tpcFV_theta_yz, "geant_tpcFV_theta_yz[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_mom", geant_tpcFV_mom, "geant_tpcFV_mom[geant_list_size_in_tpcFV]/D");
  fTree->Branch("geant_tpcFV_len", geant_tpcFV_len, "geant_tpcFV_len[geant_list_size_in_tpcFV]/D");

}

void microboone::AnalysisTree::beginSubRun(const art::SubRun& sr)
{

  art::Handle< sumdata::POTSummary > potListHandle;
  //sr.getByLabel(fPOTModuleLabel,potListHandle);
  
  if(sr.getByLabel(fPOTModuleLabel,potListHandle))
    pot=potListHandle->totpot;
  else
    pot=0.;
  
}

void microboone::AnalysisTree::analyze(const art::Event& evt)
{
  
  ResetVars();

  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);

  std::vector< art::Handle< std::vector<recob::Track> > > trackListHandle(kNTracker);
  std::vector< std::vector<art::Ptr<recob::Track> > > tracklist(kNTracker);
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

  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();
  
  art::Timestamp ts = evt.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  evttime = tts.AsDouble();

  //copied from MergeDataPaddles.cxx
  art::Handle< raw::BeamInfo > beam;
  if (evt.getByLabel("beam",beam)){
    beamtime = (double)beam->get_t_ms();
    beamtime/=1000.; //in second
  }

//  std::cout<<detprop->NumberTimeSamples()<<" "<<detprop->ReadOutWindowSize()<<std::endl;
//  std::cout<<geom->DetHalfHeight()*2<<" "<<geom->DetHalfWidth()*2<<" "<<geom->DetLength()<<std::endl;
//  std::cout<<geom->Nwires(0)<<" "<<geom->Nwires(1)<<" "<<geom->Nwires(2)<<std::endl;
  if (evt.isRealData()){
    isdata = 1;
  }
  else isdata = 0;

  //hit information
  no_hits=hitlist.size();
  for (size_t i = 0; i<hitlist.size(); ++i){//loop over hits
    hit_channel[i] = hitlist[i]->Channel();
    hit_plane[i]   = hitlist[i]->WireID().Plane;
    hit_wire[i]    = hitlist[i]->WireID().Wire;
    hit_peakT[i]   = hitlist[i]->PeakTime();
    hit_charge[i]  = hitlist[i]->Charge();
    hit_charge[i]  = hitlist[i]->Charge(true);
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
    ntracks[it1]=tracklist[it1].size(); 
    for(unsigned int i=0; i<tracklist[it1].size();++i){//loop over tracks

      art::Ptr<recob::Track> ptrack(trackListHandle[it1], i);
      const recob::Track& track = *ptrack;
      //we need to use Bezier methods for Bezier tracks
      if(fTrackModuleLabel[it1].compare("beziertracker")==0){
  	trkf::BezierTrack btrack(*ptrack);
  	int ntraj = btrack.NSegments();
  	if(ntraj > 0) {
  	  double xyz[3];
  	  btrack.GetTrackPoint(0,xyz); 
  	  TVector3 pos(xyz[0],xyz[1],xyz[2]);
  	  btrack.GetTrackDirection(0,xyz);  
  	  TVector3 dir_start(xyz[0],xyz[1],xyz[2]);
  	  btrack.GetTrackDirection(1,xyz);  
  	  TVector3 dir_end(xyz[0],xyz[1],xyz[2]);
  	  btrack.GetTrackPoint(1,xyz); 
  	  TVector3 end(xyz[0],xyz[1],xyz[2]);
	 	  
  	  double dpos	  = bdist(pos);
  	  double dend	  = bdist(end);
  	  double tlen	  = btrack.GetLength();
  	  double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
  	  double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());
  	  double mom	  = 0.;
  	  if (btrack.NumberFitMomentum() > 0)	
  	    mom = btrack.VertexMomentum();
  	  // fill bezier track reco branches
  	  trkId[it1][i] 	= i;  //bezier has some screwed up track IDs
  	  trkstartx[it1][i]	= pos.X();
  	  trkstarty[it1][i]	= pos.Y();
  	  trkstartz[it1][i]	= pos.Z();
  	  trkstartd[it1][i]	= dpos;
  	  trkendx[it1][i]	= end.X();
  	  trkendy[it1][i]	= end.Y();
  	  trkendz[it1][i]	= end.Z();
  	  trkendd[it1][i]	= dend;
  	  trktheta[it1][i]	= dir_start.Theta();
  	  trkphi[it1][i]	= dir_start.Phi();
  	  trkstartdcosx[it1][i] = dir_start.X();
  	  trkstartdcosy[it1][i] = dir_start.Y();
  	  trkstartdcosz[it1][i] = dir_start.Z();
  	  trkenddcosx[it1][i]	= dir_end.X();
  	  trkenddcosy[it1][i]	= dir_end.Y();
  	  trkenddcosz[it1][i]	= dir_end.Z();    
  	  trkthetaxz[it1][i]	= theta_xz;
  	  trkthetayz[it1][i]	= theta_yz;
  	  trkmom[it1][i]	= mom;
  	  trklen[it1][i]	= tlen; 
	}
      } 
      else {   //use the normal methods for other kinds of tracks       
        int ntraj = track.NumberTrajectoryPoints();
        if (ntraj > 0) {
  	  const TVector3& pos	   = track.Vertex();
  	  const TVector3& dir_start = track.VertexDirection();
  	  const TVector3& dir_end   = track.EndDirection();
  	  const TVector3& end	   = track.End();
 
   	  double dpos	 = bdist(pos);
  	  double dend	 = bdist(end);
  	  double tlen	 = length(track);
  	  double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
  	  double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());
  	  double mom	 = 0.;
  	  if(track.NumberFitMomentum() > 0)
  	    mom = track.VertexMomentum();
  	  // fill non-bezier-track reco branches
  	  trkId[it1][i]       	     = track.ID();
  	  trkstartx[it1][i]	     = pos.X();
  	  trkstarty[it1][i]	     = pos.Y();
  	  trkstartz[it1][i]	     = pos.Z();
  	  trkstartd[it1][i]	     = dpos;
  	  trkendx[it1][i]     	     = end.X();
  	  trkendy[it1][i]            = end.Y();
  	  trkendz[it1][i]      	     = end.Z();
  	  trkendd[it1][i]            = dend;
  	  trktheta[it1][i]           = dir_start.Theta();
  	  trkphi[it1][i]             = dir_start.Phi();
  	  trkstartdcosx[it1][i]      = dir_start.X();
  	  trkstartdcosy[it1][i]      = dir_start.Y();
  	  trkstartdcosz[it1][i]      = dir_start.Z();	  
  	  trkenddcosx[it1][i]	     = dir_end.X();
  	  trkenddcosy[it1][i]	     = dir_end.Y();
  	  trkenddcosz[it1][i]	     = dir_end.Z();	 
  	  trkthetaxz[it1][i]	     = theta_xz;
  	  trkthetayz[it1][i]	     = theta_yz;
  	  trkmom[it1][i]             = mom;
  	  trklen[it1][i]             = tlen;
        } 
      }
      
      art::FindMany<anab::Calorimetry> fmcal(trackListHandle[it1], evt, fCalorimetryModuleLabel[it1]);
      if (fmcal.isValid()){
        std::vector<const anab::Calorimetry*> calos = fmcal.at(i);
        //std::cout<<"calo size "<<calos.size()<<std::endl;
        for (size_t j = 0; j<calos.size(); ++j){
	  trkke[it1][i][j]    = calos[j]->KineticEnergy();
  	  trkrange[it1][i][j] = calos[j]->Range();
	  trkpitchc[it1][i][j]= calos[j] -> TrkPitchC();
	  ntrkhits[it1][i][j] = calos[j] -> dEdx().size();
	  for(int k = 0; k < ntrkhits[it1][i][j]; ++k) {
	    trkdedx[it1][i][j][k]  = (calos[j] -> dEdx())[k];
	    trkdqdx[it1][i][j][k]  = (calos[j] -> dQdx())[k];
	    trkresrg[it1][i][j][k] = (calos[j] -> ResidualRange())[k];
	    trkxyz[it1][i][j][k][0] = (calos[j] -> XYZ())[k].X();
	    trkxyz[it1][i][j][k][1] = (calos[j] -> XYZ())[k].Y();
	    trkxyz[it1][i][j][k][2] = (calos[j] -> XYZ())[k].Z();
	  }
        }
      }
    
      //track truth information
      if (!isdata){
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
  	  HitsPurity(hits[ipl],trkidtruth[it1][i][ipl],trkpurtruth[it1][i][ipl],maxe);
	  //std::cout<<"\n"<<it1<<"\t"<<i<<"\t"<<ipl<<"\t"<<trkidtruth[it1][i][ipl]<<"\t"<<trkpurtruth[it1][i][ipl]<<"\t"<<maxe;
  	  if (trkidtruth[it1][i][ipl]>0){
  	    const simb::MCParticle *particle = bt->TrackIDToParticle(trkidtruth[it1][i][ipl]);
  	    const std::vector<sim::IDE> vide = bt->TrackIDToSimIDE(trkidtruth[it1][i][ipl]);
  	    double tote = 0;
  	    for (size_t iide = 0; iide<vide.size(); ++iide){
  	      tote += vide[iide].energy;
  	    }
  	    trkpdgtruth[it1][i][ipl] = particle->PdgCode();
  	    trkefftruth[it1][i][ipl] = maxe/(tote/kNplanes); //tote include both induction and collection energies
	    //std::cout<<"\n"<<trkpdgtruth[it1][i][ipl]<<"\t"<<trkefftruth[it1][i][ipl];
  	  }	    
  	}
      }//end if (!isdata)
    }//end loop over track
  }//end loop over track module labels  
  
  //mc truth information
  if (!isdata){//is MC
    const sim::ParticleList& plist = bt->ParticleList();    
    //save neutrino interaction information
    mcevts_truth = mclist.size();
    if (mcevts_truth){//at least one mc record
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
	nuPDG_truth = mctruth->GetNeutrino().Nu().PdgCode();
	ccnc_truth = mctruth->GetNeutrino().CCNC();
	mode_truth = mctruth->GetNeutrino().Mode();
	Q2_truth = mctruth->GetNeutrino().QSqr();
	W_truth = mctruth->GetNeutrino().W();
	hitnuc_truth = mctruth->GetNeutrino().HitNuc();
	enu_truth = mctruth->GetNeutrino().Nu().E();
	nuvtxx_truth = mctruth->GetNeutrino().Nu().Vx();
	nuvtxy_truth = mctruth->GetNeutrino().Nu().Vy();
	nuvtxz_truth = mctruth->GetNeutrino().Nu().Vz();
	if (mctruth->GetNeutrino().Nu().P()){
	  nu_dcosx_truth = mctruth->GetNeutrino().Nu().Px()/mctruth->GetNeutrino().Nu().P();
	  nu_dcosy_truth = mctruth->GetNeutrino().Nu().Py()/mctruth->GetNeutrino().Nu().P();
	  nu_dcosz_truth = mctruth->GetNeutrino().Nu().Pz()/mctruth->GetNeutrino().Nu().P();
	}
	lep_mom_truth = mctruth->GetNeutrino().Lepton().P();
	if (mctruth->GetNeutrino().Lepton().P()){
	  lep_dcosx_truth = mctruth->GetNeutrino().Lepton().Px()/mctruth->GetNeutrino().Lepton().P();
	  lep_dcosy_truth = mctruth->GetNeutrino().Lepton().Py()/mctruth->GetNeutrino().Lepton().P();
	  lep_dcosz_truth = mctruth->GetNeutrino().Lepton().Pz()/mctruth->GetNeutrino().Lepton().P();
	}
	//flux information
	art::Ptr<simb::MCFlux>  mcflux = fluxlist[imc];
	tpx_flux = mcflux->ftpx;
	tpy_flux = mcflux->ftpy;
	tpz_flux = mcflux->ftpz;
	tptype_flux = mcflux->ftptype;

	//genie particles information
	genie_no_primaries=mctruth->NParticles();
      
	for(int j = 0; j < mctruth->NParticles(); ++j){
	  simb::MCParticle part(mctruth->GetParticle(j));
	  
	  genie_primaries_pdg[j]=part.PdgCode();
	  genie_Eng[j]=part.E();
	  genie_Px[j]=part.Px();
	  genie_Py[j]=part.Py();
	  genie_Pz[j]=part.Pz();
	  genie_P[j]=part.Px();
	  genie_status_code[j]=part.StatusCode();
	  genie_mass[j]=part.Mass();
	  genie_trackID[j]=part.TrackId();
	  genie_ND[j]=part.NumberDaughters();
	  genie_mother[j]=part.Mother();
	}
      }
      
      //GEANT particles information
      std::vector<const simb::MCParticle* > geant_part;
      int i = 0;
      geant_list_size_in_tpcFV = 0;
      for(size_t p = 0; p < plist.size(); ++p) 
      	{
        geant_part.push_back(plist.Particle(p)); 
	assert(plist.Particle(p) != 0);
        TVector3 mcstart, mcend;
        double plen = length(*plist.Particle(p), mcstart, mcend);
        if ( (plen==0) || plist.Particle(p)->PdgCode() > 10000) continue;
	else{
	  geant_tpcFV_pdg[i] = plist.Particle(p)->PdgCode();
	  double mctheta_xz = std::atan2(plist.Particle(p)->Px(), plist.Particle(p)->Pz());
	  double mctheta_yz = std::atan2(plist.Particle(p)->Py(), plist.Particle(p)->Pz());

	  geant_tpcFV_trackId[i] = plist.Particle(p)->TrackId();
	  geant_tpcFV_status[i]  = plist.Particle(p)->StatusCode();
	  //
	  geant_tpcFV_orig_E[i]	     = plist.Particle(p)->E();
	  geant_tpcFV_orig_px[i]     = plist.Particle(p)->Px();
	  geant_tpcFV_orig_py[i]     = plist.Particle(p)->Py();
	  geant_tpcFV_orig_pz[i]     = plist.Particle(p)->Pz();
	  geant_tpcFV_orig_startx[i] = plist.Particle(p)->Vx();
	  geant_tpcFV_orig_starty[i] = plist.Particle(p)->Vy();
	  geant_tpcFV_orig_startz[i] = plist.Particle(p)->Vz();
	  geant_tpcFV_orig_startt[i] = plist.Particle(p)->T();
	  geant_tpcFV_orig_endx[i]   = plist.Particle(p)->EndX();
	  geant_tpcFV_orig_endy[i]   = plist.Particle(p)->EndY();
	  geant_tpcFV_orig_endz[i]   = plist.Particle(p)->EndZ();
	  geant_tpcFV_orig_endt[i]   = plist.Particle(p)->EndT();
	  //
	  geant_tpcFV_startx[i]   = mcstart.X();
	  geant_tpcFV_starty[i]   = mcstart.Y();
	  geant_tpcFV_startz[i]   = mcstart.Z();
	  geant_tpcFV_endx[i]	  = mcend.X();
	  geant_tpcFV_endy[i]	  = mcend.Y();
	  geant_tpcFV_endz[i]	  = mcend.Z();
	  geant_tpcFV_theta[i]	  = plist.Particle(p)->Momentum().Theta();
	  geant_tpcFV_phi[i]	  = plist.Particle(p)->Momentum().Phi();
	  geant_tpcFV_theta_xz[i] = mctheta_xz;
	  geant_tpcFV_theta_yz[i] = mctheta_yz;
	  geant_tpcFV_mom[i]	  = plist.Particle(p)->Momentum().Vect().Mag();
	  geant_tpcFV_len[i]	  = plen;    
	  i++;  
        }
      }
      geant_list_size_in_tpcFV = i;				
	
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
      
      no_primaries=primary;
      geant_list_size=geant_particle;
      //std::cout<<"\n"<<geant_list_size<<"\n";
      
      for( unsigned int i = 0; i < geant_part.size(); ++i ){
	
	if(geant_part[i]->Process()==pri){
	  process_primary[i]=1;
	}
	else{
	  process_primary[i]=0;
	}
	
	Mother[i]=geant_part[i]->Mother();
	TrackId[i]=geant_part[i]->TrackId();
	pdg[i]=geant_part[i]->PdgCode();
	Eng[i]=geant_part[i]->E();
	Px[i]=geant_part[i]->Px();
	Py[i]=geant_part[i]->Py();
	Pz[i]=geant_part[i]->Pz();
	StartPointx[i]=geant_part[i]->Vx();
	StartPointy[i]=geant_part[i]->Vy();
	StartPointz[i]=geant_part[i]->Vz();
	EndPointx[i]=geant_part[i]->EndPosition()[0];
	EndPointy[i]=geant_part[i]->EndPosition()[1];
	EndPointz[i]=geant_part[i]->EndPosition()[2];
	NumberDaughters[i]=geant_part[i]->NumberDaughters();
      }
    
    int currentMergedId = 1;
    int currentMotherId = 0;
    int currentMotherTrackId = 0;
    for(int i = 0; i < geant_list_size; i++) {
       MergedId[i] = 0;
    }

    for(int i = geant_list_size - 1; i >= 0; i--) {
       if(MergedId[i] == 0) {
          MergedId[i] = currentMergedId;
          currentMotherId = Mother[i];
          currentMotherTrackId = -1;
          while(currentMotherId > 0) {
             for(int j = 0; j < geant_list_size; j++) {
                if(TrackId[j] == currentMotherId) currentMotherTrackId = j;
             }
             if(pdg[i] == pdg[currentMotherTrackId]) {
                MergedId[currentMotherTrackId] = currentMergedId;
                currentMotherId = Mother[currentMotherTrackId];
             }
             else currentMotherId = 0;
          }
          currentMergedId++;
       }
    }    
   }//if (mcevts_truth){//at least one mc record
  }//if (!isdata){  
  taulife = LArProp->ElectronLifetime();
  fTree->Fill();
}

void microboone::AnalysisTree::HitsPurity(std::vector< art::Ptr<recob::Hit> > const& hits, int& trackid, double& purity, double& maxe){

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
  
  double d1 = pos.X();  			   // Distance to right side (wires).
  double d2 = 2.*geom->DetHalfWidth() - pos.X();   // Distance to left side (cathode).
  double d3 = pos.Y() + geom->DetHalfHeight();     // Distance to bottom.
  double d4 = geom->DetHalfHeight() - pos.Y();     // Distance to top.
  double d5 = pos.Z();  			   // Distance to front.
  double d6 = geom->DetLength() - pos.Z();	   // Distance to back.

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

void microboone::AnalysisTree::ResetVars(){

  run = -99999;
  subrun = -99999;
  event = -99999;
  evttime = -99999;
  beamtime = -99999;
  isdata = -99999;
  taulife = -99999;

  no_hits = 0;
  for (int i = 0; i<kMaxHits; ++i){
    hit_plane[i] = -99999;
    hit_wire[i] = -99999;
    hit_channel[i] = -99999;
    hit_peakT[i] = -99999;
    hit_charge[i] = -99999;
    hit_ph[i] = -99999;
  }

  for (int i = 0; i < kMaxTrackers; ++i){
    ntracks[i] = 0; 
    for (int k=0; k< kMaxHits; k++){
      hit_trkid[i][k] = -99999;
    }    
    for (int j = 0; j < kMaxTrack; ++j){
      trkId[i][j]     = -99999;
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
        ntrkhits[i][j][k]    = -99999;
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
  
}

namespace microboone{

  DEFINE_ART_MODULE(AnalysisTree)
  
} 

#endif

