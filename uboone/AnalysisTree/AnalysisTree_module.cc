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
    int    hit_trkid[kMaxHits];      //is this hit associated with a reco track?

    //track information
    int    ntracks_reco;             //number of reconstructed tracks
    double trkvtxx[kMaxTrack];       //track vertex x
    double trkvtxy[kMaxTrack];       //track vertex y
    double trkvtxz[kMaxTrack];       //track vertex z
    double trkendx[kMaxTrack];       //track end x
    double trkendy[kMaxTrack];       //track end y
    double trkendz[kMaxTrack];       //track end z
    double trkstartdcosx[kMaxTrack]; //track start dcosx
    double trkstartdcosy[kMaxTrack]; //track start dcosy
    double trkstartdcosz[kMaxTrack]; //track start dcosz
    double trkenddcosx[kMaxTrack];   //track end dcosx
    double trkenddcosy[kMaxTrack];   //track end dcosy
    double trkenddcosz[kMaxTrack];   //track end dcosz
    double trkzenith[kMaxTrack];     //track zenith angle
    double trkazimuth[kMaxTrack];    //track azimuth angle
    double trkke[kMaxTrack][kNplanes];
    double trkrange[kMaxTrack][kNplanes];
    int    trkidtruth[kMaxTrack][kNplanes]; //true geant trackid
    double trkpdgtruth[kMaxTrack][kNplanes]; //true pdg code
    double trkefftruth[kMaxTrack][kNplanes]; //completeness
    double trkpurtruth[kMaxTrack][kNplanes]; //purity of track
    
    //*****
    // reconstructed info
    int    reco_trackId[kMaxTrack];
    double reco_startx[kMaxTrack];      // starting x position.
    double reco_starty[kMaxTrack];      // starting y position.
    double reco_startz[kMaxTrack];      // starting z position.
    double reco_startd[kMaxTrack];      // starting distance to boundary.
    double reco_endx[kMaxTrack];        // ending x position.
    double reco_endy[kMaxTrack];        // ending y position.
    double reco_endz[kMaxTrack];        // ending z position.
    double reco_endd[kMaxTrack];        // ending distance to boundary.
    double reco_theta[kMaxTrack];       // theta.
    double reco_phi[kMaxTrack];	      // phi.
    double reco_theta_xz[kMaxTrack];    // theta_xz.
    double reco_theta_yz[kMaxTrack];    // theta_yz.
    double reco_mom[kMaxTrack];	      // momentum.
    double reco_len[kMaxTrack];	      // length.
    //*****
       
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
    
    // **** more geant info
    int   mc_status[kMaxPrimaries];
    int   mc_trackId[kMaxPrimaries];
    int   mc_pdg[kMaxPrimaries];
    
    double mc_p_E[kMaxPrimaries];
    double mc_p_px[kMaxPrimaries];
    double mc_p_py[kMaxPrimaries];
    double mc_p_pz[kMaxPrimaries];
    double mc_p_startx[kMaxPrimaries];
    double mc_p_starty[kMaxPrimaries];
    double mc_p_startz[kMaxPrimaries];
    double mc_p_startt[kMaxPrimaries];
    double mc_p_endx[kMaxPrimaries];
    double mc_p_endy[kMaxPrimaries];
    double mc_p_endz[kMaxPrimaries];
    double mc_p_endt[kMaxPrimaries];
    
    double mc_startx[kMaxPrimaries];	  // starting x position.
    double mc_starty[kMaxPrimaries];	  // starting y position.
    double mc_startz[kMaxPrimaries];	  // starting z position.
    double mc_startd[kMaxPrimaries];	  // starting distance to boundary.
    double mc_endx[kMaxPrimaries];	  // ending x position.
    double mc_endy[kMaxPrimaries];	  // ending y position.
    double mc_endz[kMaxPrimaries];	  // ending z position.
    double mc_endd[kMaxPrimaries];	  // ending distance to boundary.
    double mc_theta[kMaxPrimaries];	  // theta.
    double mc_phi[kMaxPrimaries];	  // phi.
    double mc_theta_xz[kMaxPrimaries];    // theta_xz.
    double mc_theta_yz[kMaxPrimaries];    // theta_yz.
    double mc_mom[kMaxPrimaries];         // momentum.
    double mc_len[kMaxPrimaries];         // length.


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
    std::string fTrackModuleLabel;
    std::string fEndPoint2DModuleLabel;
    std::string fVertexModuleLabel;
    std::string fPOTModuleLabel;
    std::string fCalorimetryModuleLabel;
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
  fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel")        ),
  fEndPoint2DModuleLabel    (pset.get< std::string >("EndPoint2DModuleLabel")   ),
  fVertexModuleLabel        (pset.get< std::string >("VertexModuleLabel")       ),
  fPOTModuleLabel           (pset.get< std::string >("POTModuleLabel")          ),
  fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel")  ),
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
  fTree->Branch("hit_trkid",hit_trkid,"hit_trkid[no_hits]/I");

  fTree->Branch("ntracks_reco",&ntracks_reco,"ntracks_reco/I");
  fTree->Branch("trkvtxx",trkvtxx,"trkvtxx[ntracks_reco]/D");
  fTree->Branch("trkvtxy",trkvtxy,"trkvtxy[ntracks_reco]/D");
  fTree->Branch("trkvtxz",trkvtxz,"trkvtxz[ntracks_reco]/D");
  fTree->Branch("trkendx",trkendx,"trkendx[ntracks_reco]/D");
  fTree->Branch("trkendy",trkendy,"trkendy[ntracks_reco]/D");
  fTree->Branch("trkendz",trkendz,"trkendz[ntracks_reco]/D");
  fTree->Branch("trkstartdcosx",trkstartdcosx,"trkstartdcosx[ntracks_reco]/D");
  fTree->Branch("trkstartdcosy",trkstartdcosy,"trkstartdcosy[ntracks_reco]/D");
  fTree->Branch("trkstartdcosz",trkstartdcosz,"trkstartdcosz[ntracks_reco]/D");
  fTree->Branch("trkenddcosx",trkenddcosx,"trkenddcosx[ntracks_reco]/D");
  fTree->Branch("trkenddcosy",trkenddcosy,"trkenddcosy[ntracks_reco]/D");
  fTree->Branch("trkenddcosz",trkenddcosz,"trkenddcosz[ntracks_reco]/D");
  fTree->Branch("trkzenith",trkzenith,"trkzenith[ntracks_reco]/D");
  fTree->Branch("trkazimuth",trkazimuth,"trkazimuth[ntracks_reco]/D");
  fTree->Branch("trkke",trkke,"trkke[ntracks_reco][3]/D");
  fTree->Branch("trkrange",trkrange,"trkrange[ntracks_reco][3]/D");
  fTree->Branch("trkidtruth",trkidtruth,"trkidtruth[ntracks_reco][3]/I");
  fTree->Branch("trkpdgtruth",trkpdgtruth,"trkpdgtruth[ntracks_reco][3]/D");
  fTree->Branch("trkefftruth",trkefftruth,"trkefftruth[ntracks_reco][3]/D");
  fTree->Branch("trkpurtruth",trkpurtruth,"trkpurtruth[ntracks_reco][3]/D");

  //*****more track information  
  fTree->Branch("reco_trackId", reco_trackId, "reco_trackId[ntracks_reco]/I");
  fTree->Branch("reco_startx", reco_startx, "reco_startx[ntracks_reco]/D");
  fTree->Branch("reco_starty", reco_starty, "reco_starty[ntracks_reco]/D");
  fTree->Branch("reco_startz", reco_startz, "reco_startz[ntracks_reco]/D");
  fTree->Branch("reco_startd", reco_startd, "reco_startd[ntracks_reco]/D");
  fTree->Branch("reco_endx", reco_endx, "reco_endx[ntracks_reco]/D");
  fTree->Branch("reco_endy", reco_endy, "reco_endy[ntracks_reco]/D");
  fTree->Branch("reco_endz", reco_endz, "reco_endz[ntracks_reco]/D");
  fTree->Branch("reco_endd", reco_endd, "reco_endd[ntracks_reco]/D");
  fTree->Branch("reco_theta", reco_theta, "reco_theta[ntracks_reco]/D");
  fTree->Branch("reco_phi", reco_phi, "reco_phi[ntracks_reco]/D");
  fTree->Branch("reco_theta_xz", reco_theta_xz, "reco_theta_xz[ntracks_reco]/D");
  fTree->Branch("treco_theta_yz", reco_theta_yz, "reco_theta_yz[ntracks_reco]/D");
  fTree->Branch("reco_mom", reco_mom, "reco_mom[ntracks_reco]/D");
  fTree->Branch("reco_len", reco_len, "reco_len[ntracks_reco]/D");  
  //*****  

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

  fTree->Branch("mc_pdg", mc_pdg, "mc_pdg[geant_list_size_in_tpcFV]/I");
  fTree->Branch("mc_status", mc_status, "mc_status[geant_list_size_in_tpcFV]/I");
  fTree->Branch("mc_trackId", mc_trackId, "mc_trackId[geant_list_size_in_tpcFV]/I");
  fTree->Branch("mc_p_E", mc_p_E, "mc_p_E[geant_list_size_in_tpcFV]/F");
  fTree->Branch("mc_p_px", mc_p_px, "mc_p_px[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_p_py", mc_p_py, "mc_p_py[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_p_pz", mc_p_pz, "mc_p_pz[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_p_startx", mc_p_startx, "mc_p_startx[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_p_starty", mc_p_starty, "mc_p_starty[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_p_startz", mc_p_startz, "mc_p_startz[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_p_startt", mc_p_startt, "mc_p_startt[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_p_endx", mc_p_endx, "mc_p_endx[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_p_endy", mc_p_endy, "mc_p_endy[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_p_endz", mc_p_endz, "mc_p_endz[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_p_endt", mc_p_endt, "mc_p_endt[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_startx", mc_startx, "mc_startx[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_starty", mc_starty, "mc_starty[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_startz", mc_startz, "mc_startz[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_startd", mc_startd, "mc_startd[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_endx", mc_endx, "mc_endx[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_endy", mc_endy, "mc_endy[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_endz", mc_endz, "mc_endz[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_endd", mc_endd, "mc_endd[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_theta", mc_theta, "mc_theta[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_phi", mc_phi, "mc_phi[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_theta_xz", mc_theta_xz, "mc_theta_xz[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_theta_yz", mc_theta_yz, "mc_theta_yz[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_mom", mc_mom, "mc_mom[geant_list_size_in_tpcFV]/D");
  fTree->Branch("mc_len", mc_len, "mc_len[geant_list_size_in_tpcFV]/D");

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

  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);

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

  //associations
  art::FindMany<recob::Track>     fmtk(hitListHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Hit>      fmht(trackListHandle, evt, fTrackModuleLabel);
  art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);

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
    if (fmtk.at(i).size()!=0){
      hit_trkid[i] = fmtk.at(i)[0]->ID();
    }
  }

  //track information

  ntracks_reco=tracklist.size();
  double larStart[3];
  double larEnd[3];
  std::vector<double> trackStart;
  std::vector<double> trackEnd;

  for(unsigned int i=0; i<tracklist.size();++i){//loop over tracks
    trackStart.clear();
    trackEnd.clear();
    memset(larStart, 0, 3);
    memset(larEnd, 0, 3);
    tracklist[i]->Extent(trackStart,trackEnd); 
    tracklist[i]->Direction(larStart,larEnd);
    trkvtxx[i]        = trackStart[0];
    trkvtxy[i]        = trackStart[1];
    trkvtxz[i]        = trackStart[2];
    trkendx[i]        = trackEnd[0];
    trkendy[i]        = trackEnd[1];
    trkendz[i]        = trackEnd[2];
    trkstartdcosx[i]  = larStart[0];
    trkstartdcosy[i]  = larStart[1];
    trkstartdcosz[i]  = larStart[2];
    trkenddcosx[i]    = larEnd[0];
    trkenddcosy[i]    = larEnd[1];
    trkenddcosz[i]    = larEnd[2];
    
    std::vector< art::Ptr<anab::Calorimetry> > calos = fmcal.at(i);
    //std::cout<<"calo size "<<calos.size()<<std::endl;
    for (size_t j = 0; j<calos.size(); ++j){
      trkke[i][j] = calos[j]->KineticEnergy();
      trkrange[i][j] = calos[j]->Range();
    }
    //track truth information
    if (!isdata){
      // get the hits on each plane
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
	HitsPurity(hits[ipl],trkidtruth[i][ipl],trkpurtruth[i][ipl],maxe);
	if (trkidtruth[i][ipl]>0){
	  const simb::MCParticle *particle = bt->TrackIDToParticle(trkidtruth[i][ipl]);
	  const std::vector<sim::IDE> vide = bt->TrackIDToSimIDE(trkidtruth[i][ipl]);
	  double tote = 0;
	  for (size_t iide = 0; iide<vide.size(); ++iide){
	    tote += vide[iide].energy;
	  }
	  trkpdgtruth[i][ipl] = particle->PdgCode();
	  trkefftruth[i][ipl] = maxe/(tote/kNplanes); //I believe tote include both induction and collection energies
	}	  
      }
    }
  //*****sowjanya*****track information//


  art::Ptr<recob::Track> ptrack(trackListHandle, i);
    const recob::Track& track = *ptrack;    
    //we need to use Bezier methods for Bezier tracks

    /**/
    if(fTrackModuleLabel.compare("beziertracker")==0){
      trkf::BezierTrack btrack(*ptrack);
      int ntraj = btrack.NSegments();
      if(ntraj > 0) {
    	double xyz[3];
    	btrack.GetTrackPoint(0,xyz); 
	TVector3 pos(xyz[0],xyz[1],xyz[2]);
    	btrack.GetTrackDirection(0,xyz);  
	TVector3 dir(xyz[0],xyz[1],xyz[2]);
    	btrack.GetTrackPoint(1,xyz); 
	TVector3 end(xyz[0],xyz[1],xyz[2]);
    	
    	double dpos	= bdist(pos);
    	double dend	= bdist(end);
    	double tlen	= btrack.GetLength();
    	double theta_xz = std::atan2(dir.X(), dir.Z());
    	double theta_yz = std::atan2(dir.Y(), dir.Z());
    	double mom	= 0.;
    	if (btrack.NumberFitMomentum() > 0)   
    	  mom = btrack.VertexMomentum();
    	
    	// fill bezier track reco branches
    	reco_trackId[i]  = i;  //bezier has some screwed up track IDs
    	reco_startx[i]   = pos.X();
    	reco_starty[i]   = pos.Y();
    	reco_startz[i]   = pos.Z();
    	reco_startd[i]   = dpos;
    	reco_endx[i]     = end.X();
    	reco_endy[i]     = end.Y();
    	reco_endz[i]     = end.Z();
    	reco_startd[i]   = dend;
    	reco_theta[i]    = dir.Theta();
    	reco_phi[i]      = dir.Phi();
    	reco_theta_xz[i] = theta_xz;
    	reco_theta_yz[i] = theta_yz;
    	reco_mom[i]      = mom;
    	reco_len[i]      = tlen;  	
      }
   } 
   else {   //use the normal methods for other kinds of tracks       
     /**/
     //{
     int ntraj = track.NumberTrajectoryPoints();
     if (ntraj > 0) {
       const TVector3& pos = track.Vertex();
       const TVector3& dir = track.VertexDirection();
       const TVector3& end = track.End();

       double dpos     = bdist(pos);
       double dend     = bdist(end);
       double tlen     = length(track);
       double theta_xz = std::atan2(dir.X(), dir.Z());
       double theta_yz = std::atan2(dir.Y(), dir.Z());
       double mom      = 0.;
       if(track.NumberFitMomentum() > 0)
         mom = track.VertexMomentum();

       // fill non-bezier-track reco branches
       reco_trackId[i]  = track.ID();
       reco_startx[i]	= pos.X();
       reco_starty[i]	= pos.Y();
       reco_startz[i]	= pos.Z();
       reco_startd[i]	= dpos;
       reco_endx[i]	= end.X();
       reco_endy[i]	= end.Y();
       reco_endz[i]	= end.Z();
       reco_startd[i]	= dend;
       reco_theta[i]	= dir.Theta();
       reco_phi[i]	= dir.Phi();
       reco_theta_xz[i] = theta_xz;
       reco_theta_yz[i] = theta_yz;
       reco_mom[i]	= mom;
       reco_len[i]	= tlen;
     } 
   }    
   //*****sowjanya*****track information//
 }//end loop over track

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
      for(size_t p = 0; p < plist.size(); ++p) 
        geant_part.push_back(plist.Particle(p)); 		
	
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
      //*******
    
      std::vector<const simb::MCParticle* > geant_part_intpcFV;      
       for(size_t p = 0; p < plist.size(); ++p) { 
        assert(plist.Particle(p) != 0);
        TVector3 mcstart, mcend;
        double plen = length(*plist.Particle(p), mcstart, mcend);
        if ( (plen==0) || plist.Particle(p)->PdgCode() > 10000) continue;
		//std::cout<<plist.Particle(p)->PdgCode()<<"\n";
	geant_part_intpcFV.push_back(plist.Particle(p));
    }  
    
    geant_list_size_in_tpcFV = geant_part_intpcFV.size();      
     

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


      for( unsigned int i = 0; i < geant_part_intpcFV.size(); ++i ){ 
      TVector3 mcstart;
      TVector3 mcend;
      double plen = length(*geant_part_intpcFV[i], mcstart, mcend);

      mc_pdg[i] = geant_part_intpcFV[i]->PdgCode();

      double mctheta_xz = std::atan2(geant_part_intpcFV[i]->Px(), geant_part_intpcFV[i]->Pz());
      double mctheta_yz = std::atan2(geant_part_intpcFV[i]->Py(), geant_part_intpcFV[i]->Pz());

      mc_trackId[i] = geant_part_intpcFV[i]->TrackId();
      mc_status[i]  = geant_part_intpcFV[i]->StatusCode();
      //
      mc_p_E[i]      = geant_part_intpcFV[i]->E();
      mc_p_px[i]     = geant_part_intpcFV[i]->Px();
      mc_p_py[i]     = geant_part_intpcFV[i]->Py();
      mc_p_pz[i]     = geant_part_intpcFV[i]->Pz();
      mc_p_startx[i] = geant_part_intpcFV[i]->Vx();
      mc_p_starty[i] = geant_part_intpcFV[i]->Vy();
      mc_p_startz[i] = geant_part_intpcFV[i]->Vz();
      mc_p_startt[i] = geant_part_intpcFV[i]->T();
      mc_p_endx[i]   = geant_part_intpcFV[i]->EndX();
      mc_p_endy[i]   = geant_part_intpcFV[i]->EndY();
      mc_p_endz[i]   = geant_part_intpcFV[i]->EndZ();
      mc_p_endt[i]   = geant_part_intpcFV[i]->EndT();
      //
      mc_startx[i]   = mcstart.X();
      mc_starty[i]   = mcstart.Y();
      mc_startz[i]   = mcstart.Z();
      mc_endx[i]     = mcend.X();
      mc_endy[i]     = mcend.Y();
      mc_endz[i]     = mcend.Z();
      mc_theta[i]    = geant_part_intpcFV[i]->Momentum().Theta();
      mc_phi[i]      = geant_part_intpcFV[i]->Momentum().Phi();
      mc_theta_xz[i] = mctheta_xz;
      mc_theta_yz[i] = mctheta_yz;
      mc_mom[i]      = geant_part_intpcFV[i]->Momentum().Vect().Mag();
      mc_len[i]      = plen;      
    }    
    
    //***** 
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
    hit_trkid[i] = -99999;
  }

  ntracks_reco = 0;
  for (int i = 0; i < kMaxTrack; ++i){
    trkvtxx[i] = -99999;
    trkvtxy[i] = -99999;
    trkvtxz[i] = -99999;
    trkendx[i] = -99999;
    trkendy[i] = -99999;
    trkendz[i] = -99999;
    trkstartdcosx[i] = -99999;
    trkstartdcosy[i] = -99999;
    trkstartdcosz[i] = -99999;
    trkenddcosx[i] = -99999;
    trkenddcosy[i] = -99999;
    trkenddcosz[i] = -99999;
    trkzenith[i] = -99999;
    trkazimuth[i] = -99999;
    
    reco_trackId[i] = -99999;
    reco_startx[i] = -99999;   
    reco_starty[i] = -99999;   
    reco_startz[i] = -99999;   
    reco_startd[i] = -99999;   
    reco_endx[i] = -99999;     
    reco_endy[i] = -99999;     
    reco_endz[i] = -99999;     
    reco_endd[i] = -99999;     
    reco_theta[i] = -99999;    
    reco_phi[i] = -99999;      
    reco_theta_xz[i] = -99999; 
    reco_theta_yz[i] = -99999; 
    reco_mom[i] = -99999;      
    reco_len[i] = -99999;      
    
    for (int j = 0; j < kNplanes; ++j){
      trkke[i][j] = -99999;
      trkrange[i][j] = -99999;
      trkidtruth[i][j] = -99999;
      trkpdgtruth[i][j] = -99999;
      trkefftruth[i][j] = -99999;
      trkpurtruth[i][j] = -99999;
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
    
    mc_status[i] = -99999;
    mc_trackId[i] = -99999;
    mc_pdg[i] = -99999;

    mc_p_E[i] = -99999;
    mc_p_px[i] = -99999;
    mc_p_py[i] = -99999;
    mc_p_pz[i] = -99999;
    mc_p_startx[i] = -99999;
    mc_p_starty[i] = -99999;
    mc_p_startz[i] = -99999;
    mc_p_startt[i] = -99999;
    mc_p_endx[i] = -99999;
    mc_p_endy[i] = -99999;
    mc_p_endz[i] = -99999;
    mc_p_endt[i] = -99999;

    mc_startx[i] = -99999;  
    mc_starty[i] = -99999;  
    mc_startz[i] = -99999;  
    mc_startd[i] = -99999;  
    mc_endx[i] = -99999;    
    mc_endy[i] = -99999;    
    mc_endz[i] = -99999;    
    mc_endd[i] = -99999;    
    mc_theta[i] = -99999;   
    mc_phi[i] = -99999;     
    mc_theta_xz[i] = -99999;
    mc_theta_yz[i] = -99999;
    mc_mom[i] = -99999;     
    mc_len[i] = -99999;     
  }    
  
}

namespace microboone{

  DEFINE_ART_MODULE(AnalysisTree)
  
} 

#endif

