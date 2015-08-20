////////////////////////////////////////////////////////////////////////
//
// module to create a TTree for Diffusion analysis
//
// \author sowjanyag@phys.ksu.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef DIFFUSION_H
#define DIFFUSION_H

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
#include "RawData/raw.h"
#include "RawData/BeamInfo.h"
#include "Utilities/LArProperties.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"
#include "SummaryData/POTSummary.h"
#include "MCCheater/BackTracker.h"
#include "RecoBase/Track.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Wire.h"
#include "RecoBase/EndPoint2D.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/OpFlash.h"
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

const int kNplanes         = 3;     //number of wire planes
const int kMaxTrack        = 1000;  //maximum number of tracks
const int kMaxCluster      = 1500;  //maximum number of tracks
const int kMaxHits         = 25000; //maximum number of hits;
const int kMaxPrimaries    = 20000;  //maximum number of primary particles
const int kMaxTrackHits    = 2000;  //maximum number of hits on a track
const int kMaxClusterHits  = 2000;  //maximum number of hits on a track

namespace microboone {
   
  class Diffusion : public art::EDAnalyzer {

  public:
          
    explicit Diffusion(fhicl::ParameterSet const& pset); 
    virtual ~Diffusion();
 
    void analyze(const art::Event& evt);
    void beginJob();
    void beginSubRun(const art::SubRun& sr);

  private:
    
    void   HitsPurity(std::vector< art::Ptr<recob::Hit> > const& hits, Int_t& trackid, Float_t& purity, double& maxe);
    void   ResetVars();
    double length(const recob::Track& track);
    double length(const simb::MCParticle& part, TVector3& start, TVector3& end);

    TTree* fTree;
    //run information
    Int_t    run;                  //run number
    Int_t    subrun;               //subrun number
    Int_t    event;                //event number
    Double_t evttime;              //event time in sec
    Double_t beamtime;             //beam time
    Double_t pot;                  //protons on target
    Double_t taulife;              //electron lifetime
    Char_t   isdata;               //flag, 0=MC 1=data

    // Hit information 
    Int_t    no_hits;                  //number of hits
    Short_t  hit_plane[kMaxHits];      //plane number
    Short_t  hit_wire[kMaxHits];       //wire number
    Short_t  hit_channel[kMaxHits];    //channel ID
    Float_t  hit_peakT[kMaxHits];      //peak time
    Float_t  hit_charge[kMaxHits];     //charge (area)
    Float_t  hit_ph[kMaxHits];         //amplitude
    Float_t  hit_startT[kMaxHits];     //hit start time
    Float_t  hit_endT[kMaxHits];       //hit end time
    Float_t  hit_starttick[kMaxHits];
    Float_t  hit_endtick[kMaxHits];
    Float_t  hit_RMS[kMaxHits];
    Float_t  hit_nelec[kMaxHits];     //hit number of electrons
    Float_t  hit_energy[kMaxHits];       //hit energy
    Short_t  hit_trkid[kMaxHits];      //is this hit associated with a reco track?
    Short_t  hit_clusterid[kMaxHits];  //is this hit associated with a reco cluster?
    Float_t hit_trueX[kMaxHits];
    Int_t hit_ID[kMaxHits]; 

    Int_t rawD_pk[kMaxHits];  
    Int_t rawD_t[kMaxHits];  
    Int_t rawD_ch[kMaxHits];  
    Int_t rawD_fwhh[kMaxHits];  
    Double_t rawD_rms[kMaxHits]; 
    //wire (CalWire) information 
    Double_t wire_pk[kMaxHits];  
    Int_t wire_t[kMaxHits];  
    Double_t wire_ch[kMaxHits];  
    Double_t wire_rms[kMaxHits]; 
    //nelectrons (SimChannel) information
    Double_t sim_pk[kMaxHits];  
    Int_t sim_tdc[kMaxHits];  
    Double_t sim_ch[kMaxHits];  
    Double_t sim_rms[kMaxHits]; 

    //Cluster Information
    Short_t nclusters;
    Short_t clusterId[kMaxCluster];
    Short_t clusterView[kMaxCluster];
    Int_t   cluster_isValid[kMaxCluster];
    Float_t cluster_StartCharge[kMaxCluster];
    Float_t cluster_StartAngle[kMaxCluster];
    Float_t cluster_EndCharge[kMaxCluster];
    Float_t cluster_EndAngle[kMaxCluster];
    Float_t cluster_Integral[kMaxCluster];
    Float_t cluster_IntegralAverage[kMaxCluster];
    Float_t cluster_SummedADC[kMaxCluster];
    Float_t cluster_SummedADCaverage[kMaxCluster];
    Float_t cluster_MultipleHitDensity[kMaxCluster];
    Float_t cluster_Width[kMaxCluster];
    Short_t cluster_NHits[kMaxCluster];
    Short_t cluster_StartWire[kMaxCluster];
    Short_t cluster_StartTick[kMaxCluster];
    Short_t cluster_EndWire[kMaxCluster];
    Short_t cluster_EndTick[kMaxCluster];  
    
    Short_t  ncluhits[kMaxCluster][kNplanes];
    
    Short_t  cluhit_plane[kMaxCluster][kNplanes][kMaxClusterHits];     
    Short_t  cluhit_wire[kMaxCluster][kNplanes][kMaxClusterHits];      
    Short_t  cluhit_channel[kMaxCluster][kNplanes][kMaxClusterHits];   
    Float_t  cluhit_peakT[kMaxCluster][kNplanes][kMaxClusterHits];     
    Float_t  cluhit_charge[kMaxCluster][kNplanes][kMaxClusterHits];    
    Float_t  cluhit_ph[kMaxCluster][kNplanes][kMaxClusterHits];        
    Float_t  cluhit_startT[kMaxCluster][kNplanes][kMaxClusterHits];    
    Float_t  cluhit_endT[kMaxCluster][kNplanes][kMaxClusterHits];      
    Float_t  cluhit_starttick[kMaxCluster][kNplanes][kMaxClusterHits];
    Float_t  cluhit_endtick[kMaxCluster][kNplanes][kMaxClusterHits];
    Float_t  cluhit_RMS[kMaxCluster][kNplanes][kMaxClusterHits];
    Float_t  cluhit_nelec[kMaxCluster][kNplanes][kMaxClusterHits];     
    Float_t  cluhit_energy[kMaxCluster][kNplanes][kMaxClusterHits];    
    Short_t  cluhit_trkid[kMaxCluster][kNplanes][kMaxClusterHits];     
    Short_t  cluhit_clusterid[kMaxCluster][kNplanes][kMaxClusterHits]; 
    Int_t   cluhit_ID[kMaxCluster][kNplanes][kMaxClusterHits]; 

    //Track information
    //Track plane data
    Short_t    ntracks;               
    Float_t    trkke[kMaxTrack][kNplanes];
    Float_t    trkrange[kMaxTrack][kNplanes];
    Int_t      trkidtruth[kMaxTrack][kNplanes];  //true geant trackid
    Int_t      trkpdgtruth[kMaxTrack][kNplanes]; //true pdg code
    Float_t    trkpitchc[kMaxTrack][kNplanes];
    Short_t    ntrkhits[kMaxTrack][kNplanes];
    //Track Plane Hit data
    Float_t	  trkdedx[kMaxTrack][kNplanes][kMaxTrackHits];
    Float_t	  trkdqdx[kMaxTrack][kNplanes][kMaxTrackHits];
    Float_t	  trkresrg[kMaxTrack][kNplanes][kMaxTrackHits];
    //Track Plane Hit Coordinate data
    Float_t       trkxyz[kMaxTrack][kNplanes][kMaxTrackHits][3];

    // more track info
    Short_t trkId[kMaxTrack];
    Float_t trkstartx[kMaxTrack];     // starting x position.
    Float_t trkstarty[kMaxTrack];     // starting y position.
    Float_t trkstartz[kMaxTrack];     // starting z position.
    Float_t trkendx[kMaxTrack];	   // ending x position.
    Float_t trkendy[kMaxTrack];	   // ending y position.
    Float_t trkendz[kMaxTrack];	   // ending z position.
    Float_t trktheta[kMaxTrack];	   // theta.
    Float_t trkphi[kMaxTrack];	   // phi.
    Float_t trkstartdcosx[kMaxTrack];
    Float_t trkstartdcosy[kMaxTrack];
    Float_t trkstartdcosz[kMaxTrack];
    Float_t trkenddcosx[kMaxTrack];
    Float_t trkenddcosy[kMaxTrack];
    Float_t trkenddcosz[kMaxTrack];
    Float_t trkthetaxz[kMaxTrack];    // theta_xz.
    Float_t trkthetayz[kMaxTrack];    // theta_yz.
    Float_t trkmom[kMaxTrack];	   // momentum.
    Float_t trklen[kMaxTrack];	   // length.
    Float_t trkmomrange[kMaxTrack];    // track momentum from range using CSDA tables
    Float_t trkmommschi2[kMaxTrack];   // track momentum from multiple scattering Chi2 method
    Float_t trkmommsllhd[kMaxTrack];   // track momentum from multiple scattering LLHD method
    Int_t   trkpidpdg[kMaxTrack][kNplanes];	    // particle PID pdg code
    Float_t trkpidchi[kMaxTrack][kNplanes];
    Float_t trkpidpida[kMaxTrack][kNplanes];    // particle PIDA
    Short_t trkpidbestplane[kMaxTrack]; // this is defined as the plane with most hits	
    
    //Geant information
    Int_t    no_primaries;      //number of primary geant particles
    Int_t    geant_list_size;  //number of all geant particles
    Int_t    geant_list_size_in_tpcAV;
    Int_t    pdg[kMaxPrimaries];
    Int_t    status[kMaxPrimaries];	 
    Float_t  Eng[kMaxPrimaries];
    Float_t  EndE[kMaxPrimaries];
    Float_t  Mass[kMaxPrimaries];
    Float_t  Px[kMaxPrimaries];
    Float_t  Py[kMaxPrimaries];
    Float_t  Pz[kMaxPrimaries];
    Float_t  P[kMaxPrimaries];
    Float_t  StartPointx[kMaxPrimaries];
    Float_t  StartPointy[kMaxPrimaries];
    Float_t  StartPointz[kMaxPrimaries];
    Float_t  StartT[kMaxPrimaries];  
    Float_t  EndT[kMaxPrimaries];	    
    Float_t  EndPointx[kMaxPrimaries];
    Float_t  EndPointy[kMaxPrimaries];
    Float_t  EndPointz[kMaxPrimaries];
    Float_t  theta[kMaxPrimaries];    
    Float_t  phi[kMaxPrimaries];    
    Float_t  theta_xz[kMaxPrimaries];    
    Float_t  theta_yz[kMaxPrimaries];    
    Float_t  pathlen[kMaxPrimaries];	 
    Int_t    inTPCActive[kMaxPrimaries];    
    Float_t  StartPointx_tpcAV[kMaxPrimaries];
    Float_t  StartPointy_tpcAV[kMaxPrimaries];
    Float_t  StartPointz_tpcAV[kMaxPrimaries];
    Float_t  EndPointx_tpcAV[kMaxPrimaries];
    Float_t  EndPointy_tpcAV[kMaxPrimaries];
    Float_t  EndPointz_tpcAV[kMaxPrimaries];
    Int_t    NumberDaughters[kMaxPrimaries];
    Int_t    TrackId[kMaxPrimaries];
    Int_t    Mother[kMaxPrimaries];
    Int_t    process_primary[kMaxPrimaries];
    std::string processname[kMaxPrimaries];
    
    std::string fDigitModuleLabel;
    std::string fHitsModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fTrackModuleLabel;
    std::string fPOTModuleLabel;
    std::string fCalorimetryModuleLabel;
    std::string fParticleIDModuleLabel;
    std::string fClusterModuleLabel;
    bool  fSaveCaloInfo;
    bool  fSaveTrackInfo;
    bool  fSaveClusterInfo;
    bool  fSaveClusterHitInfo;
    float fG4minE;
  };
}

microboone::Diffusion::Diffusion(fhicl::ParameterSet const& pset) : 
  EDAnalyzer(pset),
  fDigitModuleLabel         (pset.get< std::string >("DigitModuleLabel")        ),
  fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel")         ),
  fLArG4ModuleLabel         (pset.get< std::string >("LArGeantModuleLabel")     ),
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel")     ),  
  fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel")        ),
  fPOTModuleLabel           (pset.get< std::string >("POTModuleLabel")          ),
  fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel")  ),
  fParticleIDModuleLabel    (pset.get< std::string >("ParticleIDModuleLabel")   ),
  fClusterModuleLabel       (pset.get< std::string >("ClusterModuleLabel")      ),  
  fSaveCaloInfo             (pset.get< bool>("SaveCaloInfo",false)),
  fSaveTrackInfo            (pset.get< bool>("SaveTrackInfo",false)),  
  fSaveClusterInfo          (pset.get< bool>("SaveClusterInfo",false)),
  fSaveClusterHitInfo       (pset.get< bool>("SaveClusterHitInfo",false)),
  fG4minE                   (pset.get< float>("G4minE",0.01))  
{
  if (fSaveTrackInfo == false) fSaveCaloInfo = false;
}

//-------------------------------------------------
microboone::Diffusion::~Diffusion()
{
}

void microboone::Diffusion::beginJob(){

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("diffusiontree","diffusion");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("evttime",&evttime,"evttime/D");
  fTree->Branch("beamtime",&beamtime,"beamtime/D");
  fTree->Branch("pot",&pot,"pot/D");
  fTree->Branch("isdata",&isdata,"isdata/B");
  fTree->Branch("taulife",&taulife,"taulife/D");

  fTree->Branch("no_hits",&no_hits,"no_hits/I");
  fTree->Branch("hit_plane",hit_plane,"hit_plane[no_hits]/S");
  fTree->Branch("hit_wire",hit_wire,"hit_wire[no_hits]/S");
  fTree->Branch("hit_channel",hit_channel,"hit_channel[no_hits]/S");
  fTree->Branch("hit_peakT",hit_peakT,"hit_peakT[no_hits]/F");
  fTree->Branch("hit_charge",hit_charge,"hit_charge[no_hits]/F");
  fTree->Branch("hit_ph",hit_ph,"hit_ph[no_hits]/F");
  fTree->Branch("hit_startT",hit_startT,"hit_startT[no_hits]/F");
  fTree->Branch("hit_endT",hit_endT,"hit_endT[no_hits]/F");
  fTree->Branch("hit_starttick",hit_starttick,"hit_starttick[no_hits]/F");
  fTree->Branch("hit_endtick",hit_endtick,"hit_endtick[no_hits]/F");
  fTree->Branch("hit_RMS",hit_RMS,"hit_RMS[no_hits]/F");
  fTree->Branch("hit_ID",hit_ID,"hit_ID[no_hits]/I");  
  fTree->Branch("hit_trueX",hit_trueX,"hit_trueX[no_hits]/F");   
  fTree->Branch("hit_trkid",hit_trkid,"hit_trkid[no_hits]/S");
  fTree->Branch("hit_clusterid",hit_clusterid,"hit_clusterid[no_hits]/S");
  fTree->Branch("hit_nelec",hit_nelec,"hit_nelec[no_hits]/F");
  fTree->Branch("hit_energy",hit_energy,"hit_energy[no_hits]/F");

  fTree->Branch("rawD_pk",rawD_pk,"rawD_pk[no_hits]/I");  
  fTree->Branch("rawD_t",rawD_t,"rawD_t[no_hits]/I");
  fTree->Branch("rawD_ch",rawD_ch,"rawD_ch[no_hits]/I");
  fTree->Branch("rawD_fwhh",rawD_fwhh,"rawD_fwhh[no_hits]/I");
  fTree->Branch("rawD_rms",rawD_rms,"rawD_rms[no_hits]/D");  

  fTree->Branch("wire_pk",wire_pk,"wire_pk[no_hits]/D");
  fTree->Branch("wire_t",wire_t,"wire_t[no_hits]/I");
  fTree->Branch("wire_ch",wire_ch,"wire_ch[no_hits]/D");
  fTree->Branch("wire_rms",wire_rms,"wire_rms[no_hits]/D");

  fTree->Branch("sim_pk",sim_pk,"sim_pk[no_hits]/D");
  fTree->Branch("sim_tdc",sim_tdc,"sim_tdc[no_hits]/I");
  fTree->Branch("sim_ch",sim_ch,"sim_ch[no_hits]/D");
  fTree->Branch("sim_rms",sim_rms,"sim_rms[no_hits]/D");
  
  
  fTree->Branch("nclusters",&nclusters,"nclusters/S");
  fTree->Branch("clusterId", clusterId, "clusterId[nclusters]/S");
  fTree->Branch("clusterView", clusterView, "clusterView[nclusters]/S");
  fTree->Branch("cluster_StartCharge", cluster_StartCharge, "cluster_StartCharge[nclusters]/F");
  fTree->Branch("cluster_StartAngle", cluster_StartAngle, "cluster_StartAngle[nclusters]/F");
  fTree->Branch("cluster_EndCharge", cluster_EndCharge, "cluster_EndCharge[nclusters]/F");
  fTree->Branch("cluster_EndAngle", cluster_EndAngle, "cluster_EndAngle[nclusters]/F");
  fTree->Branch("cluster_Integral", cluster_Integral, "cluster_Integral[nclusters]/F");
  fTree->Branch("cluster_IntegralAverage", cluster_IntegralAverage, "cluster_IntegralAverage[nclusters]/F");
  fTree->Branch("cluster_SummedADC", cluster_SummedADC, "cluster_SummedADC[nclusters]/F");
  fTree->Branch("cluster_SummedADCaverage", cluster_SummedADCaverage, "cluster_SummedADCaverage[nclusters]/F");
  fTree->Branch("cluster_MultipleHitDensity", cluster_MultipleHitDensity, "cluster_MultipleHitDensity[nclusters]/F");
  fTree->Branch("cluster_Width", cluster_Width, "cluster_Width[nclusters]/F");
  fTree->Branch("cluster_NHits", cluster_NHits, "cluster_NHits[nclusters]/S");
  fTree->Branch("cluster_StartWire", cluster_StartWire, "cluster_StartWire[nclusters]/S");
  fTree->Branch("cluster_StartTick", cluster_StartTick, "cluster_StartTick[nclusters]/S");
  fTree->Branch("cluster_EndWire", cluster_EndWire, "cluster_EndWire[nclusters]/S");
  fTree->Branch("cluster_EndTick", cluster_EndTick, "cluster_EndTick[nclusters]/S");
  
  fTree->Branch("ncluhits",&ncluhits,"ncluhits[nclusters][3]/S");
  if (fSaveClusterHitInfo){
    fTree->Branch("cluhit_plane",cluhit_plane,"cluhit_plane[nclusters][3][2000]/S");
    fTree->Branch("cluhit_wire",cluhit_wire,"cluhit_wire[nclusters][3][2000]/S");
    fTree->Branch("cluhit_channel",cluhit_channel,"cluhit_channel[nclusters][3][2000]/S");
    fTree->Branch("cluhit_peakT",cluhit_peakT,"cluhit_peakT[nclusters][3][2000]/F");
    fTree->Branch("cluhit_charge",cluhit_charge,"cluhit_charge[nclusters][3][2000]/F");
    fTree->Branch("cluhit_ph",cluhit_ph,"cluhit_ph[nclusters][3][2000]/F");
    fTree->Branch("cluhit_startT",cluhit_startT,"cluhit_startT[nclusters][3][2000]/F");
    fTree->Branch("cluhit_endT",cluhit_endT,"cluhit_endT[nclusters][3][2000]/F");
    fTree->Branch("cluhit_starttick",cluhit_starttick,"cluhit_starttick[nclusters][3][2000]/F");
    fTree->Branch("cluhit_endtick",cluhit_endtick,"cluhit_endtick[nclusters][3][2000]/F");
    fTree->Branch("cluhit_RMS",cluhit_RMS,"cluhit_RMS[nclusters][3][2000]/F");
    fTree->Branch("cluhit_ID",cluhit_ID,"cluhit_ID[nclusters][3][2000]/I");
  }  

  fTree->Branch("ntracks",&ntracks,"ntracks/S");
  fTree->Branch("trkId", trkId, "trkId[ntracks]/S");
  fTree->Branch("trkke",trkke,"trkke[ntracks][3]/F");
  fTree->Branch("trkrange",trkrange,"trkrange[ntracks][3]/F");
  fTree->Branch("trkidtruth",trkidtruth,"trkidtruth[ntracks][3]/I");
  fTree->Branch("trkpdgtruth",trkpdgtruth,"trkpdgtruth[ntracks][3]/I");
  fTree->Branch("trkpitchc",trkpitchc,"trkpitchc[ntracks][3]/F");
  fTree->Branch("ntrkhits",ntrkhits,"ntrkhits[ntracks][3]/S");
  if (fSaveCaloInfo){
    fTree->Branch("trkdedx",trkdedx,"trkdedx[ntracks][3][2000]/F");
    fTree->Branch("trkdqdx",trkdqdx,"trkdqdx[ntracks][3][2000]/F");
    fTree->Branch("trkxyz",trkxyz,"trkxyz[ntracks][3][2000][3]/F");  
    fTree->Branch("trkresrg",trkresrg,"trkresrg[ntracks][3][2000]/F");
  }
  fTree->Branch("trkstartx", trkstartx, "trkstartx[ntracks]/F");
  fTree->Branch("trkstarty", trkstarty, "trkstarty[ntracks]/F");
  fTree->Branch("trkstartz", trkstartz, "trkstartz[ntracks]/F");
  fTree->Branch("trkendx", trkendx, "trkendx[ntracks]/F");
  fTree->Branch("trkendy", trkendy, "trkendy[ntracks]/F");
  fTree->Branch("trkendz", trkendz, "trkendz[ntracks]/F");
  fTree->Branch("trktheta", trktheta, "trktheta[ntracks]/F");
  fTree->Branch("trkphi", trkphi, "trkphi[ntracks]/F");
  fTree->Branch("trkstartdcosx", trkstartdcosx,"trkstartdcosx[ntracks]/F");
  fTree->Branch("trkstartdcosy", trkstartdcosy,"trkstartdcosy[ntracks]/F");
  fTree->Branch("trkstartdcosz", trkstartdcosz,"trkstartdcosz[ntracks]/F");
  fTree->Branch("trkenddcosx", trkenddcosx, "trkenddcosx[ntracks]/F");
  fTree->Branch("trkenddcosy", trkenddcosy, "trkenddcosy[ntracks]/F");
  fTree->Branch("trkenddcosz", trkenddcosz, "trkenddcosz[ntracks]/F");
  fTree->Branch("trkthetaxz", trkthetaxz, "trkthetaxz[ntracks]/F");
  fTree->Branch("trkthetayz", trkthetayz, "trkthetayz[ntracks]/F");
  fTree->Branch("trkmom", trkmom, "trkmom[ntracks]/F");
  fTree->Branch("trkmomrange", trkmomrange, "trkmomrange[ntracks]/F");
  fTree->Branch("trkmommschi2", trkmommschi2, "trkmommschi2[ntracks]/F");
  fTree->Branch("trkmommsllhd", trkmommsllhd, "trkmommsllhd[ntracks]/F");  
  fTree->Branch("trklen", trklen, "trklen[ntracks]/F");   
  fTree->Branch("trkpidpdg", trkpidpdg, "trkpidpdg[ntracks][3]/F");  
  fTree->Branch("trkpidchi", trkpidchi, "trkpidchi[ntracks][3]/F");  
  fTree->Branch("trkpidpida", trkpidpida, "trkpidpida[ntracks][3]/F");  
  fTree->Branch("trkpidbestplane", trkpidbestplane, "trkpidbestplane[ntracks]/F");  

  fTree->Branch("no_primaries",&no_primaries,"no_primaries/I");
  fTree->Branch("geant_list_size",&geant_list_size,"geant_list_size/I");
  fTree->Branch("geant_list_size_in_tpcAV",&geant_list_size_in_tpcAV,"geant_list_size_in_tpcAV/I");  
  fTree->Branch("pdg",pdg,"pdg[geant_list_size]/I");
  fTree->Branch("status",status,"status[geant_list_size]/I");
  fTree->Branch("Mass",Mass,"Mass[geant_list_size]/F");
  fTree->Branch("Eng",Eng,"Eng[geant_list_size]/F");
  fTree->Branch("EndE",EndE,"EndE[geant_list_size]/F");
  fTree->Branch("Px",Px,"Px[geant_list_size]/F");
  fTree->Branch("Py",Py,"Py[geant_list_size]/F");
  fTree->Branch("Pz",Pz,"Pz[geant_list_size]/F");
  fTree->Branch("P",P,"P[geant_list_size]/F");
  fTree->Branch("StartPointx",StartPointx,"StartPointx[geant_list_size]/F");
  fTree->Branch("StartPointy",StartPointy,"StartPointy[geant_list_size]/F");
  fTree->Branch("StartPointz",StartPointz,"StartPointz[geant_list_size]/F");
  fTree->Branch("StartT",StartT,"StartT[geant_list_size]/F");
  fTree->Branch("EndPointx",EndPointx,"EndPointx[geant_list_size]/F");
  fTree->Branch("EndPointy",EndPointy,"EndPointy[geant_list_size]/F");
  fTree->Branch("EndPointz",EndPointz,"EndPointz[geant_list_size]/F");
  fTree->Branch("EndT",EndT,"EndT[geant_list_size]/F");
  fTree->Branch("theta",theta,"theta[geant_list_size]/F");
  fTree->Branch("phi",phi,"phi[geant_list_size]/F");
  fTree->Branch("theta_xz",theta_xz,"theta_xz[geant_list_size]/F");
  fTree->Branch("theta_yz",theta_yz,"theta_yz[geant_list_size]/F");
  fTree->Branch("pathlen",pathlen,"pathlen[geant_list_size]/F");
  fTree->Branch("inTPCActive",inTPCActive,"inTPCActive[geant_list_size]/I");  
  fTree->Branch("StartPointx_tpcAV",StartPointx_tpcAV,"StartPointx_tpcAV[geant_list_size]/F");
  fTree->Branch("StartPointy_tpcAV",StartPointy_tpcAV,"StartPointy_tpcAV[geant_list_size]/F");
  fTree->Branch("StartPointz_tpcAV",StartPointz_tpcAV,"StartPointz_tpcAV[geant_list_size]/F");
  fTree->Branch("EndPointx_tpcAV",EndPointx_tpcAV,"EndPointx_tpcAV[geant_list_size]/F");
  fTree->Branch("EndPointy_tpcAV",EndPointy_tpcAV,"EndPointy_tpcAV[geant_list_size]/F");
  fTree->Branch("EndPointz_tpcAV",EndPointz_tpcAV,"EndPointz_tpcAV[geant_list_size]/F");
  fTree->Branch("NumberDaughters",NumberDaughters,"NumberDaughters[geant_list_size]/I");
  fTree->Branch("Mother",Mother,"Mother[geant_list_size]/I");
  fTree->Branch("TrackId",TrackId,"TrackId[geant_list_size]/I");
  fTree->Branch("process_primary",process_primary,"process_primary[geant_list_size]/I");
  fTree->Branch("processname", processname);
}

void microboone::Diffusion::beginSubRun(const art::SubRun& sr)
{
  art::Handle< sumdata::POTSummary > potListHandle;
  
  if(sr.getByLabel(fPOTModuleLabel,potListHandle))
    pot=potListHandle->totpot;
  else
    pot=0.;  
}

void microboone::Diffusion::analyze(const art::Event& evt)
{  
  ResetVars();
  
  bool isMC = !evt.isRealData();

  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);

  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  std::vector<art::Ptr<recob::Cluster> > clusterlist;
  if (fSaveClusterInfo){
      if (evt.getByLabel(fClusterModuleLabel,clusterListHandle))
      art::fill_ptr_vector(clusterlist, clusterListHandle);
  }      

  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (fSaveTrackInfo){
    if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
       art::fill_ptr_vector(tracklist, trackListHandle);
  }     
    
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if (isMC){
     if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
        art::fill_ptr_vector(mclist, mctruthListHandle);
  }  

  //services
  art::ServiceHandle<geo::Geometry> geom;  
  art::ServiceHandle<cheat::BackTracker> bt;
  art::ServiceHandle<util::DetectorProperties> detprop;
  art::ServiceHandle<util::LArProperties> LArProp;
  
  //associations
  if (fSaveTrackInfo){
    art::FindManyP<recob::Track>      fmtk(hitListHandle, evt, fTrackModuleLabel);
    art::FindManyP<recob::Hit>        fmht(trackListHandle, evt, fTrackModuleLabel);
    art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
  }
  
  std::vector<const sim::SimChannel*> fSimChannels;
  if (isMC)
     evt.getView(fLArG4ModuleLabel, fSimChannels);

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

  if (evt.isRealData()){
    isdata = 1;
  }
  else isdata = 0;
  
  //hit information
  const size_t NHits     = hitlist.size(); // number of hits
  no_hits=(int) NHits;
  if (NHits > kMaxHits){
    // got this error? consider increasing kMaxHits
    mf::LogError("Diffusion:limits") << "event has " << NHits
    << " hits, only kMaxHits=" << kMaxHits << " stored in tree";
  }
  
  for (size_t i = 0; i < NHits && i < kMaxHits ; ++i){//loop over hits
    hit_channel[i] = hitlist[i]->Channel();
    hit_plane[i]   = hitlist[i]->WireID().Plane;
    hit_wire[i]    = hitlist[i]->WireID().Wire;
    hit_peakT[i]   = hitlist[i]->PeakTime();
    hit_charge[i]  = hitlist[i]->Integral();
    hit_ph[i]  = hitlist[i]->PeakAmplitude();
    hit_startT[i] = hitlist[i]->PeakTimeMinusRMS();
    hit_endT[i] = hitlist[i]->PeakTimePlusRMS();
    hit_starttick[i] = hitlist[i]->StartTick();
    hit_endtick[i] = hitlist[i]->EndTick();
    hit_RMS[i] = hitlist[i]->RMS();	
    //give a unique id to the hit to map to the full hit collection
    if (hitlist[i]->Channel()>=0 || hitlist[i]->PeakTime()>=0)
        hit_ID[i]  = (hitlist[i]->Channel() % 200000) * 10000 + (int(std::abs(hitlist[i]->PeakTime())) % 10000);
    //std::vector<double> xyz = bt->HitToXYZ(hitlist[i]);
    //when the size of simIDEs is zero, the above function throws an exception
    //and crashes, so check that the simIDEs have non-zero size before 
    //extracting hit true XYZ from simIDEs
    if (isMC){
      std::vector<sim::IDE> ides;
      bt->HitToSimIDEs(hitlist[i], ides);
      if (ides.size()>0){
         std::vector<double> xyz = bt->SimIDEsToXYZ(ides);
         hit_trueX[i] = xyz[0];
      }
    }   
    	    
    //Hit to CalWire information	    
    art::FindManyP<recob::Wire> fmwire(hitListHandle,evt,fHitsModuleLabel);
     
    int dataSize = fmwire.at(i)[0]->NSignal();
    int t0 = 0;
    int t1 = dataSize-1;
    std::vector<float> signal(fmwire.at(i)[0]->Signal());
    
    wire_pk[i] = -1.0;
    wire_t[i] = -1.0;
    for (int j = t0; j<=t1; ++j){
       if (signal[j]>wire_pk[i]){
    	  wire_pk[i] = signal[j];
    	  wire_t[i]  = j;
       }   
    }

    wire_ch[i] = 0.0;
    double mean_t = 0.0;
    double mean_t2 = 0.0;
    for (int j = t0; j<=t1; ++j){
       if (signal[j]>=0.1*wire_pk[i]){
    	   wire_ch[i] += signal[j];
    	   mean_t  += double(j)*(signal[j]);
    	   mean_t2 += double(j)*double(j)*(signal[j]);
       }
       
    }
    mean_t/=wire_ch[i];
    mean_t2/=wire_ch[i];
    wire_rms[i] = sqrt(mean_t2-mean_t*mean_t);
    	
    //Hit to RawDigit information    	 
    art::FindManyP<raw::RawDigit> fmrd(hitListHandle,evt,fHitsModuleLabel);

    dataSize = fmrd.at(i)[0]->Samples();
    short ped = fmrd.at(i)[0]->GetPedestal();
    std::vector<short> rawadc(dataSize);
    raw::Uncompress(fmrd.at(i)[0]->ADCs(), rawadc, fmrd.at(i)[0]->Compression());
    t0 = 0;
    t1 = dataSize-1;
    rawD_pk[i] = -1;
    rawD_t[i] = -1;
    for (int j = t0; j<=t1; ++j){
      if (rawadc[j]-ped>rawD_pk[i]){
    	rawD_pk[i] = rawadc[j]-ped;
    	rawD_t[i] = j;
      }
    }
    rawD_ch[i] = 0;
    rawD_fwhh[i] = 0;
    mean_t = 0.0;
    mean_t2 = 0.0;
    for (int j = t0; j<=t1; ++j){
      if (rawadc[j]-ped>=0.5*rawD_pk[i]){
    	++rawD_fwhh[i];
      }
      if (rawadc[j]-ped>=0.1*rawD_pk[i]){
    	rawD_ch[i] += rawadc[j]-ped;
    	mean_t += (double)j*(rawadc[j]-ped);
    	mean_t2 += (double)j*(double)j*(rawadc[j]-ped);
      }
    }
    mean_t/=rawD_ch[i];
    mean_t2/=rawD_ch[i];
    rawD_rms[i] = sqrt(mean_t2-mean_t*mean_t);
    
    //Hit to SimChannel information
    if (isMC){
      hit_nelec[i] = 0;
      hit_energy[i] = 0;
      const sim::SimChannel* chan = 0;
      for(size_t sc = 0; sc < fSimChannels.size(); ++sc){
        if(fSimChannels[sc]->Channel() == hitlist[i]->Channel()) 
    	   chan = fSimChannels[sc];
      }     
      if (chan){
        const std::map<unsigned short, std::vector<sim::IDE> >& tdcidemap = chan->TDCIDEMap();
        int k=-1;
        double elec[tdcidemap.size()];
        int tdc[tdcidemap.size()];
        for(auto mapitr = tdcidemap.begin(); mapitr != tdcidemap.end(); mapitr++){
    	   k++;
    	   tdc[k]=(*mapitr).first;
    	   const std::vector<sim::IDE> idevec = (*mapitr).second;
   	   double nelec=0;
    	   for(size_t iv = 0; iv < idevec.size(); ++iv){
    	      nelec += idevec[iv].numElectrons;
	      hit_nelec[i] += idevec[iv].numElectrons;
	      hit_energy[i] += idevec[iv].energy;
    	   }
    	   elec[k] = nelec;
        }
        sim_pk[i] = -1;
        sim_tdc[i] = -1;
        for(unsigned int f=0;f<tdcidemap.size();f++){
    	  if (elec[f]>sim_pk[i]){
    	    sim_pk[i] = elec[f];
    	    sim_tdc[i] = tdc[f]; 
    	  }	    
        }
        sim_ch[i] = 0;
    	mean_t = 0;
    	mean_t2 = 0;
    	for (unsigned int f = 0; f<tdcidemap.size();f++){
    	   if (elec[f]>=0.1*sim_pk[i]){
    	       sim_ch[i]+= elec[f];
    	       mean_t+= double(tdc[f])*elec[f];
    	       mean_t2+= double(tdc[f])*double(tdc[f])*elec[f];
    	   }
    	}
    	mean_t/=sim_ch[i];
    	mean_t2/=sim_ch[i];
    	sim_rms[i] = sqrt(mean_t2-mean_t*mean_t);
      }
    }
  }     

  /*if (evt.getByLabel(fHitsModuleLabel,hitListHandle)){
    //Find tracks associated with hits
    art::FindManyP<recob::Track> fmtk(hitListHandle,evt,fTrackModuleLabel);
    for (size_t i = 0; i < NHits && i < kMaxHits ; ++i){//loop over hits
      if (fmtk.isValid()){
  	if (fmtk.at(i).size()!=0){
  	  hit_trkid[i] = fmtk.at(i)[0]->ID();
  	}
  	else
  	  hit_trkid[i] = -1;
      }
    }
  }*/
  
  //doesn't work for whatever reason!!! will resolve later!
  /*if (evt.getByLabel(fHitsModuleLabel,hitListHandle)){
    //Find clusters associated with hits
    art::FindManyP<recob::Cluster> fmcl(hitListHandle,evt,fClusterModuleLabel);
    for (size_t i = 0; i < NHits && i < kMaxHits ; ++i){//loop over hits
      if (fmcl.isValid()){
        std::cout<<"\n"<<i<<","<<fmcl.at(i).size();
  	if (fmcl.at(i).size()!=0){
  	  hit_clusterid[i] = fmcl.at(i)[0]->ID();
  	}
  	else
  	  hit_clusterid[i] = -1;
      }
    }
  }*/
  
  if (fSaveClusterInfo){
     size_t NClusters = clusterlist.size();
     nclusters = (int) NClusters;
     for(unsigned int ic=0; ic<NClusters;++ic){//loop over clusters
         art::Ptr<recob::Cluster> clusterholder(clusterListHandle, ic);
	 const recob::Cluster& cluster = *clusterholder;
	 clusterId[ic] = cluster.ID();
	 clusterView[ic] = cluster.View();
	 cluster_isValid[ic] = cluster.isValid();
	 cluster_StartCharge[ic] = cluster.StartCharge();
	 cluster_StartAngle[ic] = cluster.StartAngle();
	 cluster_EndCharge[ic] = cluster.EndCharge();
	 cluster_EndAngle[ic] = cluster.EndAngle();
	 cluster_Integral[ic] = cluster.Integral();
	 cluster_IntegralAverage[ic] = cluster.IntegralAverage();
	 cluster_SummedADC[ic] = cluster.SummedADC();
	 cluster_SummedADCaverage[ic] = cluster.SummedADCaverage();
	 cluster_MultipleHitDensity[ic] = cluster.MultipleHitDensity();
	 cluster_Width[ic] = cluster.Width();
	 cluster_NHits[ic] = cluster.NHits();
	 cluster_StartWire[ic] = cluster.StartWire();
	 cluster_StartTick[ic] = cluster.StartTick();
	 cluster_EndWire[ic] = cluster.EndWire();
	 cluster_EndTick[ic] = cluster.EndTick();
	 
	 if (cluster.View() == 0) ncluhits[ic][0] = cluster.NHits();
         if (cluster.View() == 1) ncluhits[ic][1] = cluster.NHits();
         if (cluster.View() == 2) ncluhits[ic][2] = cluster.NHits();
	 	 
	 if (cluster.NHits() > kMaxClusterHits){
     	   // if you get this error, you'll have to increase kMaxClusterHits
     	   mf::LogError("Diffusion:limits")
     	     << "the " << fClusterModuleLabel << " Cluster #" << ic
     	     << " has " << cluster.NHits() << " hits on plane #" << cluster.View()
     	     <<", only "
     	     << kMaxClusterHits << " stored in tree";
     	 }
	 
	 if (fSaveClusterHitInfo){	 
	 art::FindManyP<recob::Hit> fmhc(clusterListHandle, evt, fClusterModuleLabel);
	 if (fmhc.isValid()){
	     std::vector< art::Ptr<recob::Hit> > cluhits = fmhc.at(ic);
	     //note:cluhits.size() is same as cluster.NHits() above!
	     for (Int_t pl=0;pl<3;pl++){
	        for (Int_t iht = 0; iht<ncluhits[ic][pl]; ++iht){
		  if (clusterView[ic] == pl){
		    //std::cout<<ncluhits[ic][pl]<<","<<iht;
	            cluhit_channel[ic][pl][iht]	= cluhits[iht]->Channel();
	            cluhit_plane[ic][pl][iht]	= cluhits[iht]->WireID().Plane;
	            cluhit_wire[ic][pl][iht]	= cluhits[iht]->WireID().Wire;
	            cluhit_peakT[ic][pl][iht]	= cluhits[iht]->PeakTime();
	            cluhit_charge[ic][pl][iht]	= cluhits[iht]->Integral();
	            cluhit_ph[ic][pl][iht]	= cluhits[iht]->PeakAmplitude();
	            cluhit_startT[ic][pl][iht]	= cluhits[iht]->PeakTimeMinusRMS();
	            cluhit_endT[ic][pl][iht]	= cluhits[iht]->PeakTimePlusRMS();
	            cluhit_starttick[ic][pl][iht] = cluhits[iht]->StartTick();
	            cluhit_endtick[ic][pl][iht]	= cluhits[iht]->EndTick();
	            cluhit_RMS[ic][pl][iht]	= cluhits[iht]->RMS(); 
		    //give a unique id to the hit to map to the full hit collection
		    if (cluhits[iht]->Channel() >=0 || cluhits[iht]->PeakTime()>=0)
		       cluhit_ID[ic][pl][iht]  = (cluhits[iht]->Channel() % 200000) * 10000 + (int(std::abs(cluhits[iht]->PeakTime())) % 10000);
		 }      
	       } //end loop over cluster hits
	    }//end loop over planes
	  }//end fmhc.isValid()  
	 }//end fSaveClusterHitInfo
	  
     }//end loop over clusters
  }//end fSaveClusterInfo
	 
  //Track information
  if (fSaveTrackInfo){
  size_t NTracks = tracklist.size(); 
  ntracks= (int) NTracks;
 
  //call the track momentum algorithm that gives you momentum based on track range
  trkf::TrackMomentumCalculator trkm;

  for(unsigned int i=0; i<NTracks;++i){//loop over tracks
    art::Ptr<recob::Track> ptrack(trackListHandle, i);
    const recob::Track& track = *ptrack;   
    
    TVector3 pos, dir_start, dir_end, end;        

    double tlen = 0., mom = 0.;
    int TrackID = -1; 
    int ntraj = 0;

    //we need to use Bezier methods for Bezier tracks
    if (fTrackModuleLabel.find("beziertracker")!=std::string::npos) {
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
         TrackID = i;  //bezier has some screwed up track IDs
       }
     }  	 
     else {   //use the normal methods for other kinds of tracks
        ntraj = track.NumberTrajectoryPoints();
        if (ntraj > 0) {
     	  pos	   = track.Vertex();
     	  dir_start = track.VertexDirection();
     	  dir_end   = track.EndDirection();
     	  end	   = track.End();

     	  tlen	     = length(track);
     	  if(track.NumberFitMomentum() > 0)
     	     mom = track.VertexMomentum();
     	  //fill non-bezier-track reco branches
     	  TrackID = track.ID();
        }
     }
     
     if (ntraj > 0) {
       double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
       double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());
       trkId[i]		   = TrackID;
       trkstartx[i]  	   = pos.X();
       trkstarty[i]  	   = pos.Y();
       trkstartz[i]  	   = pos.Z();
       trkendx[i]	   = end.X();
       trkendy[i]	   = end.Y();
       trkendz[i]	   = end.Z();
       trktheta[i]	   = dir_start.Theta();
       trkphi[i]	   = dir_start.Phi();
       trkstartdcosx[i]	   = dir_start.X();
       trkstartdcosy[i]	   = dir_start.Y();
       trkstartdcosz[i]	   = dir_start.Z();
       trkenddcosx[i]	   = dir_end.X();
       trkenddcosy[i]	   = dir_end.Y();
       trkenddcosz[i]	   = dir_end.Z();
       trkthetaxz[i] 	   = theta_xz;
       trkthetayz[i] 	   = theta_yz;
       trkmom[i]	   = mom;
       trklen[i]	   = tlen;
       trkmomrange[i]	   = trkm.GetTrackMomentum(tlen,13);
       trkmommschi2[i]	   = trkm.GetMomentumMultiScatterChi2(ptrack);
       trkmommsllhd[i]	   = trkm.GetMomentumMultiScatterLLHD(ptrack);
     } // if we have trajectory
     
     // find particle ID info
     art::FindMany<anab::ParticleID> fmpid(trackListHandle, evt, fParticleIDModuleLabel);
     if(fmpid.isValid()) {
       std::vector<const anab::ParticleID*> pids = fmpid.at(i);
       for (size_t ipid = 0; ipid < pids.size(); ++ipid){
     	 if (!pids[ipid]->PlaneID().isValid) continue;
     	 int planenum = pids[ipid]->PlaneID().Plane;
     	 if (planenum<0||planenum>2) continue;
     	 trkpidpdg[i][planenum]  = pids[ipid]->Pdg();
     	 trkpidchi[i][planenum]  = pids[ipid]->MinChi2();
     	 trkpidpida[i][planenum] = pids[ipid]->PIDA();
       }
     } // fmpid.isValid()
     
     art::FindMany<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
     //find Calorimetry info
     if (fmcal.isValid()){
       std::vector<const anab::Calorimetry*> calos = fmcal.at(i);
       if (calos.size() > 3) {
     	 // if you get this message, there is probably a bug somewhere since
     	 // the calorimetry planes should be 3.
     	 mf::LogError("Diffusion:limits")
     	   << "the " << fTrackModuleLabel << " track #" << i
     	   << " has " << calos.size() << " planes for calorimetry , only 3"
     	   << " stored in tree";
       }
       for (size_t ical = 0; ical<calos.size(); ++ical){
     	 if (!calos[ical]) continue;
     	 if (!calos[ical]->PlaneID().isValid) continue;
     	 int planenum = calos[ical]->PlaneID().Plane;
     	 if (planenum<0||planenum>2) continue;
     	 trkke[i][planenum]    = calos[ical]->KineticEnergy();
     	 trkrange[i][planenum] = calos[ical]->Range();
     	 //For now make the second argument as 13 for muons. 
     	 trkpitchc[i][planenum]= calos[ical] -> TrkPitchC();
     	 const size_t NHits = calos[ical] -> dEdx().size();
     	 ntrkhits[i][planenum] = (int) NHits;
     	 if (NHits > kMaxTrackHits){
     	   // if you get this error, you'll have to increase kMaxTrackHits
     	   mf::LogError("Diffusion:limits")
     	     << "the " << fTrackModuleLabel << " track #" << i
     	     << " has " << NHits << " hits on calorimetry plane #" << planenum
     	     <<", only "
     	     << kMaxTrackHits << " stored in tree";
     	 }
	 if (fSaveCaloInfo){
     	   for(size_t iHit = 0; iHit < NHits && iHit < kMaxTrackHits; ++iHit) {
     	     trkdedx[i][planenum][iHit]  = (calos[ical] -> dEdx())[iHit];
     	     trkdqdx[i][planenum][iHit]  = (calos[ical] -> dQdx())[iHit];
     	     trkresrg[i][planenum][iHit] = (calos[ical] -> ResidualRange())[iHit];
     	     const auto& TrkPos = (calos[ical] -> XYZ())[iHit];
     	     auto& TrkXYZ = trkxyz[i][planenum][iHit];
     	     TrkXYZ[0] = TrkPos.X();
     	     TrkXYZ[1] = TrkPos.Y();
     	     TrkXYZ[2] = TrkPos.Z();
     	   } // for track hits
	 }// if (fSaveCaloInfo) 
       } // for calorimetry info
       
       if(ntrkhits[i][0] > ntrkhits[i][1] && ntrkhits[i][0] > ntrkhits[i][2]) trkpidbestplane[i] = 0;
       else if(ntrkhits[i][1] > ntrkhits[i][0] && ntrkhits[i][1] > ntrkhits[i][2]) trkpidbestplane[i] = 1;
       else if(ntrkhits[i][2] > ntrkhits[i][0] && ntrkhits[i][2] > ntrkhits[i][1]) trkpidbestplane[i] = 2;
       else if(ntrkhits[i][2] == ntrkhits[i][0] && ntrkhits[i][2] > ntrkhits[i][1]) trkpidbestplane[i] = 2;
       else if(ntrkhits[i][2] == ntrkhits[i][1] && ntrkhits[i][2] > ntrkhits[i][0]) trkpidbestplane[i] = 2;
       else if(ntrkhits[i][1] == ntrkhits[i][0] && ntrkhits[i][1] > ntrkhits[i][2]) trkpidbestplane[i] = 0;
       else if(ntrkhits[i][1] == ntrkhits[i][0] && ntrkhits[i][1] == ntrkhits[i][2]) trkpidbestplane[i] = 2;
     } // if has calorimetry info
     
     //track truth information
     if (isMC){
       //get the hits on each plane
       art::FindManyP<recob::Hit>      fmht(trackListHandle, evt, fTrackModuleLabel);
       std::vector< art::Ptr<recob::Hit> > allHits = fmht.at(i);
       std::vector< art::Ptr<recob::Hit> > hits[kNplanes];
       for(size_t ah = 0; ah < allHits.size(); ++ah){
     	 if (allHits[ah]->WireID().Plane <  3){
     	   hits[allHits[ah]->WireID().Plane].push_back(allHits[ah]);
     	 }
       }       
       for (size_t ipl = 0; ipl < 3; ++ipl){
     	 double maxe = 0;
	 float purity;
     	 HitsPurity(hits[ipl],trkidtruth[i][ipl],purity,maxe);
     	 if (trkidtruth[i][ipl]>0){
     	   const simb::MCParticle *particle = bt->TrackIDToParticle(trkidtruth[i][ipl]);
     	   trkpdgtruth[i][ipl] = particle->PdgCode();
     	 }
       }
     }     
  }//end loop over tracks
  }//if fSaveTrackInfo

  //GEANT information
  if (!isdata){//is MC
    Int_t mcevts_truth = mclist.size();
    if (mcevts_truth){//at least one mc record 
      //GEANT particles information
      const sim::ParticleList& plist = bt->ParticleList();    
      
      std::string pri("primary");
      int primary=0;
      int active = 0;
      int geant_particle=0;
      sim::ParticleList::const_iterator itPart = plist.begin(),
        pend = plist.end(); // iterator to pairs (track id, particle)
        	
      for(size_t iPart = 0; (iPart < plist.size()) && (itPart != pend); ++iPart){
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

         if (pPart->E()>fG4minE||isPrimary){
          process_primary[iPart] = int(isPrimary);
          processname[iPart]= pPart->Process();
          Mother[iPart]=pPart->Mother();
          TrackId[iPart]=TrackID;
          pdg[iPart]=pPart->PdgCode();
          status[iPart] = pPart->StatusCode();
          Eng[iPart]=pPart->E();
          EndE[iPart]=pPart->EndE();
          Mass[iPart]=pPart->Mass();
          Px[iPart]=pPart->Px();
          Py[iPart]=pPart->Py();
          Pz[iPart]=pPart->Pz();
          P[iPart]=pPart->Momentum().Vect().Mag();
          StartPointx[iPart]=pPart->Vx();
          StartPointy[iPart]=pPart->Vy();
          StartPointz[iPart]=pPart->Vz();
          StartT[iPart] = pPart->T();
          EndPointx[iPart]=pPart->EndPosition()[0];
          EndPointy[iPart]=pPart->EndPosition()[1];
          EndPointz[iPart]=pPart->EndPosition()[2];
          EndT[iPart] = pPart->EndT();
          theta[iPart] = pPart->Momentum().Theta();
          phi[iPart] = pPart->Momentum().Phi();
          theta_xz[iPart] = std::atan2(pPart->Px(), pPart->Pz());
          theta_yz[iPart] = std::atan2(pPart->Py(), pPart->Pz());
          pathlen[iPart]  = plen;
          NumberDaughters[iPart]=pPart->NumberDaughters();
          inTPCActive[iPart] = int(isActive);
          if (isActive){	
            StartPointx_tpcAV[iPart] = mcstart.X();
            StartPointy_tpcAV[iPart] = mcstart.Y();
            StartPointz_tpcAV[iPart] = mcstart.Z();
            EndPointx_tpcAV[iPart] = mcend.X();
            EndPointy_tpcAV[iPart] = mcend.Y();
            EndPointz_tpcAV[iPart] = mcend.Z();
          }
	} 		     
     }//plist loop    
     
     geant_list_size_in_tpcAV = active;
     no_primaries = primary;
     geant_list_size = geant_particle;
   }//if (mcevts_truth)
  }//if (!isdata) 
  
  taulife = LArProp->ElectronLifetime();
  fTree->Fill();
}

void microboone::Diffusion::HitsPurity(std::vector< art::Ptr<recob::Hit> > const& hits, Int_t& trackid, Float_t& purity, double& maxe){

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

// Length of reconstructed track, trajectory by trajectory.
double microboone::Diffusion::length(const recob::Track& track)
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

// Length of MC particle, tracjectory by tracjectory.
double microboone::Diffusion::length(const simb::MCParticle& part, TVector3& start, TVector3& end)
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
      // Doing some manual shifting to account for an interaction not occuring with the beam dump
      // we will reconstruct an x distance different from where the particle actually passed to to the time
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

void microboone::Diffusion::ResetVars(){

  run = -99999;
  subrun = -99999;
  event = -99999;
  evttime = -99999;
  beamtime = -99999;
  isdata = -99;
  taulife = -99999;

  no_hits = 0;
  for (int i = 0; i<kMaxHits; ++i){
    hit_plane[i] = -9999;
    hit_wire[i] = -9999;
    hit_channel[i] = -9999;
    hit_peakT[i] = -99999;
    hit_charge[i] = -99999;
    hit_ph[i] = -99999;
    hit_startT[i] = -99999.; 
    hit_endT[i] = -99999.; 
    hit_starttick[i] = -99999.;
    hit_endtick[i] = -99999.; 
    hit_RMS[i] = -99999.;
    hit_nelec[i] = -99999.;
    hit_energy[i] = -99999.;
    hit_trueX[i] = -99999.;
    hit_trkid[i] = -9999; 
    hit_clusterid[i] = -9999;   
    hit_ID[i] = -9999; 
    //
    rawD_pk[i] = -9999;
    rawD_t[i] = -9999;
    rawD_ch[i] = -9999;
    rawD_fwhh[i] = -9999;
    rawD_rms[i] = -99999.;
    wire_pk[i] = -99999.;
    wire_t[i] = -9999;
    wire_ch[i] = -99999.;
    wire_rms[i] = -99999.;
    sim_pk[i] = -99999.;
    sim_tdc[i] = -9999;
    sim_ch[i] = -99999.;
    sim_rms[i] = -99999.;
          
  }

  nclusters=0;
  for (int i = 0; i<kMaxCluster; ++i){
    clusterId[i] = -9999;
    clusterView[i] = -9999;
    cluster_isValid[i] = -1;
    cluster_StartCharge[i] = -99999.;
    cluster_StartAngle[i] = -99999.;
    cluster_EndCharge[i] = -99999.;
    cluster_EndAngle[i] = -99999.;
    cluster_Integral[i] = -99999.;
    cluster_IntegralAverage[i] = -99999.;
    cluster_SummedADC[i] = -99999.;
    cluster_SummedADCaverage[i] = -99999.;
    cluster_MultipleHitDensity[i] = -99999.;
    cluster_Width[i] = -99999.;
    cluster_NHits[i] = -9999;
    cluster_StartWire[i] = -9999;
    cluster_StartTick[i] = -9999;
    cluster_EndWire[i] = -9999;
    cluster_EndTick[i] = -9999;
    
    for (int j = 0; j < kNplanes; ++j){
      ncluhits[i][j] = -9999;
      if (fSaveClusterHitInfo){
        for(int k = 0; k < kMaxClusterHits; k++) {
    	  cluhit_plane[i][j][k]   = -9999;
   	  cluhit_wire[i][j][k]    = -9999;
   	  cluhit_channel[i][j][k] = -9999;
   	  cluhit_peakT[i][j][k]   = -99999;
          cluhit_charge[i][j][k]  = -99999;
          cluhit_ph[i][j][k]      = -99999;
    	  cluhit_startT[i][j][k]  = -99999.; 
    	  cluhit_endT[i][j][k]    = -99999.; 
    	  cluhit_starttick[i][j][k] = -99999.;
    	  cluhit_endtick[i][j][k]   = -99999.; 
    	  cluhit_RMS[i][j][k]       = -99999.;
	  cluhit_ID[i][j][k]       = -9999;
	}  
      }	
    } 
  }  

  ntracks = 0;
  for (int i = 0; i < kMaxTrack; ++i){
    trkId[i] = -9999;
    trkstartx[i] = -99999.;   
    trkstarty[i] = -99999.;   
    trkstartz[i] = -99999.;   
    trkendx[i] = -99999.;     
    trkendy[i] = -99999.;     
    trkendz[i] = -99999.;     
    trktheta[i] = -99999.;    
    trkphi[i] = -99999.; 
    trkstartdcosx[i] = -99999.;     
    trkstartdcosy[i] = -99999.;     
    trkstartdcosz[i] = -99999.;  
    trkenddcosx[i] = -99999.;
    trkenddcosy[i] = -99999.;
    trkenddcosz[i] = -99999.;   
    trkthetaxz[i] = -99999.; 
    trkthetayz[i] = -99999.; 
    trkmom[i] = -99999.;      
    trkmomrange[i] = -99999.;      
    trkmommschi2[i] = -99999.;      
    trkmommsllhd[i] = -99999.;      
    trklen[i] = -99999.;
    trkpidbestplane[i]= -1;
    
    for (int j = 0; j < kNplanes; ++j){
      trkke[i][j] = -99999.;
      trkrange[i][j] = -99999.;
      trkidtruth[i][j] = -99999;
      trkpdgtruth[i][j] = -99999;
      trkpitchc[i][j] = -99999.;
      ntrkhits[i][j] = -9999;
      trkpidpdg[i][j] = -1;
      trkpidchi[i][j] = -99999.;
      trkpidpida[i][j] = -99999.;

      if (fSaveCaloInfo){
        for(int k = 0; k < kMaxTrackHits; k++) {
          trkdedx[i][j][k]  = 0.;
          trkdqdx[i][j][k]  = 0.;
          trkresrg[i][j][k] = 0.;
	  for(int l=0;l < 3; l++)
	    trkxyz[i][j][k][l] = 0.;
        }
      }	
    }
  }

  no_primaries = 0;
  geant_list_size=0;
  geant_list_size_in_tpcAV = 0;
  for (int i = 0; i<kMaxPrimaries; ++i){
    pdg[i] = -99999;
    status[i] = -99999;
    Mass[i] = -99999.;
    Eng[i] = -99999.;
    EndE[i] = -99999.;
    Px[i] = -99999.;
    Py[i] = -99999.;
    Pz[i] = -99999.;
    P[i] = -99999.;
    StartPointx[i] = -99999.;
    StartPointy[i] = -99999.;
    StartPointz[i] = -99999.;
    StartT[i] = -99999.;
    EndT[i] = -99999.;    
    EndPointx[i] = -99999.;
    EndPointy[i] = -99999.;
    EndPointz[i] = -99999.;
    EndT[i] = -99999.;
    theta[i] = -99999.;
    phi[i] = -99999.;
    theta_xz[i] = -99999.;
    theta_yz[i] = -99999.;
    pathlen[i] = -99999.;
    inTPCActive[i] = -99999;
    StartPointx_tpcAV[i] = -99999.;
    StartPointy_tpcAV[i] = -99999.;
    StartPointz_tpcAV[i] = -99999.;
    EndPointx_tpcAV[i] = -99999.;
    EndPointy_tpcAV[i] = -99999.;
    EndPointz_tpcAV[i] = -99999.;  
    NumberDaughters[i] = -99999;
    Mother[i] = -99999;
    TrackId[i] = -99999;
    process_primary[i] = -99999;
    processname[i] = "noname";
  }    
}

namespace microboone{

  DEFINE_ART_MODULE(Diffusion)
  
} 

#endif


