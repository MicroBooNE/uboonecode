////////////////////////////////////////////////
// Lifetime module to be used for nearline monitoring
// Author: Sowjanya Gollapinni, first commit: July 11, 2016
// 
// The module selects tracks that cross anode and cathode
// and applies some angular cuts and produces dQ/ds vs drifttime
// scatter plot. The actual lifetime analysis is done in a 
// post processing script to extract QA/QC value.
// Associated fcl file is lifetime.fcl
///////////////////////////////////////////////

#ifndef Lifetime_Module
#define Lifetime_Module

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/RecoBase/Track.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/Utilities/AssociationUtil.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLine.h"
#include "TAxis.h"
#include <vector>
#include <fstream>
#include "TPaveStats.h"

const int nbin = 22; //split the total drift time into 22 bins
const double binsize = 100; //us
const int kMaxTrack  = 1000;  //maximum number of tracks
const int kNplanes   = 3;     //number of wire planes

using namespace std;

namespace {

// Local functions.

//========================================================================
// Length of reconstructed track, trajectory by trajectory.
double length(const recob::Track& track)
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

} // end namespace

//========================================================================

namespace microboone{

class Lifetime : public art::EDAnalyzer {
public:

    explicit Lifetime(fhicl::ParameterSet const& pset);
    virtual ~Lifetime();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);

    void reconfigure(fhicl::ParameterSet const& pset);
    void reset();

private:

    // the parameters we'll read from the .fcl
    std::string fTrackModuleLabel;
    std::string fCalorimetryModuleLabel;
    bool  fSaveTrackInfo;
     
    TTree *fEventTree;
    TH2D *dqdstime;    
    TH1D *dqdx[nbin];
 
    // Event 
    Int_t    run;                  //run number
    Int_t    subrun;               //subrun number
    Int_t    event;                //event number
    Double_t evttime; 		   //event time in sec
    
    //Track information
    //Track plane data
    Short_t    ntracks;               
    Short_t    ntracks_cross;   
    Short_t    ntracks_selec;                               
    Short_t    ntrkhits_selec[kMaxTrack][kNplanes];
    
    // more track info
    Float_t trklen[kMaxTrack];
    Float_t trklenx[kMaxTrack];
    Float_t trklenx_cross[kMaxTrack];
    Float_t trklenx_selec[kMaxTrack];
    Float_t trkstartx_cross[kMaxTrack];     // starting x position.
    Float_t trkstarty_cross[kMaxTrack];     // starting y position.
    Float_t trkstartz_cross[kMaxTrack];     // starting z position.
    Float_t trkendx_cross[kMaxTrack];	   // ending x position.
    Float_t trkendy_cross[kMaxTrack];	   // ending y position.
    Float_t trkendz_cross[kMaxTrack];	   // ending z position.
    Float_t trktheta_cross[kMaxTrack];	   // theta.
    Float_t trkphi_cross[kMaxTrack];	   // phi.
    Float_t trkstartdcosx_cross[kMaxTrack];
    Float_t trkstartdcosy_cross[kMaxTrack];
    Float_t trkstartdcosz_cross[kMaxTrack];
    Float_t trkenddcosx_cross[kMaxTrack];
    Float_t trkenddcosy_cross[kMaxTrack];
    Float_t trkenddcosz_cross[kMaxTrack];
    Float_t trkthetaxz[kMaxTrack];    // theta_xz.
    Float_t trkthetayz[kMaxTrack];    // theta_xz.    
    Float_t trkthetaxz_cross[kMaxTrack];    // theta_xz.
    Float_t trkthetayz_cross[kMaxTrack];    // theta_yz.
    Float_t trkthetaxz_selec[kMaxTrack];    // theta_xz.
    Float_t trkthetayz_selec[kMaxTrack];    // theta_xz.    
    Float_t trkmom_cross[kMaxTrack];	   // momentum.
    Float_t trklen_cross[kMaxTrack];	   // length.
    
    detinfo::DetectorProperties const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double XDriftVelocity = detprop->DriftVelocity()*1e-3; //cm/ns
    art::ServiceHandle<geo::Geometry> geom;

 }; // class Lifetime


//========================================================================
Lifetime::Lifetime(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
{
    reconfigure(parameterSet);
}
//========================================================================
Lifetime::~Lifetime(){
  //destructor
}
//========================================================================
void Lifetime::reconfigure(fhicl::ParameterSet const& p){

    fTrackModuleLabel         = p.get<std::string>("TrackModuleLabel");
    fCalorimetryModuleLabel   = p.get<std::string>("CalorimetryModuleLabel");
    fSaveTrackInfo            = p.get< bool>("SaveTrackInfo",false);  
 
}
//========================================================================
void Lifetime::beginJob(){
  std::cout<<"job begin..."<<std::endl;

  art::ServiceHandle<art::TFileService> tfs;
  
  dqdstime = tfs->make<TH2D>("dqdstime","",100,0,2200,100,0,1000);
  dqdstime->GetXaxis()->SetTitle("Drift time (#mus)");
  dqdstime->GetYaxis()->SetTitle("dQ/ds (ADC/cm)");
  
  for (int i = 0; i<nbin; ++i){
      dqdx[i] = tfs->make<TH1D>(Form("dqdx_%d",i),Form("dqdx_%d",i),100,0,1500);
      dqdx[i]->GetXaxis()->SetTitle("dq/dx (ADC/cm)");
      dqdx[i]->Sumw2();
    }

    fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");
    fEventTree->Branch("event", &event);
    fEventTree->Branch("run", &run);
    fEventTree->Branch("subrun", &subrun); 
    fEventTree->Branch("evttime", &evttime);
    if (fSaveTrackInfo){
      fEventTree->Branch("ntracks",&ntracks,"ntracks/S");
      fEventTree->Branch("ntracks_cross",&ntracks_cross,"ntracks_cross/S");
      fEventTree->Branch("ntracks_selec",&ntracks_selec,"ntracks_selec/S");	        	  
      fEventTree->Branch("trklen", trklen, "trklen[ntracks]/F");
      fEventTree->Branch("trklenx", trklenx, "trklenx[ntracks]/F");
      fEventTree->Branch("trkthetaxz", trkthetaxz, "trkthetaxz[ntracks]/F");
      fEventTree->Branch("trkthetaxz_cross", trkthetaxz_cross, "trkthetaxz_cross[ntracks_cross]/F");
      fEventTree->Branch("trkthetaxz_selec", trkthetaxz_selec, "trkthetaxz_selec[ntracks_selec]/F");
      fEventTree->Branch("trkthetayz", trkthetayz, "trkthetayz[ntracks]/F");
      fEventTree->Branch("trkthetayz_cross", trkthetayz_cross, "trkthetayz_cross[ntracks_cross]/F");
      fEventTree->Branch("trkthetayz_selec", trkthetayz_selec, "trkthetayz_selec[ntracks_selec]/F");
      fEventTree->Branch("trkstartx_cross", trkstartx_cross, "trkstartx_cross[ntracks_cross]/F");
      fEventTree->Branch("trkstarty_cross", trkstarty_cross, "trkstarty_cross[ntracks_cross]/F");
      fEventTree->Branch("trkstartz_cross", trkstartz_cross, "trkstartz_cross[ntracks_cross]/F");
      fEventTree->Branch("trkendx_cross", trkendx_cross, "trkendx_cross[ntracks_cross]/F");
      fEventTree->Branch("trkendy_cross", trkendy_cross, "trkendy_cross[ntracks_cross]/F");
      fEventTree->Branch("trkendz_cross", trkendz_cross, "trkendz_cross[ntracks_cross]/F");
      fEventTree->Branch("trktheta_cross", trktheta_cross, "trktheta_cross[ntracks_cross]/F");
      fEventTree->Branch("trkphi_cross", trkphi_cross, "trkphi_cross[ntracks_cross]/F");
      fEventTree->Branch("trkstartdcosx_cross", trkstartdcosx_cross,"trkstartdcosx_cross[ntracks_cross]/F");
      fEventTree->Branch("trkstartdcosy_cross", trkstartdcosy_cross,"trkstartdcosy_cross[ntracks_cross]/F");
      fEventTree->Branch("trkstartdcosz_cross", trkstartdcosz_cross,"trkstartdcosz_cross[ntracks_cross]/F");
      fEventTree->Branch("trkenddcosx_cross", trkenddcosx_cross, "trkenddcosx_cross[ntracks_cross]/F");
      fEventTree->Branch("trkenddcosy_cross", trkenddcosy_cross, "trkenddcosy_cross[ntracks_cross]/F");
      fEventTree->Branch("trkenddcosz_cross", trkenddcosz_cross, "trkenddcosz_cross[ntracks_cross]/F");
      fEventTree->Branch("trkmom_cross", trkmom_cross, "trkmom_cross[ntracks_cross]/F");
      fEventTree->Branch("trklen_cross", trklen_cross, "trklen_cross[ntracks_cross]/F");
      fEventTree->Branch("trklenx_cross", trklenx_cross, "trklenx_cross[ntracks_cross]/F");
      fEventTree->Branch("trklenx_selec", trklenx_selec, "trklenx_selec[ntracks_selec]/F");

      fEventTree->Branch("ntrkhits_selec",ntrkhits_selec,"ntrkhits_selec[ntracks_selec][3]/S"); 
   }//fSaveTrackInfo     
}

//========================================================================
void Lifetime::endJob(){     

}

//========================================================================
void Lifetime::beginRun(const art::Run& /*run*/){
  mf::LogInfo("Lifetime")<<"begin run..."<<std::endl;
}
//========================================================================
void Lifetime::analyze( const art::Event& evt ){
    if (!evt.isRealData()) return;
    reset();

    event  = evt.id().event(); 
    run    = evt.run();
    subrun = evt.subRun();
    
    art::Timestamp ts = evt.time();   
    TTimeStamp tts(ts.timeHigh(), ts.timeLow());
    evttime = tts.AsDouble();        

    art::Handle< std::vector<recob::Track> > trackListHandle;
    std::vector<art::Ptr<recob::Track> > tracklist;
    if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
       art::fill_ptr_vector(tracklist, trackListHandle);
       
    size_t NTracks = tracklist.size(); 
    ntracks = (int) NTracks;
    
    size_t cross=0, selec=0;
    for(unsigned int i=0; i<NTracks;++i){//loop over tracks
      art::Ptr<recob::Track> ptrack(trackListHandle, i);
      const recob::Track& track = *ptrack;    
      
      TVector3 pos, dir_start, dir_end, end;  
      
      double tlen = 0.,mom = 0.;
      int ntraj = 0;	
      ntraj = track.NumberTrajectoryPoints();

      if (ntraj > 0) {
        pos	 = track.Vertex();
        dir_start = track.VertexDirection();
        dir_end   = track.EndDirection();
        end	 = track.End();
        tlen	 = length(track);
	if(track.NumberFitMomentum() > 0)
     	     mom = track.VertexMomentum();
	trklen[i] = tlen;
	double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
        double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());	
	//get the X projected length
     	float X = std::abs(pos.X()-end.X());
	trklenx[i] = (X);
	trkthetaxz[i] = theta_xz;
	trkthetayz[i] = theta_yz;
     	//check if it is a crossing track
        if (X>250 && X<270){
	  cross++;
	  //save track information of crossing tracks
	  trklenx_cross[i] = X;
	  trkthetaxz_cross[i] = theta_xz;
	  trkthetayz_cross[i] = theta_yz;
	  trkstartx_cross[i]        = pos.X();
	  trkstarty_cross[i]        = pos.Y();
	  trkstartz_cross[i]        = pos.Z();
	  trkendx_cross[i]	    = end.X();
	  trkendy_cross[i]	    = end.Y();
	  trkendz_cross[i]	    = end.Z();
	  trktheta_cross[i]	    = dir_start.Theta();
	  trkphi_cross[i]	    = dir_start.Phi();
	  trkstartdcosx_cross[i]    = dir_start.X();
	  trkstartdcosy_cross[i]    = dir_start.Y();
	  trkstartdcosz_cross[i]    = dir_start.Z();
	  trkenddcosx_cross[i]      = dir_end.X();
	  trkenddcosy_cross[i]      = dir_end.Y();
	  trkenddcosz_cross[i]      = dir_end.Z();
	  trkthetaxz_cross[i]       = theta_xz;
	  trkthetayz_cross[i]       = theta_yz;
	  trkmom_cross[i]	    = mom;
	  trklen_cross[i]	    = tlen;  
	  //only accept tracks that are at a certain angle to remove track reconstruction issues
	  if (!((TMath::Abs(theta_xz) > 1.484) && (TMath::Abs(theta_xz) < 1.658))){ 
	    selec++;
	    trklenx_selec[i] = X;
	    trkthetaxz_selec[i] = theta_xz;
	    trkthetayz_selec[i] = theta_yz;	      
  	    art::FindMany<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
	    if (fmcal.isValid()){
     	      std::vector<const anab::Calorimetry*> calos = fmcal.at(i);
     	      if (calos.size() > 3) {
     	        // if you get this message, there is probably a bug somewhere since
     	        // the calorimetry planes should be 3.
     	        mf::LogError("Lifetime:limits")
     	        << "the " << fTrackModuleLabel << " track #" << i
     	        << " has " << calos.size() << " planes for calorimetry , only 3"
     	        << " stored in tree";
     	      }//if (calos.size() > 3)
    	      for (size_t ical = 0; ical<calos.size(); ++ical){
     	        if (!calos[ical]) continue;
     	        if (!calos[ical]->PlaneID().isValid) continue;
     	        int planenum = calos[ical]->PlaneID().Plane;
     	        if (planenum<0||planenum>2) continue;
	        const size_t NHits = calos[ical] -> dEdx().size();
     	        ntrkhits_selec[i][planenum] = (int) NHits;	    
	        if (planenum == 2){ //only interested in plane 2 for now
     	   	  double minx = 1e10;
     	   	  for(size_t iHit = 0; iHit < NHits; ++iHit) {	      
     	   	    const auto& TrkPos = (calos[ical] -> XYZ())[iHit];
     	   	    if (TrkPos.X()<minx)
     	   	    minx = TrkPos.X();
     	   	  }// loop NHits
     	   	  for(size_t iHit = 0; iHit < NHits; ++iHit) {				
     	   	    const auto& TrkPos1 = (calos[ical] -> XYZ())[iHit];
     	   	    double x = TrkPos1.X()-minx; //subtract the minx to get correct t0
	   	    double t = x/(XDriftVelocity*1000); //change the velocity units to cm/ns to cm/us
     	   	    dqdstime->Fill(t, (calos[ical] -> dQdx())[iHit]);
     	   	    int bin = int(t/binsize);
     	   	    if (bin>=0&&bin<nbin)
     	   	      dqdx[bin]->Fill((calos[ical] -> dQdx())[iHit]);
     	   	  }// loop NHits 
     	       } // if planenum ==2	  
     	     }// loop over ical	
     	   }// if fmcal.isValid()
     	 }// if (!((TMath::Abs(theta_xz) cut 
       }// if (X>250 && X<270)
     }// if ntraj>0
   }// loop over tracks 
   ntracks_cross = (int) cross;
   ntracks_selec = (int) selec;
   
   fEventTree->Fill();

      
}// end of analyze function




//========================================================================
void Lifetime::reset(){
  run = -99999;
  subrun = -99999;
  event = -99999;
  evttime = -99999;
  
if (fSaveTrackInfo){
  ntracks = 0;
  ntracks_cross =0;
  ntracks_selec =0;
  for (int i = 0; i < kMaxTrack; ++i){
    trkstartx_cross[i] = -99999.;   
    trkstarty_cross[i] = -99999.;   
    trkstartz_cross[i] = -99999.;   
    trkendx_cross[i] = -99999.;     
    trkendy_cross[i] = -99999.;     
    trkendz_cross[i] = -99999.;     
    trktheta_cross[i] = -99999.;    
    trkphi_cross[i] = -99999.; 
    trkstartdcosx_cross[i] = -99999.;     
    trkstartdcosy_cross[i] = -99999.;     
    trkstartdcosz_cross[i] = -99999.;  
    trkenddcosx_cross[i] = -99999.;
    trkenddcosy_cross[i] = -99999.;
    trkenddcosz_cross[i] = -99999.;
    trkthetaxz[i] = -99999.;   
    trkthetaxz_cross[i] = -99999.;
    trkthetaxz_selec[i] = -99999.;
    trkthetayz[i] = -99999.;        
    trkthetayz_cross[i] = -99999.; 
    trkthetayz_selec[i] = -99999.; 
    trkmom_cross[i] = -99999.;      
    trklen_cross[i] = -99999.;
    trklenx[i] = -99999.;
    trklenx_cross[i] = -99999.;
    trklenx_selec[i] = -99999.;
    
    for (int j = 0; j < kNplanes; ++j){
      ntrkhits_selec[i][j] = -9999;
    }
  }
}//fSaveTrackInfo

  
}

//========================================================================
DEFINE_ART_MODULE(Lifetime)

} 

#endif // Lifetime_Module
