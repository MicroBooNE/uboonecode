#ifndef Lifetime_Module
#define Lifetime_Module

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardata/RecoBase/Track.h"
#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"
#include "larcore/SimpleTypesAndConstants/geo_types.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/RawData/ExternalTrigger.h"
#include "art/Framework/Core/FindMany.h"
#include "lardata/AnalysisBase/Calorimetry.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Persistency/Common/PtrVector.h"
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
     
    TTree *fEventTree;
    TH1D *len;
    TH1D *len_cross;    
    TH2D *dqdstime;    
    TH1D *dqdx[nbin];
 
    // Event 
    int Event;
    int Run;
    int SubRun;

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

    fTrackModuleLabel    = p.get<std::string>("TrackModuleLabel");
    fCalorimetryModuleLabel    = p.get<std::string>("CalorimetryModuleLabel");

}
//========================================================================
void Lifetime::beginJob(){
  std::cout<<"job begin..."<<std::endl;

  art::ServiceHandle<art::TFileService> tfs;

  len = tfs->make<TH1D>("len","Length of tracks (cm)",100,0,1000);
  len_cross = tfs->make<TH1D>("len_cross","Length of crossing tracks (cm)",100,0,1000);
  
  dqdstime = tfs->make<TH2D>("dqdstime","",100,0,2200,100,0,1000);
  dqdstime->GetXaxis()->SetTitle("Drift time (#mus)");
  dqdstime->GetYaxis()->SetTitle("dQ/ds (ADC/cm)");
  
  for (int i = 0; i<nbin; ++i){
      dqdx[i] = tfs->make<TH1D>(Form("dqdx_%d",i),Form("dqdx_%d",i),100,0,1500);
      dqdx[i]->GetXaxis()->SetTitle("dq/dx (ADC/cm)");
      dqdx[i]->Sumw2();
    }

  len->Sumw2();
  len_cross->Sumw2();

  //if( fSaveMCTree ){
    fEventTree = tfs->make<TTree>("Event", "Event Tree from Sim & Reco");
    fEventTree->Branch("eventNo", &Event);
    fEventTree->Branch("runNo", &Run);
    fEventTree->Branch("subRunNo", &SubRun);    
  //}
}

//========================================================================
void Lifetime::endJob(){     

}

//========================================================================
void Lifetime::beginRun(const art::Run& /*run*/){
  mf::LogInfo("Lifetime")<<"begin run..."<<std::endl;
}
//========================================================================
void Lifetime::analyze( const art::Event& event ){
    if (!event.isRealData()) return;
    reset();

    Event  = event.id().event(); 
    Run    = event.run();
    SubRun = event.subRun();
    fEventTree->Fill();
    
    
    art::Handle< std::vector<recob::Track> > trackListHandle;
    std::vector<art::Ptr<recob::Track> > tracklist;
    if (event.getByLabel(fTrackModuleLabel,trackListHandle))
       art::fill_ptr_vector(tracklist, trackListHandle);
       
    size_t NTracks = tracklist.size(); 
    
    for(unsigned int i=0; i<NTracks;++i){//loop over tracks
      art::Ptr<recob::Track> ptrack(trackListHandle, i);
      const recob::Track& track = *ptrack;    
      
      TVector3 pos, dir_start, dir_end, end;  
      
      double tlen = 0.;
      int ntraj = 0;	  
      
      ntraj = track.NumberTrajectoryPoints();
      if (ntraj > 0) {
        pos	 = track.Vertex();
        dir_start = track.VertexDirection();
        dir_end   = track.EndDirection();
        end	 = track.End();
        tlen	 = length(track);
	len->Fill(tlen);
	double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
	//only accept tracks that are at a certain angle to remove track reconstruction issues
	if (!((TMath::Abs(theta_xz) > 1.484) && (TMath::Abs(theta_xz) < 1.658))){
     	  //get the X projected length
     	  float X = std::abs(pos.X()-end.X());
     	  //check if it is a crossing track
     	  if (X>250 && X<270){
     	    len_cross->Fill(X); 	      
     	    art::FindMany<anab::Calorimetry> fmcal(trackListHandle, event, fCalorimetryModuleLabel);
     	    if (fmcal.isValid()){
     	       std::vector<const anab::Calorimetry*> calos = fmcal.at(i);
     	       if (calos.size() > 3) {
     	  	  // if you get this message, there is probably a bug somewhere since
     	  	  // the calorimetry planes should be 3.
     	  	  mf::LogError("Lifetime:limits")
     	  	  << "the " << fTrackModuleLabel << " track #" << i
     	  	  << " has " << calos.size() << " planes for calorimetry , only 3"
     	  	  << " stored in tree";
     	       }
     	       for (size_t ical = 0; ical<calos.size(); ++ical){
     	  	 if (!calos[ical]) continue;
     	  	 if (!calos[ical]->PlaneID().isValid) continue;
     	  	 int planenum = calos[ical]->PlaneID().Plane;
     	  	 if (planenum<0||planenum>2) continue;
     	  	 if (planenum == 2){ //only interested in plane 2 for now
     	  	     const size_t NHits = calos[ical] -> dEdx().size();
     	  	     //std::cout<<"\n"<<NHits;
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
     	   }// if (X>250 && X<270)
	 }// if (!((TMath::Abs(theta_xz) cut 
       }// if ntraj>0
    }// loop over tracks    
}// end of analyze function


//========================================================================
void Lifetime::reset(){
//do nothing
}

//========================================================================
DEFINE_ART_MODULE(Lifetime)

} 

#endif // Lifetime_Module
