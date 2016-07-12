#ifndef LifetimeQAQC_Module
#define LifetimeQAQC_Module

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

//========================================================================
double langaufun(double *x, double *par){
   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation), 
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      //double mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      double np = 100.0;      // number of convolution steps
      double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      double xx;
      double mpc;
      double fland;
      double sum = 0.0;
      double xlow,xupp;
      double step;
      double i;

      // MP shift correction
      mpc = par[1];

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[2] * step * sum * invsq2pi / par[3]);
}

//========================================================================
TF1 *langaufit(TH1D *his, double *fitrange, double *startvalues, double *parlimitslo, double *parlimitshi, double *fitparams, double *fiterrors, double *ChiSqr, int *NDF){
   // Once again, here are the Landau * Gaussian parameters:
   //   par[0]=Width (scale) parameter of Landau density
   //   par[1]=Most Probable (MP, location) parameter of Landau density
   //   par[2]=Total area (integral -inf to inf, normalization constant)
   //   par[3]=Width (sigma) of convoluted Gaussian function
   //
   // Variables for langaufit call:
   //   his             histogram to fit
   //   fitrange[2]     lo and hi boundaries of fit range
   //   startvalues[4]  reasonable start values for the fit
   //   parlimitslo[4]  lower parameter limits
   //   parlimitshi[4]  upper parameter limits
   //   fitparams[4]    returns the final fit parameters
   //   fiterrors[4]    returns the final fit errors
   //   ChiSqr          returns the chi square
   //   NDF             returns ndf

   int i;
   char FunName[100];

   sprintf(FunName,"Fitfcn_%s",his->GetName());

   art::ServiceHandle<art::TFileService> tfs;

   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
   if (ffitold) delete ffitold;

   TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Width","MPV","Area","GSigma");
   
   for (i=0; i<4; i++) {
      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
   }

   //his->Fit(FunName,"RB0ME");   // fit within specified range, use ParLimits, do not plot
   his->Fit(FunName,"RB0");   // fit within specified range, use ParLimits, do not plot

   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<4; i++) {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();           // obtain ndf

   return (ffit);              // return fit function

}

//========================================================================
int langaupro(double *params, double &maxx, double &FWHM){
   // Seaches for the location (x value) at the maximum of the 
   // Landau-Gaussian convolute and its full width at half-maximum.
   //
   // The search is probably not very efficient, but it's a first try.

   double p,x,fy,fxr,fxl;
   double step;
   double l,lold;
   int i = 0;
   int MAXCALLS = 10000;

   // Search for maximum

   p = params[1] - 0.1 * params[0];
   step = 0.05 * params[0];
   lold = -2.0;
   l    = -1.0;

   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;
      lold = l;
      x = p + step;
      l = langaufun(&x,params);
 
      if (l < lold)
         step = -step/10;
 
      p += step;
   }

   if (i == MAXCALLS)
      return (-1);

   maxx = x;

   fy = l/2;

   // Search for right x location of fy
   p = maxx + params[0];
   step = params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;

   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;
      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);
 
      if (l > lold)
         step = -step/10;
 
      p += step;
   }

   if (i == MAXCALLS)
      return (-2);

   fxr = x;

   // Search for left x location of fy
   p = maxx - 0.5 * params[0];
   step = -params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;

   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;
      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);
 
      if (l > lold)
         step = -step/10;
 
      p += step;
   }

   if (i == MAXCALLS)
      return (-3);

   fxl = x;

   FWHM = fxr - fxl;
   return (0);
}

//========================================================================
void langaus(TH1D** dqdx) {
  art::ServiceHandle<art::TFileService> tfs;
  TH1D *hsigma = tfs->make<TH1D>("hsigma","hsigma",50,0,100);
  hsigma->Sumw2();
     
  vector<double> time;
  vector<double> etime;
  vector<double> charge;
  vector<double> echarge;

  double totalndf = 0;
  double totalchi2 = 0;

  for (int i = 0; i<nbin; ++i){
    printf("Fitting...\n");
    //TH1D *hSNR= (TH1D*)Form("dqdx_%d",i);    
    TH1D *hSNR= (TH1D*)dqdx[i];    

    double fr[2];
    double sv[4], pllo[4], plhi[4], fp[4], fpe[4];
    std::cout<<"\n"<<0.2*hSNR->GetMean()<<"\n";
    fr[0]=0.2*hSNR->GetMean();
    fr[1]=4.0*hSNR->GetMean();
    pllo[0]=10; pllo[1]=100.0; pllo[2]=1.0; pllo[3]=10;
    plhi[0]=1000.0; plhi[1]=3000.0; plhi[2]=10000000.0; plhi[3]=1000.0;
    sv[0]=0.13*hSNR->GetRMS(); sv[1]=0.8*hSNR->GetMean(); sv[2]=hSNR->GetEntries()*50; sv[3]=100.0;

    double chisqr;
    int    ndf;
    TF1 *fitsnr = langaufit(hSNR,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
    double SNRPeak, SNRFWHM;
    langaupro(fp,SNRPeak,SNRFWHM);
    std::cout<<fitsnr->GetParameter(0)<<" "<<
      fitsnr->GetParameter(1)<<" "<<
      fitsnr->GetParameter(2)<<" "<<
      fitsnr->GetParameter(3)<<endl;
      
    if (gMinuit&&hSNR->GetEntries()>100){
      TString test =  gMinuit->fCstatu.Data(); 
      if (test.EqualTo("CONVERGED ")&&fitsnr->GetChisquare()/fitsnr->GetNDF()<10&&fitsnr->GetParError(1)<1000){
	std::cout<<fitsnr->GetChisquare()<<" "<<fitsnr->GetNDF()<<" "<<fitsnr->GetChisquare()/fitsnr->GetNDF()<<endl;
	totalndf += fitsnr->GetNDF();
	totalchi2 += fitsnr->GetChisquare();
	time.push_back((i*binsize+binsize/2));
	etime.push_back(binsize/2);
	hsigma->Fill(fitsnr->GetParameter(3));
	charge.push_back(fitsnr->GetParameter(1));
	echarge.push_back(fitsnr->GetParError(1));
      }//end  if (test.EqualTo("CONVERGED ") ....
    }// end if (gMinuit&&hSNR->GetEntries()>100)
  }// if nbin 

  cout<<"Total ndf = "<<totalndf<<endl;
  cout<<"Total chi2 = "<<totalchi2<<endl;
  cout<<"\n"<<time.size();
  for(unsigned int t=0;t<time.size();t++){
  cout<<"\n"<<time[t]<<","<<etime[t]<<","<<charge[t]<<","<<echarge[t];
  }
  TGraphErrors *gr = tfs->make<TGraphErrors>(time.size(),&time[0],&charge[0],&etime[0],&echarge[0]);
  gr->GetXaxis()->SetTitle("err");
  gr->SetMarkerStyle(21);
  gr->SetMarkerColor(4);
  gr->GetXaxis()->SetTitle("Drift time (#mus)");
  gr->GetYaxis()->SetTitle("dQ/ds (ADC/cm)");
  gr->GetXaxis()->SetLabelSize(0.04);
  gr->GetYaxis()->SetLabelSize(0.04);
  gr->GetYaxis()->SetTitleOffset(1.25);
  gr->GetXaxis()->SetTitleOffset(1.25);
  
  //Fit an exponential + constant function to the TGraph
  TF1 *fit = new TF1("fit","[0] + exp([1] + [2]*x)",50,2150);
  fit->SetParameter(0,0);
  fit->SetParameter(1,5.46113);
  fit->SetParameter(2,-0.000127173);
  
  gr->Fit("fit","R"); 
  gr->Draw("AP");
  fit->Draw("same");
  gr->Write();
  
  Double_t par[3];
  Double_t par_error[3];
  Double_t QA,QC;
  Double_t lifetime;
  Double_t QA_err, QC_err, lifetime_err,QAQC_err;
  fit->GetParameters(&par[0]); 
  for (int i =0; i <3; i++ )
       par_error[i] = fit->GetParError(i);
  
  QC = par[0] + TMath::Exp(par[1]);
  QA = par[0] + TMath::Exp(par[1] + par[2]*2200);
  lifetime = -1/par[2];
  
  QC_err = TMath::Sqrt(TMath::Power(par_error[0],2) + (TMath::Exp(2*par[1]))*(TMath::Power(par_error[1],2)));
  QA_err = TMath::Sqrt(TMath::Power(par_error[0],2) + TMath::Exp(par[1] + par[2]*2200)*(TMath::Power(par_error[1],2) + (TMath::Power(par_error[2],2))*2200*2200));
  lifetime_err = TMath::Abs((1/(par[1]*par[1]))*par_error[2]);
  QAQC_err = (QA + QA_err)/(QC - QC_err) - QA/QC;
  
  std::cout << "***************Purity Information *******************" << endl;
  std::cout << "QC    : " << QC << " +/- " << QC_err << endl;
  std::cout << "QA    : " << QA << " +/- " << QA_err << endl;
  std::cout << "QA/QC : " << QA/QC << " +/- " << QAQC_err << endl;
  std::cout << "Electron Life : " << lifetime << " +/- " << lifetime_err << endl;
  
  const int n=1;
  Double_t x[n] = {0.5}, ex[n]={0};
  Double_t QaQc[n], QaQc_err[n];
  QaQc[0]=QA/QC;
  QaQc_err[0]=QAQC_err;
  
  TCanvas *qAqC = tfs->make<TCanvas>("qAqC","");
   qAqC->SetGridy(1);
   qAqC->SetTicky(1);

  TGraphErrors *qaqc = tfs->make<TGraphErrors>(1,x,QaQc,ex,QaQc_err);
  qaqc->SetMarkerStyle(20);
  qaqc->SetMarkerColor(kRed);
  qaqc->GetXaxis()->SetNdivisions(0,2,0);
  qaqc->GetYaxis()->SetRangeUser(0,1.2);
  qaqc->GetXaxis()->SetRangeUser(0.4,0.604);
  qaqc->SetTitle(" ");  
  qaqc->GetYaxis()->SetTitle("QA/QC");
  qaqc->Draw("AP");
  
  /*TLine *line = new TLine(0.4,1,0.6,1);
  line->SetLineColor(kBlack);
  line->SetLineWidth(2);
  line->Draw("same");*/

  qaqc->Write();
  qAqC->Write();
}

} // end namespace

//========================================================================

namespace microboone{

class LifetimeQAQC : public art::EDAnalyzer {
public:

    explicit LifetimeQAQC(fhicl::ParameterSet const& pset);
    virtual ~LifetimeQAQC();

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

 }; // class LifetimeQAQC


//========================================================================
LifetimeQAQC::LifetimeQAQC(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
{
    reconfigure(parameterSet);
}
//========================================================================
LifetimeQAQC::~LifetimeQAQC(){
  //destructor
}
//========================================================================
void LifetimeQAQC::reconfigure(fhicl::ParameterSet const& p){

    fTrackModuleLabel    = p.get<std::string>("TrackModuleLabel");
    fCalorimetryModuleLabel    = p.get<std::string>("CalorimetryModuleLabel");

}
//========================================================================
void LifetimeQAQC::beginJob(){
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
void LifetimeQAQC::endJob(){     

 langaus(dqdx);
}

//========================================================================
void LifetimeQAQC::beginRun(const art::Run& /*run*/){
  mf::LogInfo("LifetimeQAQC")<<"begin run..."<<std::endl;
}
//========================================================================
void LifetimeQAQC::analyze( const art::Event& event ){
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
	//double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
	//only accept tracks that are at a certain angle to remove track reconstruction issues
	//if (!((TMath::Abs(theta_xz) > 1.484) && (TMath::Abs(theta_xz) < 1.658))){
     	  //get the X projected length
     	  float X = std::abs(pos.X()-end.X());
     	  //check if it is a crossing track
     	  if (X>250 && X<270)//{
     	    len_cross->Fill(X); 	      
     	    art::FindMany<anab::Calorimetry> fmcal(trackListHandle, event, fCalorimetryModuleLabel);
     	    if (fmcal.isValid()){
     	       std::vector<const anab::Calorimetry*> calos = fmcal.at(i);
     	       if (calos.size() > 3) {
     	  	  // if you get this message, there is probably a bug somewhere since
     	  	  // the calorimetry planes should be 3.
     	  	  mf::LogError("LifetimeQAQC:limits")
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
     	   //}// if (X>250 && X<270)
	// }// if (!((TMath::Abs(theta_xz) cut 
       }// if ntraj>0
    }// loop over tracks    
}// end of analyze function


//========================================================================
void LifetimeQAQC::reset(){
//do nothing
}

//========================================================================
DEFINE_ART_MODULE(LifetimeQAQC)

} 

#endif // LifetimeQAQC_Module
