//-----------------------------------------------------------------------
//
//	Convoluted Landau and Gaussian Fitting Function for Lifetime analysis
//         (using ROOT's Landau and Gauss functions)
//
//  Based on a Fortran code by R.Fruehwirth (fruhwirth@hephy.oeaw.ac.at)
//  Adapted for C++/ROOT by H.Pernegger (Heinz.Pernegger@cern.ch) and
//   Markus Friedl (Markus.Friedl@cern.ch)
//
//  to execute this example, do:
//  root > .x Lifetime.C
// or
//  root > .x Lifetime.C++
//
//-----------------------------------------------------------------------

#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TVectorD.h"
#include <vector>

using namespace std;

Double_t langaufun(Double_t *x, Double_t *par) {

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
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;


      // MP shift correction
      //mpc = par[1] - mpshift * par[0]; 
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



TF1 *langaufit(TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
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

   Int_t i;
   Char_t FunName[100];

   sprintf(FunName,"Fitfcn_%s",his->GetName());

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


Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {

   // Seaches for the location (x value) at the maximum of the 
   // Landau-Gaussian convolute and its full width at half-maximum.
   //
   // The search is probably not very efficient, but it's a first try.

   Double_t p,x,fy,fxr,fxl;
   Double_t step;
   Double_t l,lold;
   Int_t i = 0;
   Int_t MAXCALLS = 10000;

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

void Lifetime(int itime = 0) {
   //input merged root file with dQ/ds scatter plot
   TFile file("/uboone/data/users/sowjanya/run6123_merged.root");
   //TFile file("/uboone/app/users/sowjanya/v580_lt/lifetime.root");
   
   //create on output file to write out important plots
   TFile f("lifetime_qaqc_out.root","recreate");
   
   gStyle->SetOptStat(1111);
   gStyle->SetOptFit(111);
   gStyle->SetOptStat(0);
   gStyle->SetLabelSize(0.03,"x");
   gStyle->SetLabelSize(0.03,"y");
   
   TH1D *hsigma = new TH1D("hsigma","hsigma",50,0,100);
   hsigma->Sumw2();
   
   TH2D *dqdstime = (TH2D*)file.Get("Lifetime/dqdstime");
   dqdstime->Write();

   const int tbinsize = 200;
   const int ntbins = 11;
  
   vector<double> time;
   vector<double> etime;
   vector<double> charge;
   vector<double> echarge;
  
   double totalndf = 0;
   double totalchi2 = 0;
   
   TCanvas *c = new TCanvas("c","c");

   for (int i = 0; i<ntbins; ++i){
    int bin=i*tbinsize;
    TCanvas *c = new TCanvas(Form("dqdx_fit_%dus",bin),Form("dqdx_fit_%dus",bin));
    c->cd();
    
    TH1D *hSNR = (TH1D*)file.Get(Form("Lifetime/dqdx_%d",i));

    printf("Fitting...\n");    
    Double_t fr[2];
    Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
    fr[0]=0.2*hSNR->GetMean();
    fr[1]=4.0*hSNR->GetMean();    
    pllo[0]=10; pllo[1]=100.0; pllo[2]=1.0; pllo[3]=10;
    plhi[0]=1000.0; plhi[1]=3000.0; plhi[2]=10000000.0; plhi[3]=1000.0;
    sv[0]=0.13*hSNR->GetRMS(); sv[1]=0.8*hSNR->GetMean(); sv[2]=hSNR->GetEntries()*50; sv[3]=100.0;

    Double_t chisqr;
    Int_t    ndf;
    TF1 *fitsnr = langaufit(hSNR,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
    Double_t SNRPeak, SNRFWHM;
    langaupro(fp,SNRPeak,SNRFWHM);
    cout<<fitsnr->GetParameter(0)<<" "<<
      fitsnr->GetParameter(1)<<" "<<
      fitsnr->GetParameter(2)<<" "<<
      fitsnr->GetParameter(3)<<endl;
      
    if (gMinuit&&hSNR->GetEntries()>100){
      TString test =  gMinuit->fCstatu.Data(); 
      if (test.EqualTo("CONVERGED ")&&fitsnr->GetChisquare()/fitsnr->GetNDF()<10&&fitsnr->GetParError(1)<1000){
	cout<<fitsnr->GetChisquare()<<" "<<fitsnr->GetNDF()<<" "<<fitsnr->GetChisquare()/fitsnr->GetNDF()<<endl;
	totalndf += fitsnr->GetNDF();
	totalchi2 += fitsnr->GetChisquare();
	time.push_back((i*tbinsize+tbinsize/2));
	etime.push_back(tbinsize/2);
	hsigma->Fill(fitsnr->GetParameter(3));
	charge.push_back(fitsnr->GetParameter(1));
	echarge.push_back(fitsnr->GetParError(1));
      } //end  if (test.EqualTo("CONVERGED ") ....
    } // end if (gMinuit&&hSNR->GetEntries()>100)

   //the command below doesn't work for some reason, investigate later...
   //hSNR->SetTitle(Form("%.0f<Drift Time<%.0f",i*tbinsize,(i+1)*tbinsize));
   hSNR->SetXTitle("dQ/ds (ADC/cm)");
   hSNR->SetYTitle("Number of hits");
   hSNR->GetXaxis()->SetRangeUser(0,5000);
   hSNR->GetXaxis()->SetLabelSize(0.04);
   hSNR->GetYaxis()->SetLabelSize(0.04);
   hSNR->GetYaxis()->SetTitleOffset(1.30);
   
   hSNR->DrawCopy();
   fitsnr->SetLineColor(2);
   fitsnr->SetLineWidth(2);
   fitsnr->Draw("lsame");
   hSNR->DrawCopy();
   fitsnr->Draw("lsame");
   c->Write();
  } // end ntbin

  printf("Fitting done\nPlotting results...\n");

  cout<<"Total ndf = "<<totalndf<<endl;
  cout<<"Total chi2 = "<<totalchi2<<endl;
  
  TCanvas *dqdsvsdt = new TCanvas("dqdsvsdt","dqdsvsdt");
  dqdsvsdt->cd();
  TGraphErrors *gr = new TGraphErrors(time.size(),&time[0],&charge[0],&etime[0],&echarge[0]);
  gr->SetTitle("");
  gr->GetXaxis()->SetTitle("Drift time (#mus)");
  gr->GetYaxis()->SetTitle("dQ/ds (ADC/cm)");
  gr->GetXaxis()->SetLabelSize(0.04);
  gr->GetYaxis()->SetLabelSize(0.04);
  gr->GetYaxis()->SetTitleOffset(1.25);
  gr->GetXaxis()->SetTitleOffset(1.25);

  gr->Draw("AP");   
  //gr->Write();
  dqdsvsdt->Write();
  hsigma->Write();
  
  //Fit an exponential + constant function to the TGraph
  TCanvas *dqdsvsdt_fit = new TCanvas("dqdsvsdt_fit","dqdsvsdt_fit");
  dqdsvsdt_fit->cd();
  TF1 *fit = new TF1("fit","[0] + exp([1] + [2]*x)",50,2150);
  fit->SetParameter(0,0);
  fit->SetParameter(1,5.46113);
  fit->SetParameter(2,-0.000127173);
  
  gr->Fit("fit","R"); 
  gr->Draw("AP");
  fit->Draw("same");
  dqdsvsdt_fit->Write();    
  gr->Write();
  
  //Now do expoenential (or expo+constant) fits to the gr TGraphs plot
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
  
  TCanvas *qAqC = TCanvas("qAqC","");
   qAqC->SetGridy(1);
   qAqC->SetTicky(1);

  TGraphErrors *qaqc = TGraphErrors(1,x,QaQc,ex,QaQc_err);
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

