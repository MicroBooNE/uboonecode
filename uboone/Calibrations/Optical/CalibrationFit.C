Double_t MPEfit(Double_t * x, Double_t * par) {

// Defines the multi-photoelectron (MPE) fit used to fit charge and amplitude histograms.
// Fit is comprised of a series of poisson convolved gaussians.

  double ped_mu     = par[0];     // Centroid of the gaussian fit to the pedestal
  double ped_sig    = par[1];     // FWHM of the gaussian fit to the pedestal
  double mu_pe      = par[2];     // Average number of photoelectrons seen, determines the behavior of the poisson terms
  double spe        = par[3];     // Centroid of the 1-PE gaussian peak, used to extract the SPE value
  double pe_sig     = par[4];     // Intrincic FWHM of the photoelectron peaks, actual peaks will have FWHM including the baseline spread (ped_sig)
  double gain       = par[5];     // Peak to peak distance of photoelectron peaks, i.e. distance between the n PE and n+1 PE centroids
  double mpe_norm   = par[6];     // Scaling constant for the sum of poisson convolved gaussians

  double MPE = 0;

  int NPEtoFit = 5;

  for (int npe = 1; npe <= NPEtoFit ; npe++) {

    double subgauss_sigma = TMath::Sqrt(npe*pe_sig*pe_sig + ped_sig*ped_sig);
    MPE = MPE + TMath::Exp(-mu_pe)*TMath::Power(mu_pe,npe)/TMath::Factorial(npe) *TMath::Gaus(x[0],spe + ped_mu + (npe-1)*gain,subgauss_sigma,kTRUE);

  }

  MPE = MPE*mpe_norm;

  return MPE;

    }

//--------------------------------------------------------------------------------------------------------------------------------------------
// Actual fitting code follows                                                                                                             ---
//--------------------------------------------------------------------------------------------------------------------------------------------

double CalibrationFit(TString input_filename, int channel_num) {

  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TString input_tree_name = "outtree";                                   // Whatever the analyzed tree is named
  string charge           = "charge";                                    // Whatever integrated peaks are called in the tree
  string amplitude        = "maxamp";                                    // Whatever peak maxima are called in the tree
  string baselinerms      = "baselinerms";                               // Whatever the baseline rms for the pre charge region is called in the tree
  string baselinerms2     = "baselinerms2";                              // Whatever the baseline rms for the post charge region is called in the tree
  string channel          = "opchannel";                                 // Whatever the channel number is called in the tree

  char base_ratio[100];
  sprintf(base_ratio,"%s / %s",baselinerms.c_str(),baselinerms2.c_str());


  int charge_min          = 70;                                          // Determines where the pedestal region ends for charge
  int charge_max          = 500;                                         // Determines where the region to be fit ends for charge
  int amp_min             = 10;                                          // Determines where the pedestal region ends for amplitude
  int amp_max             = 70;                                          // Determines where the region to be fit ends for amplitude
  int baseline_cutoff     = 3;                                           // Determines the maximum variability in baseline allowed
  double ratio_threshold  = 0.05;                                        // Determines the amount by which the post-charge baseline can vary from the pre-charge baseline
  

  TFile* in_file = new TFile(input_filename.Data());                     // Get input file output from swizzler
  TTree* in_tree = (TTree*)in_file->Get(input_tree_name);                // Get tree with PMT waveforms and data from input file


//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// This following segment histograms the charge and amplitudes and obtains an estimate for the position of the SPE peak to be                                  
// used to seed the actual fit later.                                                                                                                          

  TH1D* charge_histo = new TH1D("charge_histo","Charge Histogram",100,0,1000);
  TH1D* amp_histo    = new TH1D("amp_histo","Amplitude Histogram",100,0,100);

  char buffer1[300];
  char buffer2[300];
  sprintf(buffer1,"%s>%i && %s<%i && %s<%i && %s<(1+%f) && %s>(1-%f) && %s==%i",amplitude.c_str(),amp_min,baselinerms.c_str(),baseline_cutoff,baselinerms2.c_str(),baseline_cutoff,base_ratio,ratio_threshold,base_ratio,ratio_threshold,channel.c_str(),channel_num);
  sprintf(buffer2,"%s>%i && %s<%i && %s<%i && %s<(1+%f) && %s>(1-%f) && %s==%i",charge.c_str(),charge_min,baselinerms.c_str(),baseline_cutoff,baselinerms2.c_str(),baseline_cutoff,base_ratio,ratio_threshold,base_ratio,ratio_threshold,channel.c_str(),channel_num);
  in_tree->Project("amp_histo",amplitude.c_str(),buffer1);
  in_tree->Project("charge_histo",charge.c_str(),buffer2);

  double spe_charge_estimate = charge_histo->GetXaxis()->GetBinCenter(charge_histo->GetMaximumBin());
  double spe_amp_estimate = amp_histo->GetXaxis()->GetBinCenter(amp_histo->GetMaximumBin());


//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// The following couple lines gives an estimate of the mpe normalization by simply finding the size of the 1 PE hump, note the                                 
// factor of 2.5066 is to account for the sqrt(2pi) norm factor in the gaussians, doesn't account for the sigma as we dont know                                
// that until we fit anyway, simply a rough initial starting point.                                                                                            

  double charge_mpe_norm_estimate = 2.5066*charge_histo->GetMaximum();                                       
  double amp_mpe_norm_estimate = 2.5066*amp_histo->GetMaximum();


//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// This following segment histograms the charge and amplitudes in the pedestal region and obtains an estimate or the position                                  
// of the charge and amplitude pedestal means                                                                                                                 

  TH1D* charge_ped_histo = new TH1D("charge_ped_histo","Charge Ped Histogram",100,-50,50);
  TH1D* amp_ped_histo    = new TH1D("amp_ped_histo","Amplitude Ped Histogram",100,0,5);

  char buffer3[300];
  char buffer4[300];

  sprintf(buffer3,"%s < %i && %s < %i && %s < %i && %s<(1+%f) && %s>(1-%f) && %s==%i",charge.c_str(),charge_min,baselinerms.c_str(),baseline_cutoff,baselinerms2.c_str(),baseline_cutoff,base_ratio,ratio_threshold,base_ratio,ratio_threshold,channel.c_str(),channel_num);
  sprintf(buffer4,"%s < %i && %s < %i && %s < %i && %s >0 && %s<(1+%f) && %s>(1+%f) && %s==%i",amplitude.c_str(),amp_min,baselinerms.c_str(),baseline_cutoff,baselinerms2.c_str(),baseline_cutoff,amplitude.c_str(),base_ratio,ratio_threshold,base_ratio,ratio_threshold,channel.c_str(),channel_num);
  
  in_tree->Project("charge_ped_histo",charge.c_str(),buffer3);
  in_tree->Project("amp_ped_histo",amplitude.c_str(),buffer4);

  double ped_charge_estimate = charge_ped_histo->GetXaxis()->GetBinCenter(charge_ped_histo->GetMaximumBin());
  double ped_amp_estimate = amp_ped_histo->GetXaxis()->GetBinCenter(amp_ped_histo->GetMaximumBin());

//---------------------------------------------------------------------------------------------------------------------------------------------------------------
// Makes the assumption that all events smaller than the midpoint of the 1PE and pedestal peaks are zero events                                               
// and uses this to get an estimate for the mean number of photoelectrons                                                                                      ///

  char buffer5[300];
  char buffer6[300];
  char buffer7[300];

  sprintf(buffer5,"%s < %f && %s < %i && %s < %i && %s<(1+%f) && %s>(1-%f) && %s==%i",charge.c_str(),(spe_charge_estimate + ped_charge_estimate)/2.0,baselinerms.c_str(),baseline_cutoff,baselinerms2.c_str(),baseline_cutoff,base_ratio,ratio_threshold,base_ratio,ratio_threshold,channel.c_str(),channel_num);
  
  sprintf(buffer6,"%s < %f && %s < %i && %s < %i && %s<(1+%f) && %s>(1-%f) && %s==%i",amplitude.c_str(),(spe_amp_estimate + ped_amp_estimate)/2.0,baselinerms.c_str(),baseline_cutoff,baselinerms2.c_str(),baseline_cutoff,base_ratio,ratio_threshold,base_ratio,ratio_threshold,channel.c_str(),channel_num);

 sprintf(buffer7,"%s < %i && %s < %i && %s<(1+%f) && %s>(1-%f) && %s==%i",baselinerms.c_str(),baseline_cutoff,baselinerms2.c_str(),baseline_cutoff,base_ratio,ratio_threshold,base_ratio,ratio_threshold,channel.c_str(),channel_num);
 

  double charge_zero_fraction = (double(in_tree->GetEntries(buffer5)))/(in_tree->GetEntries(buffer7));
  double amp_zero_fraction = (double(in_tree->GetEntries(buffer6)))/(in_tree->GetEntries(buffer7));

  double charge_mu_pe_estimate = TMath::Log(1.0/charge_zero_fraction);
  double amp_mu_pe_estimate = TMath::Log(1.0/amp_zero_fraction);
    

//---------------------------------------------------------------------------------------------------------------------------------------------------------
// Performs the actual fitting, seeding with the estimates made above                                                                                         

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);  // Raises maximum minimization calls, decreases likelikhood of early termination

  // This first half does the fitting to the amplitude histogram
  TH1F* amp_fit_histo = new TH1F("amp_fit_histo","Amplitude Histogram With Fit",100,0,50);

  amp_fit_histo->SetTitle("Amplitude Distribution With Fit");
  amp_fit_histo->GetYaxis()->SetTitle("Counts");
  amp_fit_histo->GetXaxis()->SetTitle("Amplitude (ADCs)");

  char buffer8[300];
  sprintf(buffer8,"%s < %i && %s < %i && %s < %i && %s<(1+%f) && %s>(1-%f) && %s==%i",amplitude.c_str(),amp_max,baselinerms.c_str(),baseline_cutoff,baselinerms2.c_str(),baseline_cutoff,base_ratio,ratio_threshold,base_ratio,ratio_threshold,channel.c_str(),channel_num);

  in_tree->Project("amp_fit_histo",amplitude.c_str(),buffer8,"goff",in_tree->GetEntries());

  TF1* amp_fit  = new TF1("amp_fit",MPEfit,amp_min,amp_max,7);
  amp_fit->SetNpx(1000);
  Double_t amp_par[7] = {ped_amp_estimate,2,amp_mu_pe_estimate,spe_amp_estimate,5,5,amp_mpe_norm_estimate};
  amp_fit->SetParNames("ped_mu","ped_sig","mu_pe","spe","pe_sig","gain","mpe_norm");
  amp_fit->SetParameters(amp_par);

  amp_fit->SetParLimits(0,ped_amp_estimate*0.9,ped_amp_estimate*1.1);
  amp_fit->SetParLimits(1,0,1000);
  amp_fit->SetParLimits(2,0,1000);
  amp_fit->SetParLimits(3,0,1000);
  amp_fit->SetParLimits(4,0,1000);
  amp_fit->SetParLimits(5,0,1000);
  amp_fit->SetParLimits(6,0,1e10);

  TCanvas *c1 = new TCanvas("c1","",10,10,800,600);
  c1->cd();
  c1->SetLogy();
  amp_fit_histo->Fit("amp_fit","R");
  
  double onePE_sig_  = TMath::Sqrt(TMath::Power(amp_fit->GetParameter(4),2)+TMath::Power(amp_fit->GetParameter(1),2));
  double onePE_norm_ = amp_fit->GetParameter(6)/(onePE_sig_*TMath::Sqrt(2*3.1415926))*TMath::Exp(-amp_fit->GetParameter(2))*amp_fit->GetParameter(2);
  double onePE_cent_ = amp_fit->GetParameter(3) + amp_fit->GetParameter(0);

  double twoPE_sig_  = TMath::Sqrt(2*TMath::Power(amp_fit->GetParameter(4),2)+TMath::Power(amp_fit->GetParameter(1),2));
  double twoPE_norm_ = amp_fit->GetParameter(6)/(twoPE_sig_*TMath::Sqrt(2*3.1415926))*TMath::Exp(-amp_fit->GetParameter(2))*TMath::Power(amp_fit->GetParameter(2),2)/2.0;
  double twoPE_cent_ = amp_fit->GetParameter(3) + amp_fit->GetParameter(0)+amp_fit->GetParameter(5);

  double threePE_sig_  = TMath::Sqrt(3*TMath::Power(amp_fit->GetParameter(4),2)+TMath::Power(amp_fit->GetParameter(1),2));
  double threePE_norm_ = amp_fit->GetParameter(6)/(threePE_sig_*TMath::Sqrt(2*3.1415926))*TMath::Exp(-amp_fit->GetParameter(2))*TMath::Power(amp_fit->GetParameter(2),3)/6.0;
  double threePE_cent_ = amp_fit->GetParameter(3) + amp_fit->GetParameter(0)+2*amp_fit->GetParameter(5);

  double fourPE_sig_  = TMath::Sqrt(4*TMath::Power(amp_fit->GetParameter(4),2)+TMath::Power(amp_fit->GetParameter(1),2));
  double fourPE_norm_ = amp_fit->GetParameter(6)/(fourPE_sig_*TMath::Sqrt(2*3.1415926))*TMath::Exp(-amp_fit->GetParameter(2))*TMath::Power(amp_fit->GetParameter(2),4)/24.0;
  double fourPE_cent_ = amp_fit->GetParameter(3) + amp_fit->GetParameter(0)+3*amp_fit->GetParameter(5);

  // The following segment plots the first several PE peaks overlaid with the final fit

  TF1 * subgauss1a = new TF1("subgauss1a","gaus",0,50);
  subgauss1a->SetLineColor(4);
  subgauss1a->SetParameters(onePE_norm_,onePE_cent_,onePE_sig_);
  subgauss1a->Draw("sames");

  TF1 * subgauss2a = new TF1("subgauss2a","gaus",0,50);
  subgauss2a->SetLineColor(4);
  subgauss2a->SetParameters(twoPE_norm_,twoPE_cent_,twoPE_sig_);
  subgauss2a->Draw("sames");

  TF1 * subgauss3a = new TF1("subgauss3a","gaus",0,50);
  subgauss3a->SetLineColor(4);
  subgauss3a->SetParameters(threePE_norm_,threePE_cent_,threePE_sig_);
  subgauss3a->Draw("sames");

  TF1 * subgauss4a = new TF1("subgauss4a","gaus",0,50);
  subgauss4a->SetLineColor(4);
  subgauss4a->SetParameters(fourPE_norm_,fourPE_cent_,fourPE_sig_);
  subgauss4a->Draw("sames");
 

  //----------------------------------------------------------------------------------------------------------------------------------------------------------
  // This second half basically does the same as the first half above, except fits the charge histogram


  TH1F* charge_fit_histo = new TH1F("charge_fit_histo","Charge Histogram With Fit",100,-100,500);

  charge_fit_histo->SetTitle("Charge Distribution With Fit");
  charge_fit_histo->GetYaxis()->SetTitle("Counts");
  charge_fit_histo->GetXaxis()->SetTitle("Charge (ADC*Ticks)");

  char buffer9[300];
  sprintf(buffer9,"%s < %i && %s < %i && %s < %i && %s < (1+%f) && %s >(1-%f) && %s==%i",charge.c_str(),charge_max,baselinerms.c_str(),baseline_cutoff,baselinerms2.c_str(),baseline_cutoff,base_ratio,ratio_threshold,base_ratio,ratio_threshold,channel.c_str(),channel_num);

  in_tree->Project("charge_fit_histo",charge.c_str(),buffer9,"goff",in_tree->GetEntries());

  TF1* charge_fit  = new TF1("charge_fit",MPEfit,charge_min,charge_max,7);
  charge_fit->SetNpx(1000);
  Double_t charge_par[7] = {ped_charge_estimate,2,charge_mu_pe_estimate,spe_charge_estimate,5,5,charge_mpe_norm_estimate};
  charge_fit->SetParNames("ped_mu","ped_sig","mu_pe","spe","pe_sig","gain","mpe_norm");
  charge_fit->SetParameters(charge_par);

  charge_fit->SetParLimits(0,ped_charge_estimate*0.9,ped_charge_estimate*1.1);
  charge_fit->SetParLimits(1,0,1000);
  charge_fit->SetParLimits(2,0,1);
  charge_fit->SetParLimits(3,0,1000);
  charge_fit->SetParLimits(4,0,1000);
  charge_fit->SetParLimits(5,0,1000);
  charge_fit->SetParLimits(6,0,1e10);

  TCanvas *c2 = new TCanvas("c2","",10,10,800,600);
  c2->cd();
  c2->SetLogy();
  charge_fit_histo->Fit("charge_fit","R");

  // The following segment plots the first several PE peaks overlaid with the final fit

  double onePE_sig  = TMath::Sqrt(TMath::Power(charge_fit->GetParameter(4),2)+TMath::Power(charge_fit->GetParameter(1),2));
  double onePE_norm = charge_fit->GetParameter(6)/(onePE_sig*TMath::Sqrt(2*3.1415926))*TMath::Exp(-charge_fit->GetParameter(2))*charge_fit->GetParameter(2);
  double onePE_cent = charge_fit->GetParameter(3) + charge_fit->GetParameter(0);

  double twoPE_sig  = TMath::Sqrt(2*TMath::Power(charge_fit->GetParameter(4),2)+TMath::Power(charge_fit->GetParameter(1),2));
  double twoPE_norm = charge_fit->GetParameter(6)/(twoPE_sig*TMath::Sqrt(2*3.1415926))*TMath::Exp(-charge_fit->GetParameter(2))*TMath::Power(charge_fit->GetParameter(2),2)/2.0;
  double twoPE_cent = charge_fit->GetParameter(3) + charge_fit->GetParameter(0)+charge_fit->GetParameter(5);

  double threePE_sig  = TMath::Sqrt(3*TMath::Power(charge_fit->GetParameter(4),2)+TMath::Power(charge_fit->GetParameter(1),2));
  double threePE_norm = charge_fit->GetParameter(6)/(threePE_sig*TMath::Sqrt(2*3.1415926))*TMath::Exp(-charge_fit->GetParameter(2))*TMath::Power(charge_fit->GetParameter(2),3)/6.0;
  double threePE_cent = charge_fit->GetParameter(3) + charge_fit->GetParameter(0)+2*charge_fit->GetParameter(5);

  double fourPE_sig  = TMath::Sqrt(4*TMath::Power(charge_fit->GetParameter(4),2)+TMath::Power(charge_fit->GetParameter(1),2));
  double fourPE_norm = charge_fit->GetParameter(6)/(fourPE_sig*TMath::Sqrt(2*3.1415926))*TMath::Exp(-charge_fit->GetParameter(2))*TMath::Power(charge_fit->GetParameter(2),4)/24.0;
  double fourPE_cent = charge_fit->GetParameter(3) + charge_fit->GetParameter(0)+3*charge_fit->GetParameter(5);


  TF1 * subgauss1 = new TF1("subgauss1","gaus",-100,500);
  subgauss1->SetLineColor(4);
  subgauss1->SetParameters(onePE_norm,onePE_cent,onePE_sig);
  subgauss1->Draw("sames");

  TF1 * subgauss2 = new TF1("subgauss2","gaus",-100,500);
  subgauss2->SetLineColor(4);
  subgauss2->SetParameters(twoPE_norm,twoPE_cent,twoPE_sig);
  subgauss2->Draw("sames");

  TF1 * subgauss3 = new TF1("subgauss3","gaus",-100,500);
  subgauss3->SetLineColor(4);
  subgauss3->SetParameters(threePE_norm,threePE_cent,threePE_sig);
  subgauss3->Draw("sames");

  TF1 * subgauss4 = new TF1("subgauss4","gaus",-100,500);
  subgauss4->SetLineColor(4);
  subgauss4->SetParameters(fourPE_norm,fourPE_cent,fourPE_sig);
  subgauss4->Draw("sames");

  //-----------------------------------------------------------------------------------------------------------------------------------------------------------


  return amp_fit->GetParameter(3)

}



 
