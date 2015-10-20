//plot high voltage histograms
void plotHighVoltage(TCanvas *c1, TFile *data, TFile *ref)
{

  if(!data) return;

  // Get histograms.
  TH1* hv          = data->Get("ElectronicsMapAlg/PMThv");
  TH1* hvresid     = data->Get("ElectronicsMapAlg/PMThvResid");
  TH1* hvdisabled1 = data->Get("ElectronicsMapAlg/PMThvDisabled_crate0");
  TH1* hvdisabled2 = data->Get("ElectronicsMapAlg/PMThvDisabled_crate1");

  // Check if HV disabled PMTs exist
  bool hv_disabled = false;
  if(hvdisabled1 || hvdisabled2){
    hv_disabled = true;
  }
 
  c1->Divide(2, 1);

  // Plot HV histogram.
  c1->cd(1);
  if(hv){
    // fix Y axis limits
    hv->SetMinimum(1);
    hv->SetMaximum(100000);

    // turn plot red if HV disabled PMTs were found
    if(hv_disabled) 
      hv->SetFillColor(kRed);

    hv->GetXaxis()->CenterTitle();
    hv->SetStats(kFALSE);
    hv->SetLineWidth(2);
    hv->Draw("HIST");
    gPad->SetLogy();
  }
 
  // Plot HV residuals histogram.
  c1->cd(2);
  if(hvresid){
    // fix Y axis limits
    hvresid->SetMinimum(1);
    hvresid->SetMaximum(100000);

    // turn plot red if HV disabled PMTs were found
    if(hv_disabled) 
      hvresid->SetFillColor(kRed);

    hvresid->GetXaxis()->CenterTitle();
    hvresid->SetStats(kFALSE);
    hvresid->SetLineWidth(2);
    hvresid->Draw("HIST");
    gPad->SetLogy();
  }
  
} //end of plotHighVoltage

