//unelegant way to plot veto histograms
void PlotVetoOne(TCanvas *c1, TFile *data, TFile *ref)
{
  // Get histograms.
  TH1* PairOneBoardOne = data->Get("VetoPlotAlg/HitPerPixDiscYPairOne_0");
  TH1* PairOneBoardTwo = data->Get("VetoPlotAlg/HitPerPixDiscYPairOne_1");
  TH1* PairTwoBoardTwo = data->Get("VetoPlotAlg/HitPerPixDiscYPairTwo_1");
  TH1* TripPushes      = data->Get("VetoPlotAlg/TripPushes");
  
  c1->Divide(2, 2);

  // Plot first histogram.
  c1->cd(1);
  PairOneBoardOne->GetXaxis()->SetTitle("Channel Number");
  PairOneBoardOne->GetYaxis()->SetTitle("Number of Hits/Gate");
  PairOneBoardOne->Draw("colrz");
  PairOneBoardOne->SetStats(kFALSE);
  setPal(PairOneBoardOne);
  
  c1->cd(2);
  PairOneBoardTwo->GetXaxis()->SetTitle("Channel Number");
  PairOneBoardTwo->GetYaxis()->SetTitle("Number of Hits/Gate");
  PairOneBoardTwo->Draw("colrz");
  PairOneBoardTwo->SetStats(kFALSE);
  setPal(PairOneBoardTwo);
  
  c1->cd(3);
  PairTwoBoardTwo->GetXaxis()->SetTitle("Channel Number");
  PairTwoBoardTwo->GetYaxis()->SetTitle("Number of Hits/Gate");
  PairTwoBoardTwo->Draw("colrz");
  PairTwoBoardTwo->SetStats(kFALSE);
  setPal(PairTwoBoardTwo);
  
  
  c1->cd(4);
  for(int i=1; i<5; ++i){
     TripPushes->GetXaxis()->SetBinLabel(i, Form("Pair %d", i - 1));
  }
  for(int i=1; i<11; ++i){
     TripPushes->GetYaxis()->SetBinLabel(i, Form("%d", i - 1));
  }
  
  TripPushes->GetXaxis()->SetTitle("Trip pair. (0 and 1: FEB 1, 2 and 3: FEB 2)");
  TripPushes->GetYaxis()->SetTitle("Number of Pushes/Gate");
  TripPushes->SetMinimum(0);
  TripPushes->SetMaximum(1000);
  gPad->SetLogz(1);
  TripPushes->SetStats(kFALSE);
  TripPushes->Draw("colrz");
  setPal(TripPushes);

} //end of PlotVetoOne

void PlotVetoTwo(TCanvas *c1, TFile *data, TFile *ref)
{
  // Get histograms.
  TH1* PhysMapWallOne = data->Get("VetoPlotAlg/HitPerTube_0");
  TH1* PhysMapWallTwo = data->Get("VetoPlotAlg/HitPerTube_1");
  for(int i=2; i < 8; ++i){
    PhysMapWallTwo->GetYaxis()->SetBinLabel(i, Form("%d", 8 - i));
    PhysMapWallOne->GetYaxis()->SetBinLabel(i, Form("%d", 8 - i));
    }
  PhysMapWallTwo->GetXaxis()->SetBinLabel(2, "1");
  PhysMapWallTwo->GetXaxis()->SetBinLabel(3, "2");
  PhysMapWallOne->GetXaxis()->SetBinLabel(2, "1");
  PhysMapWallOne->GetXaxis()->SetBinLabel(3, "2");    
  c1->Divide(2, 1);

  // Plot first histogram.
  c1->cd(1);
  PhysMapWallTwo->GetXaxis()->SetTitle("PMT Number (1 = W 2 = E )");
  PhysMapWallTwo->GetYaxis()->SetTitle("Paddle Number (1 = Highest)");
  PhysMapWallTwo->SetStats(kFALSE);
  
  PhysMapWallTwo->Draw("colrz");
  
  setPal(PhysMapWallTwo);
  
  c1->cd(2);
  PhysMapWallOne->GetXaxis()->SetTitle("PMT Number (1 = W 2 = E )");
  PhysMapWallOne->GetYaxis()->SetTitle("Paddle Number (1 = Highest)");
  PhysMapWallOne->SetStats(kFALSE);
  PhysMapWallOne->Draw("colrz");
  setPal(PhysMapWallOne);
  

} //end of PlotVetoTwo

void setPal( TH1 *hist, double label_size = 0.03){
 
 gPad->Update();
 TPaletteAxis *p = (TPaletteAxis*)hist->GetListOfFunctions()->FindObject("palette");
 if(p == NULL)
    return;
 p->SetLabelSize(label_size);
 gPad->Modified();
 gPad->Update();

}

