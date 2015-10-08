// A. Fiorentini 
// Macro to run in gmbrowser. 
//
// Show occupancy plots made by nearline in GMBrowser
//==================================================================================================//
//numib plots

void numibNhitsStrip(TCanvas *canvas, TFile *dat, TFile *ref)
{
  det_plots(canvas,dat,"h_nhits_ModVsStrip","Numib NHits/Gate for Strip (y) vs Module (x)");
}

void numibNhitsBar(TCanvas *canvas, TFile *dat, TFile *ref)
{
  det_plots(canvas,dat,"h_nhits_TowerVsBar","Numib NHits/Gate for Bar (y) vs Module (x)");
}

void numibNhitsCrate0NotFired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_nhits_ChannelMap_1_0_0","Numib NHits/Gate Disc Not Fired Crate 0");
}

void numibNhitsCrate0Fired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_nhits_ChannelMap_2_0_0","Numib NHits/Gate Disc Fired Crate 0");
}

void numibNhitsCrate1NotFired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_nhits_ChannelMap_1_0_1","Numib NHits/Gate Disc Not Fired Crate 1");
}

void numibNhitsCrate1Fired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_nhits_ChannelMap_2_0_1","Numib NHits/Gate Disc Fired Crate 1");
}

void numibAvgqhiStrip(TCanvas *canvas, TFile *dat, TFile *ref)
{
  det_plots(canvas,dat,"h_qhi_ModVsStrip","Numib AvgQhi for Strip (y) vs Module (x)");
}

void numibAvgqhiBar(TCanvas *canvas, TFile *dat, TFile *ref)
{
  det_plots(canvas,dat,"h_qhi_TowerVsBar","Numib AvgQhi for Bar (y) vs Module (x)");
}

void numibAvgqhiCrate0NotFired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_qhi_ChannelMap_1_0_0","Numib AvgQhi Disc Not Fired Crate 0");
}

void numibAvgqhiCrate0Fired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_qhi_ChannelMap_2_0_0","Numib AvgQhi Disc Fired Crate 0");
}

void numibAvgqhiCrate1NotFired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_qhi_ChannelMap_1_0_1","Numib AvgQhi Disc Not Fired Crate 1");
}

void numibAvgqhiCrate1Fired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_qhi_ChannelMap_2_0_1","Numib AvgQhi Disc Fired Crate 1");
}

//==================================================================================================//
//linjc plots

void linjcNhitsStrip(TCanvas *canvas, TFile *dat, TFile *ref)
{
  det_plots(canvas,dat,"h_nhits_ModVsStrip","Linjc NHits/Gate for Strip (y) vs Module (x)");
}

void linjcNhitsBar(TCanvas *canvas, TFile *dat, TFile *ref)
{
  det_plots(canvas,dat,"h_nhits_TowerVsBar","Linjc NHits/Gate for Bar (y) vs Module (x)");
}

void linjcNhitsCrate0NotFired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_nhits_ChannelMap_1_0_0","Linjc NHits/Gate Disc Not Fired Crate 0");
}

void linjcNhitsCrate0Fired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_nhits_ChannelMap_2_0_0","Linjc NHits/Gate Disc Fired Crate 0");
}

void linjcNhitsCrate1NotFired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_nhits_ChannelMap_1_0_1","Linjc NHits/Gate Disc Not Fired Crate 1");
}

void linjcNhitsCrate1Fired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_nhits_ChannelMap_2_0_1","Linjc NHits/Gate Disc Fired Crate 1");
}

void linjcAvgqhiStrip(TCanvas *canvas, TFile *dat, TFile *ref)
{
  det_plots(canvas,dat,"h_qhi_ModVsStrip","Linjc AvgQhi for Strip (y) vs Module (x)");
}

void linjcAvgqhiBar(TCanvas *canvas, TFile *dat, TFile *ref)
{
  det_plots(canvas,dat,"h_qhi_TowerVsBar","Linjc AvgQhi for Bar (y) vs Module (x)");
}

void linjcAvgqhiCrate0NotFired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_qhi_ChannelMap_1_0_0","Linjc AvgQhi Disc Not Fired Crate 0");
}

void linjcAvgqhiCrate0Fired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_qhi_ChannelMap_2_0_0","Linjc AvgQhi Disc Fired Crate 0");
}

void linjcAvgqhiCrate1NotFired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_qhi_ChannelMap_1_0_1","Linjc AvgQhi Disc Not Fired Crate 1");
}

void linjcAvgqhiCrate1Fired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_qhi_ChannelMap_2_0_1","Linjc AvgQhi Disc Fired Crate 1");
}

//==================================================================================================//
//pdstl plots

void pdstlNhitsStrip(TCanvas *canvas, TFile *dat, TFile *ref)
{
  det_plots(canvas,dat,"h_nhits_ModVsStrip","Pdstl NHits/Gate for Strip (y) vs Module (x)");
}

void pdstlNhitsBar(TCanvas *canvas, TFile *dat, TFile *ref)
{
  det_plots(canvas,dat,"h_nhits_TowerVsBar","Pdstl NHits/Gate for Bar (y) vs Module (x)");
}

void pdstlNhitsCrate0NotFired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_nhits_ChannelMap_1_0_0","Pdstl NHits/Gate Disc Not Fired Crate 0");
}

void pdstlNhitsCrate0Fired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_nhits_ChannelMap_2_0_0","Pdstl NHits/Gate Disc Fired Crate 0");
}

void pdstlNhitsCrate1NotFired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_nhits_ChannelMap_1_0_1","Pdstl NHits/Gate Disc Not Fired Crate 1");
}

void pdstlNhitsCrate1Fired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_nhits_ChannelMap_2_0_1","Pdstl NHits/Gate Disc Fired Crate 1");
}

void pdstlAvgqhiStrip(TCanvas *canvas, TFile *dat, TFile *ref)
{
  det_plots(canvas,dat,"h_qhi_ModVsStrip","Pdstl AvgQhi for Strip (y) vs Module (x)");
}

void pdstlAvgqhiBar(TCanvas *canvas, TFile *dat, TFile *ref)
{
  det_plots(canvas,dat,"h_qhi_TowerVsBar","Pdstl AvgQhi for Bar (y) vs Module (x)");
}

void pdstlAvgqhiCrate0NotFired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_qhi_ChannelMap_1_0_0","Pdstl AvgQhi Disc Not Fired Crate 0");
}

void pdstlAvgqhiCrate0Fired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_qhi_ChannelMap_2_0_0","Pdstl AvgQhi Disc Fired Crate 0");
}

void pdstlAvgqhiCrate1NotFired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_qhi_ChannelMap_1_0_1","Pdstl AvgQhi Disc Not Fired Crate 1");
}

void pdstlAvgqhiCrate1Fired(TCanvas *canvas, TFile *dat, TFile *ref)
{
  elec_plots(canvas,dat,"h_qhi_ChannelMap_2_0_1","Pdstl AvgQhi Disc Fired Crate 1");
}

//==================================================================================================//
//drawing funtions

void elec_plots(TCanvas *canvas, TFile *dat, string histo_name, string title)
{
  // check if files are there (and open)
  if(!dat->IsOpen()){
    cerr << "Can not open file " << file << endl;
    return;
  }

  // get info from title to plot histos
  double min = 0;
  double max = 0;
  bool scale = false;
  string ngates_histo = "";//name of histogram containing number of gates 
  string histo_folder = "";//name of folder containing occupancy plots 
  if( !get_info(title,min,max,scale,histo_folder,ngates_histo) ) return;

  // get number of gates from this histogram
  int ngates = 0;
  if(scale){
    TH1D* h = (TH1D*)dat->Get( ngates_histo.c_str() );
    if( !h ) return;
    ngates = h->GetEntries();
    if( !ngates ) return;
  }

  TH2D* histo;

  canvas->Divide(2,4);

  char temp[100];
  for(int pad=1; pad<=8; pad += 2){
    //****************
    //*** PAD #pad ***
    //****************
    canvas->cd(pad);
    sprintf(temp,"%s/%s_%i",histo_folder.c_str(),histo_name.c_str(),pad);
    histo = (TH2D*)dat->Get(temp);
    if( !histo ) return;
    sprintf(temp,"%s Croc %i;board;channel",title.c_str(),pad);
    histo->SetTitle(temp);
    if( scale ) histo->Scale(1./ngates);
    histo->SetStats(kFALSE);
    histo->SetTitleOffset(0.6,"X");
    histo->SetTitleSize(0.075,"X");
    histo->SetTitleOffset(0.4,"Y");
    histo->SetTitleSize(0.075,"Y");
    histo->SetMinimum(min);
    histo->SetMaximum(max);
    histo->Draw("colz");
  }

   for(int pad=2; pad<=8; pad += 2){
    //****************
    //*** PAD #pad ***
    //****************
    canvas->cd(pad);
    sprintf(temp,"%s/%s_%i",histo_folder.c_str(),histo_name.c_str(),pad);
    histo = (TH2D*)dat->Get(temp);
    if( !histo ) return;
    sprintf(temp,"%s Croc %i;board;channel",title.c_str(),pad);
    histo->SetTitle(temp);
    if( scale ) histo->Scale(1./ngates);
    histo->SetStats(kFALSE);
    histo->SetTitleOffset(0.6,"X");
    histo->SetTitleSize(0.075,"X");
    histo->SetTitleOffset(0.4,"Y");
    histo->SetTitleSize(0.075,"Y");
    histo->SetMinimum(min);
    histo->SetMaximum(max);
    histo->Draw("colz");
  }

  canvas->Update();
}

void det_plots(TCanvas *canvas, TFile *dat, string histo_name, string title)
{
  // check if files are there (and open)
  if (!dat->IsOpen())
      cerr << "Can not open file " << file << endl;

  // get info from title to plot histos
  double min = 0;
  double max = 0;
  bool scale = false;
  string ngates_histo = "";//name of histogram containing number of gates 
  string histo_folder = "";//name of folder containing occupancy plots 
  if( !get_info(title,min,max,scale,histo_folder,ngates_histo) ) return;

  TH2D* histo;

  // get number of gates from this histogram
  int ngates = 0;
  if(scale){
    TH1D* h = (TH1D*)dat->Get( ngates_histo.c_str() );
    if( !h ) return;
    ngates = h->GetEntries();
    if( !ngates ) return;
  }

  canvas->Divide(1,2);

  char temp[100];
  for(int pad=1; pad<=2; pad++){
    //****************
    //*** PAD #pad ***
    //****************
    canvas->cd(pad);
    sprintf(temp,"%s/%s_%i",histo_folder.c_str(),histo_name.c_str(),pad); 
    histo = (TH2D*)dat->Get(temp);
    if( !histo ) return;
    histo->SetTitle(title.c_str());
    histo->GetXaxis()->CenterTitle();
    histo->GetYaxis()->CenterTitle();
    if( scale ) histo->Scale(1./ngates);
    histo->SetStats(kFALSE);
    histo->SetMinimum(min);
    histo->SetMaximum(max);
    histo->Draw("colz");
  }

  canvas->Update();
}

bool get_info(string title_s,double& min,double& max, bool& scale, string& histo_folder, string& ngates_histo )
{

  TString title = title_s;

  // scale NHits histograms only
  bool nhits = title.Contains("NHits") ? true : false;
  bool avgqhi = title.Contains("AvgQhi") ? true : false;
  bool scale = nhits;

  if( title.Contains("Numib") ){

    if(nhits) min = 0.0;
    if(nhits) max = 0.25; // raise to .25
    if(avgqhi) min = 0.0;
    if(avgqhi) max = 1500.0;
    ngates_histo = "NumibDAQHeaderHistos/h_nADCFrames";
    histo_folder = "NumibOccupancyPlots.MnvDetPlottingTool";

  } else if( title.Contains("Linjc") ) {

    if(nhits) min = 0.0;
    if(nhits) max = 1.0;
    if(avgqhi) min = 0.0;
    if(avgqhi) max = 1500.0;
    ngates_histo = "LinjcDAQHeaderHistos/h_nADCFrames";
    histo_folder = "LinjcOccupancyPlots.MnvDetPlottingTool";

  } else if( title.Contains("Pdstl") ) {

    if(nhits) min = 0.0;
    if(nhits) max = 1.5;
    if(avgqhi) min = 0.0;
    if(avgqhi) max = 600.0;
    ngates_histo = "PdstlDAQHeaderHistos/h_nADCFrames";
    histo_folder = "PdstlOccupancyPlots.MnvDetPlottingTool";

  } else return false;

  return true;

}
