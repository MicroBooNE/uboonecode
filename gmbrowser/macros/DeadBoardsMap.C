void DeadBoardsMap(TCanvas *c, TFile *data, TFile *ref){

  // Get and add board maps for linjc and numib.
  TH2D* boardMapCrate0Numib = data->Get("NumibOccupancyPlots.MnvDetPlottingTool/h_nhits_BoardMap_2_0_0");
  TH2D* boardMapCrate1Numib = data->Get("NumibOccupancyPlots.MnvDetPlottingTool/h_nhits_BoardMap_2_0_1");
  TH2D* boardMapCrate0Linjc = data->Get("LinjcOccupancyPlots.MnvDetPlottingTool/h_nhits_BoardMap_2_0_0");
  TH2D* boardMapCrate1Linjc = data->Get("LinjcOccupancyPlots.MnvDetPlottingTool/h_nhits_BoardMap_2_0_1");

  TH2D* boardMapCrate0;
  if( boardMapCrate0Numib && boardMapCrate0Linjc ){
    boardMapCrate0 = (TH2D*)boardMapCrate0Numib->Clone("boardMapCrate0Total");
    boardMapCrate0->Add(boardMapCrate0Numib,boardMapCrate0Linjc);
  } else if( boardMapCrate0Numib ){
    boardMapCrate0 = (TH2D*)boardMapCrate0Numib->Clone("boardMapCrate0Total");
  } else if( boardMapCrate0Linjc ){
    boardMapCrate0 = (TH2D*)boardMapCrate0Linjc->Clone("boardMapCrate0Total");
  }

  TH2D* boardMapCrate1;
  if( boardMapCrate1Numib && boardMapCrate1Linjc ){
    boardMapCrate1 = (TH2D*)boardMapCrate1Numib->Clone("boardMapCrate1Total");
    boardMapCrate1->Add(boardMapCrate1Numib,boardMapCrate1Linjc);
  } else if( boardMapCrate1Numib ){
    boardMapCrate1 = (TH2D*)boardMapCrate1Numib->Clone("boardMapCrate1Total");
  } else if( boardMapCrate1Linjc ){
    boardMapCrate1 = (TH2D*)boardMapCrate1Linjc->Clone("boardMapCrate1Total");
  }

  // Get histogram with number of events.
  TProfile* h_nevents_numib = data->Get("NumibDAQHeaderHistos/Run");
  int nevents_numib = (bool)h_nevents_numib ? h_nevents_numib->GetEntries() : 0;

  TProfile* h_nevents_linjc = data->Get("LinjcDAQHeaderHistos/Run");
  int nevents_linjc = (bool)h_nevents_linjc ? h_nevents_linjc->GetEntries() : 0;

  // create 3 pads inside canvas
  TPad* pad;
  double x1,y1,x2,y2;

  x1=0.01; y1=0.1+0.01; x2=0.5-0.01; y2=1-0.01;
  pad = new TPad("c1","c1",x1,y1,x2,y2,0);
  pad->SetNumber(1);
  pad->Draw();

  x1=0.5+0.01; y1=0.1+0.01; x2=1-0.01; y2=1-0.01;
  pad = new TPad("c2","c2",x1,y1,x2,y2,0);
  pad->SetNumber(2);
  pad->Draw();

  x1=0.01; y1=0.01; x2=1-0.01; y2=0.1-0.01;
  pad = new TPad("c3","c3",x1,y1,x2,y2,0);
  pad->SetNumber(3);
  pad->Draw();

  // Plot Crate 0
  c->cd(1);
  int nBoardsInCroc0[32] = {10,10,10,6,10,10,9,5,10,10,10,10,9,9,9,9,10,10,10,10,9,9,9,9,10,10,10,10,9,9,9,9}; 
  if( boardMapCrate0 ){
    boardMapCrate0->SetTitle("Dead Boards Map (Crate 0);Board Number;Croc Number");
    refillMap( boardMapCrate0, nBoardsInCroc0 );
    drawMap( boardMapCrate0 );
  }

  // Plot Crate 1
  c->cd(2);
  int nBoardsInCroc1[32] = {10,10,10,10,9,9,9,9,10,10,6,6,9,9,5,5,6,6,6,2,5,5,5,0,10,10,10,10,0,0,0,0}; 
  if( boardMapCrate1 ){
    boardMapCrate1->SetTitle("Dead Boards Map (Crate 1);Board Number;Croc Number");
    refillMap( boardMapCrate1, nBoardsInCroc1 );
    drawMap( boardMapCrate1 );
  }

  // Plot Number of events
  c->cd(3);
  TLatex* text_numib = new TLatex(0.1,0.9,Form("Number of Numi Beam Gates = %d",nevents_numib));
  text_numib->SetTextAlign(13);
  text_numib->SetTextSize(0.4);
  text_numib->Draw();

  TLatex* text_linjc = new TLatex(0.1,0.1,Form("Number of Light Injection Gates = %d",nevents_linjc));
  text_linjc->SetTextAlign(11);
  text_linjc->SetTextSize(0.4);
  text_linjc->Draw();

}

void refillMap( TH2D* h, int* nBoardsInCroc ){

  //refill histograms
  int nbinsx = h->GetNbinsX();
  int nbinsy = h->GetNbinsY();
  for(int binx=1; binx<=nbinsx; binx++){
    for(int biny=1; biny<=nbinsy; biny++){

      //avoid not conected boards
      int bin = h->GetBin(binx,biny);
      if( binx>nBoardsInCroc[biny-1] ){
        h->SetBinContent(bin,0);
        continue;
      }

      //if board has no hits give it a value of 1 (red)
      if( h->GetBinContent(bin)==0 ) h->SetBinContent(bin,1);
      //if board has hits give it a value of 0.1 (blue) 
      else h->SetBinContent(bin,0.15);

    }
  }

}

void drawMap( TH2D* h ){

  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetNdivisions(12,kFALSE);
  h->GetXaxis()->CenterLabels();

  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetLabelColor(0);
  h->GetYaxis()->SetNdivisions(408,kFALSE);

  h->GetZaxis()->SetRangeUser(0,1);

  h->SetStats(0);
  h->Draw("col");

  //new axis for croc numbers
  TF1 *f=new TF1("f","x/8",0,8);
  TGaxis *axis = new TGaxis(1,0,1,32,1,9,408,"B");
  axis->CenterLabels();
  axis->Draw();

  gPad->SetGrid();

}

