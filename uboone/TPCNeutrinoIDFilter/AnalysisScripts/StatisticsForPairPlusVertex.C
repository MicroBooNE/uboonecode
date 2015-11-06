//**************************************
// This root Macro will use analysis tree files and search for trios of two tracks and one reco vertex. It evaluates the distances between these objects and calculates passing rates based on a maximum distance requirement in order to identify neutrino events. It runs on cosmics, bnb files, bnb + cosmics, and experimental data.
// Author: aschu@fnal.gov
// ************************************

using namespace std;

double GetDist(double x1, double y1, double z1, double x2, double y2, double z2) {
   return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
}

int StatisticsForPairPlusVertex() {

   TChain *tree = new TChain("analysistree/anatree");

   //Set the path to an analysistree file here.

   //tree -> Add("/Users/aschu/Desktop/ccinclusive/prodgenie_bnb_nu_cosmic_uboone/anahist.root"); // this is the same pre-MCC6.1, but containes pandora vertex as well
   //tree -> Add("/Users/aschu/Desktop/ccinclusive/prodgenie_bnb_nu_uboone/anahist.root"); //
   //tree -> Add("/Users/aschu/Desktop/ccinclusive/prodcosmics_uboone/anahist.root");
   //tree -> Add("/Users/aschu/Desktop/FirstNeutrinoInteraction/anahist-run1712.root");
   //tree -> Add("/Users/aschu/Desktop/FirstNeutrinoInteraction/anahist-vertices-bnbonly.root");
   //tree -> Add("/uboone/data/users/aschu/CosmicTagging/anahist-vertices-cosmiconly.root");
   //tree -> Add("/uboone/data/users/aschu/CosmicTagging/anahist-cosmic-runset_0.root");
   //tree -> Add("/uboone/data/users/aschu/CosmicTagging/anahist-cosmic-runset_1.root");
   tree -> Add("/uboone/data/users/aschu/CosmicTagging/anahist-bnbonly-70kV-1us-v04_21_01.root");

   //********  Fill this!  *************************************************************************************
   // it will make event rates scale to 1h running at 5 Hz.
   //
   int mode = 0;// 0 = BNB only, 1 = Cosmic Only Simulation, 2 = BNB + cosmic simulation, 3 = Cosmic DATA
   //
   //***********************************************************************************************************

   const int maxentries = 50000;
   const int maxtracks = 100;
   const int maxvtx = 300;
   const int maxnu = 10;

   Int_t           event;
   Int_t           ccnc_truth;

   Int_t           ntracks_reco;
   Float_t         trkstartx[maxtracks];
   Float_t         trkstarty[maxtracks];
   Float_t         trkstartz[maxtracks];
   Float_t         trkendx[maxtracks];
   Float_t         trkendy[maxtracks];
   Float_t         trkendz[maxtracks];
   Float_t         trkcosmicscore_tagger[maxtracks];
   Float_t         trkcosmicscore_flashmatch[maxtracks];

   Short_t         nvtx;
   Float_t         vtxx[maxvtx];
   Float_t         vtxy[maxvtx];
   Float_t         vtxz[maxvtx];

   //Note that you have to change the branch names if you want to test different reconstruction algorithms.
   tree -> SetBranchAddress("event", &event);
   tree -> SetBranchAddress("ntracks_trackkalmanhit", &ntracks_reco);
   tree -> SetBranchAddress("trkstartx_trackkalmanhit", trkstartx);
   tree -> SetBranchAddress("trkstarty_trackkalmanhit", trkstarty);
   tree -> SetBranchAddress("trkstartz_trackkalmanhit", trkstartz);
   tree -> SetBranchAddress("trkendx_trackkalmanhit", trkendx);
   tree -> SetBranchAddress("trkendy_trackkalmanhit", trkendy);
   tree -> SetBranchAddress("trkendz_trackkalmanhit", trkendz);
   tree -> SetBranchAddress("trkcosmicscore_tagger_trackkalmanhit", trkcosmicscore_tagger);
   tree -> SetBranchAddress("trkcosmicscore_flashmatch_trackkalmanhit", trkcosmicscore_flashmatch);
   tree -> SetBranchAddress("nvtx_pandoraNu", &nvtx);
   tree -> SetBranchAddress("vtxx_pandoraNu", vtxx);
   tree -> SetBranchAddress("vtxy_pandoraNu", vtxy);
   tree -> SetBranchAddress("vtxz_pandoraNu", vtxz);

   //Loop through MC file
   int Size = tree -> GetEntries();
   cout << "number of events used is: " << Size << endl;

   TH1D *h_triangles = new TH1D("h_triangles", "Max distance for each triangle found", 2000, 0, 1000);
   TH1D *h_event = new TH1D("h_event", "Max distance for the best triangle in the event", 2000, 0, 1000);

   double dist = 0;
   double dist1 = 0;
   double dist2 = 0;
   double diststart = 0;
   double distend = 0;
   double temp = 0;
   double maxdist = 0;
   double mindist = 0;

   for(int i = 0; i < Size; i++) {
      if(i%1000 == 0) cout << "\t... " << i << endl;

      mindist = 100000;
      tree -> GetEntry(i);

      //loop over all vertices
      if(nvtx > 0) {
         for(int v = 0; v < nvtx; v++) {
           //loop over tracks
            for(int j = 0; j < ntracks_reco; j++) {
               if((trkcosmicscore_tagger[j] < 0.4)) {
                  diststart = GetDist(trkstartx[j], trkstarty[j], trkstartz[j], vtxx[v], vtxy[v], vtxz[v]);
                  distend = GetDist(trkendx[j], trkendy[j], trkendz[j], vtxx[v], vtxy[v], vtxz[v]);
                  if(distend < diststart) {
                     temp = trkstartx[j]; trkstartx[j] = trkendx[j]; trkendx[j] = temp;
                     temp = trkstarty[j]; trkstarty[j] = trkendy[j]; trkendy[j] = temp;
                     temp = trkstartz[j]; trkstartz[j] = trkendz[j]; trkendz[j] = temp;
		     dist1 = distend;
                  }
                  else dist1 = diststart;
                  for(int k = j+1; k < ntracks_reco; k++) {
                     if((trkcosmicscore_tagger[k] < 0.4)) {
                        diststart = GetDist(trkstartx[k], trkstarty[k], trkstartz[k], vtxx[v], vtxy[v], vtxz[v]);
                        distend = GetDist(trkendx[k], trkendy[k], trkendz[k], vtxx[v], vtxy[v], vtxz[v]);
                        if(distend < diststart) {
                           temp = trkstartx[k]; trkstartx[k] = trkendx[k]; trkendx[k] = temp;
                           temp = trkstarty[k]; trkstarty[k] = trkendy[k]; trkendy[k] = temp;
                           temp = trkstartz[k]; trkstartz[k] = trkendz[k]; trkendz[k] = temp;
                           dist2 = distend;
                        }
                        else dist2 = diststart;
                        dist = GetDist(trkstartx[j], trkstarty[j], trkstartz[j], trkstartx[k], trkstarty[k], trkstartz[k]);

                        //get the maximum distance of the triangle
                        if(dist >= dist1 && dist >= dist2) maxdist = dist;
                        else if(dist1 >= dist && dist1 >= dist2) maxdist = dist1;
                        else if(dist2 >= dist && dist2 >= dist1) maxdist = dist2;
                        h_triangles -> Fill(maxdist);

                        //get the best matching triangle in the event
                        if(maxdist < mindist) mindist = maxdist;
                     }//if neutrino tag
                  }//loop over inner track
               }//if neutrino tag
            }//loop over outer track
         }//end loop over vertices
      }//end vertex if

      if(mindist < 100000) h_event -> Fill(mindist);
   }

   //Scaling to events per hour
   double scalefactor = 1;
   if(mode == 0) scalefactor = 20000/((double) Size) / 5.3e19 * 4e12 * 5 * 3600;
   if(mode == 1) scalefactor = 3600 * 5 / ((double) Size);
   if(mode == 2) scalefactor = 20000/((double) Size) / 5.3e19 * 4e12 * 5 * 3600;
   if(mode == 3) scalefactor = 3600 * 5 / ((double) Size);
   cout << "The scale factor is: " << scalefactor << endl;

   TCanvas *c2 = new TCanvas("c2", "c2", 800, 800);
   c2 -> SetGridx(1);
   c2 -> SetGridy(1);
   c2 -> SetTicks(1, 1);
   c2 -> SetLogy(1);

   h_triangles -> Scale(scalefactor);
   (h_triangles -> GetYaxis()) -> SetRangeUser(0.1, 10);
   (h_triangles -> GetYaxis()) -> SetTitle("triangles per hour @ 5Hz");
   (h_triangles -> GetYaxis()) -> SetTitleOffset(1.5);
   (h_triangles -> GetXaxis()) -> SetTitle("#Delta [cm]");
   (h_triangles -> GetXaxis()) -> SetTitleOffset(1.2);
   h_triangles -> SetLineColor(2);
   h_triangles -> SetStats(0);
   h_triangles -> SetLineWidth(3);
   h_triangles -> Draw();

   TCanvas *c3 = new TCanvas("c3", "c3", 800, 800);
   c3 -> SetGridx(1);
   c3 -> SetGridy(1);
   c3 -> SetTicks(1, 1);
   c3 -> SetLogy(1);

   h_event -> Scale(scalefactor);
   (h_event -> GetYaxis()) -> SetRangeUser(0.1, 10);
   (h_event -> GetXaxis()) -> SetTitle("#Delta [cm]");
   (h_event -> GetXaxis()) -> SetTitleOffset(1.2);
   (h_event -> GetYaxis()) -> SetTitle("# events per hour @ 5Hz");
   (h_event -> GetYaxis()) -> SetTitleOffset(1.5);
   h_event -> SetLineColor(4);
   h_event -> SetStats(0);
   h_event -> SetLineWidth(3);
   h_event -> Draw();

   cout << "0.5\t1\t1.5\t2\t2.5\t3\t3.5\t4\t4.5\t5" << endl;
   double n05 = h_event -> GetBinContent(1);
   double n10 = h_event -> GetBinContent(2) + n05;
   double n15 = h_event -> GetBinContent(3) + n10;
   double n20 = h_event -> GetBinContent(4) + n15;
   double n25 = h_event -> GetBinContent(5) + n20;
   double n30 = h_event -> GetBinContent(6) + n25;
   double n35 = h_event -> GetBinContent(7) + n30;
   double n40 = h_event -> GetBinContent(8) + n35;
   double n45 = h_event -> GetBinContent(9) + n40;
   double n50 = h_event -> GetBinContent(10) + n45;
   cout << n05 << "\t" << n10 << "\t" << n15 << "\t" << n20 << "\t" << n25 << "\t" << n30 << "\t" << n35 << "\t" << n40 << "\t" << n45 << "\t" << n50 << endl;

   return 0;

}
