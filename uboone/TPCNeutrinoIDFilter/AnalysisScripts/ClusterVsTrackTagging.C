using namespace std;
#include <iostream>
#include <string>
#include "TChain.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TAxis.h"

int ClusterVsTrackTagging() {

   TChain *tree = new TChain("analysistree/anatree");

   //tree -> Add("/uboone/data/users/aschu/CosmicTagging/CosmicTaggingComparison/anahist-cosmic.root");
   tree -> Add("/uboone/data/users/aschu/CosmicTagging/CosmicTaggingComparison/anahist-data-v04_22_00-fixed-RunSet0.root");

   //tree -> Add("/pnfs/uboone/scratch/users/bcarls/neutrino_id/v04_21_01/mergeana/prod_data_run_1712_uboone/anahist.root");
   //tree -> Add("/uboone/data/users/aschu/CosmicTagging/anahist-cosmic-1712-1713-1715-1716-1717.root");
   //tree -> Add("/uboone/data/users/aschu/CosmicTagging/anahist-cosmic-runset_0.root");
 
   const int maxtracks = 100;
   const int maxvtx = 500;
   const int kMaxClusters = 1200;

   Int_t           run;
   Int_t           subrun;
   Int_t           event;

   Short_t nclusters;     //number of clusters in a  given event
   Short_t cluncosmictags_tagger[kMaxClusters];      //No. of cosmic tags associated to this cluster
   Float_t clucosmicscore_tagger[kMaxClusters];      //Cosmic score associated to this cluster. In the case of more than one tag, the first one is associated.
   Short_t clucosmictype_tagger[kMaxClusters];       //Cosmic tag type for this cluster.
   Short_t cluster_StartWire[kMaxClusters]; //wire coordinate of the start of the cluster 
   Short_t cluster_StartTick[kMaxClusters];  //tick coordinate of the start of the cluster in time ticks
   Short_t cluster_EndWire[kMaxClusters];  //wire coordinate of the end of the cluster
   Short_t cluster_EndTick[kMaxClusters];  //tick coordinate of the end of the cluster in time ticks
   Short_t clusterView[kMaxClusters];   //which plane this cluster belongs to        

   Short_t         ntracks_reco;
   Float_t         trkke[maxtracks][3];   //[ntracks_trackkalmanhit]
   Float_t         trklen[maxtracks];   //[ntracks_trackkalmanhit]
   Short_t         ntrkhits[maxtracks][3];   //[ntracks_trackkalmanhit]
   Float_t         trkefftruth[maxtracks][3];
   Short_t         trkorigin[maxtracks][3];
   Short_t         trkpidbestplane[maxtracks];
   Float_t         trkstartx[maxtracks];
   Float_t         trkstarty[maxtracks];
   Float_t         trkstartz[maxtracks];
   Float_t         trkendx[maxtracks];
   Float_t         trkendy[maxtracks];
   Float_t         trkendz[maxtracks];
   Float_t         trktheta[maxtracks];
   Float_t         trkphi[maxtracks];
   Float_t         trkpidpida[maxtracks][3];
   Float_t         trkcosmicscore_tagger[maxtracks];
   Float_t         trkcosmicscore_flashmatch[maxtracks];

   tree -> SetBranchAddress("run", &run);
   tree -> SetBranchAddress("subrun", &subrun);
   tree -> SetBranchAddress("event", &event);
   tree -> SetBranchAddress("cluncosmictags_tagger", cluncosmictags_tagger);
   tree -> SetBranchAddress("clucosmicscore_tagger", clucosmicscore_tagger);
   tree -> SetBranchAddress("nclusters", &nclusters);
   tree -> SetBranchAddress("cluster_StartWire", cluster_StartWire);
   tree -> SetBranchAddress("cluster_EndWire", cluster_EndWire);
   tree -> SetBranchAddress("cluster_StartTick", cluster_StartTick);
   tree -> SetBranchAddress("cluster_EndTick", cluster_EndTick);
   tree -> SetBranchAddress("clusterView", clusterView);
   tree -> SetBranchAddress("trkstartx_trackkalmanhit", trkstartx);
   tree -> SetBranchAddress("trkstarty_trackkalmanhit", trkstarty);
   tree -> SetBranchAddress("trkstartz_trackkalmanhit", trkstartz);
   tree -> SetBranchAddress("trkendx_trackkalmanhit", trkendx);
   tree -> SetBranchAddress("trkendy_trackkalmanhit", trkendy);
   tree -> SetBranchAddress("trkendz_trackkalmanhit", trkendz);
   tree -> SetBranchAddress("ntracks_trackkalmanhit", &ntracks_reco);
   tree -> SetBranchAddress("trkke_trackkalmanhit", trkke);
   tree -> SetBranchAddress("trklen_trackkalmanhit", trklen);
   tree -> SetBranchAddress("ntrkhits_trackkalmanhit", ntrkhits);
   tree -> SetBranchAddress("trktheta_trackkalmanhit", trktheta);
   tree -> SetBranchAddress("trkphi_trackkalmanhit", trkphi);
   tree -> SetBranchAddress("trkpidpida_trackkalmanhit", trkpidpida);
   tree -> SetBranchAddress("trkpidbestplane_trackkalmanhit", trkpidbestplane);
   tree -> SetBranchAddress("trkorigin_trackkalmanhit", trkorigin);
   tree -> SetBranchAddress("trkcosmicscore_tagger_trackkalmanhit", trkcosmicscore_tagger);
   tree -> SetBranchAddress("trkcosmicscore_flashmatch_trackkalmanhit", trkcosmicscore_flashmatch);

   TCanvas *c = new TCanvas("c", "c", 600, 600);
   c -> SetGridx(1);
   c -> SetGridy(1);
   c -> SetTicks(1, 1);

   TH1D *h_cluster = new TH1D("cluster", "", 11, -0.05, 1.05);
   TH1D *h_track = new TH1D("track", "", 11, -0.05, 1.05);

   TH1D *h_cluster_ev = new TH1D("cluster_ev", "", 11, -0.05, 1.05);
   TH1D *h_track_ev = new TH1D("track_ev", "", 11, -0.05, 1.05);

   int evtrack = 0;
   int evcluster = 0;

   int maxnclusters = 0;
   int maxntracks_reco = 0;

   int Size = tree -> GetEntries();
   cout << Size << endl;

   for(int i = 0; i < Size; i++) {
      tree -> GetEntry(i);

//      if(nclusters > maxnclusters) maxnclusters = nclusters;
//      if(ntracks_reco > maxntracks_reco) maxntracks_reco = ntracks_reco;

      evtrack = 1;
      evcluster = 1;
      
      if(nclusters > 0) {
         for(int j = 0; j < nclusters; j++) {
            h_cluster -> Fill(clucosmicscore_tagger[j]);
            if(clucosmicscore_tagger[j] == 0) evcluster = 0;
         }
         if(evcluster == 0) h_cluster_ev -> Fill(0);
         else h_cluster_ev -> Fill(1);
      }

      if(ntracks_reco > 0) {
         for(int j = 0; j < ntracks_reco; j++) {
            h_track -> Fill(trkcosmicscore_tagger[j]);
            if(trkcosmicscore_tagger[j] == 0) evtrack = 0;
         }
         if(evtrack == 0) h_track_ev -> Fill(0);
         else h_track_ev -> Fill(1);
      }
   }

   //cout << maxnclusters << "\t" << maxntracks_reco << endl;

   (h_cluster -> GetXaxis()) -> SetTitle("cosmicscore");
   (h_cluster -> GetYaxis()) -> SetTitle("# clusters / # tracks");
   h_cluster -> SetLineWidth(2);
   h_cluster -> SetLineColor(4);
   h_cluster -> Scale(1./h_cluster -> Integral());
   h_cluster -> Draw();

   h_track -> SetLineWidth(2);
   h_track -> SetLineColor(2);
   h_track -> Scale(1./h_track -> Integral());
   h_track -> Draw("same");

   TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
   c2 -> SetGridx(1);
   c2 -> SetGridy(1);
   c2 -> SetTicks(1, 1);

   (h_cluster_ev -> GetXaxis()) -> SetTitle("cosmicscore");
   (h_cluster_ev -> GetYaxis()) -> SetTitle("# events tagged");
   h_cluster_ev -> SetLineWidth(2);
   h_cluster_ev -> SetLineColor(4);
   //h_cluster_ev -> Draw();

   h_track_ev -> SetLineWidth(2);
   h_track_ev -> SetLineColor(2);
   h_track_ev -> Draw("");

   return 0;

}
