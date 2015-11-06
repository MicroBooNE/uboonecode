using namespace std;
#include <iostream>
#include <string>
#include <algorithm>
#include "TChain.h"
#include "TCanvas.h"
#include "TMath.h"

int clusterPlots_2clust_JKKJedits() {

   int   plane       = 2;   //2 = Y plane
   int   minClusters = 2;   //number of clusters
   int   maxDistance = 10;   //difference in wires
   int   maxTime     = 30;   //difference in time ticks
   float maxAngle    = .5;   //start Angle max
   int   maxLength   = 500;   
   int   maxLengthcut= 200;   
   int   minHits     = 10;   //number of hits in cluster
   int   minWire     = 5;
   int   maxWire     = 3420;
   int   minTick     = 3210;
   int   maxTick     = 6370;
   float maxScore    = 0.4; //cosmic score
   int   minLen      = 30; // minimum good cluster len in z wires
   int   minTime     = 30; // minimum good cluster len in time ticks (x)
  

   std::ofstream myFile("features-newdir8.csv");
   myFile << "type,event,minangle,openangle,charge,maxzlen,totalZ,difflen" << std::endl;

   TFile f                 = TFile("clustHists-newdir8.root","new");
   TH1F* hOpenAngle_nu     = new TH1F("Opening_Angle_Nu",";Opening Angle [rad];",30,0.,1.6);
   TH1F* hStartAngle_nu    = new TH1F("Start_Angle_Nu",";Start Angle [rad];",30,0.,2.);
   TH1F* hStartCharge_nu   = new TH1F("Start_Charge_Nu",";Starting Charge [ADC];",50,0,2048);
   TH1F* hHitIntensity_nu  = new TH1F("Hit_Intensity_Nu",";Hit Intensity;",50,0,10000);
   TH1F* hHitEvtIntensity_nu  = new TH1F("Hit_Evt_Intensity_Nu",";Hit Intensity;",50,0,30000);
   TH1F* hNVertices_nu     = new TH1F("Number_of_Vertices_Nu",";Number of Vertices;",11,0,10);
   TH1F* hTotalZ_nu        = new TH1F("Sum_of_Length_in_Z_Nu",";Total Z [wires];",100,0,2000);
   TH1F* hMaxLen_nu        = new TH1F("Max_of_Length_in_Z_Nu",";Max Length [wires];",100,0,2000);
   TH1F* hDiffLen_nu       = new TH1F("Diff_of_Length_in_Z_Nu",";Length Diff. [wires';",100,0,2000);
   TH1F* hTotalHits_nu     = new TH1F("Total_Number_of_Hits_Nu",";Number of Hits;",100,0,500);
   TH2F* hZLenVAngle_nu    = new TH2F("Z_Length_vs_Angle_Nu",";Z Length [wires];",100,0,2000,100,0.,2.);
   TH2F* hZLenVCharge_nu   = new TH2F("Z_Length_vs_Charge_Nu",";Z Length [wires];",100,0,2000,100,0.,2048);

   TH1F* hOpenAngle_cs     = new TH1F("Opening_Angle_Cosmic",";Opening Angle [rad];",30,0.,1.6);
   TH1F* hStartAngle_cs    = new TH1F("Start_Angle_Cosmic",";Start Angle [rad];",30,0.,2.);
   TH1F* hStartCharge_cs   = new TH1F("Start_Charge_Cosmic",";Starting Charge [ADC];",50,0,2048);
   TH1F* hHitIntensity_cs  = new TH1F("Hit_Intensity_Cosmic",";Hit Intensity;",50,0,10000);
   TH1F* hHitEvtIntensity_cs  = new TH1F("Hit_Evt_Intensity_Cosmic",";Hit Intensity;",50,0,30000);
   TH1F* hNVertices_cs     = new TH1F("Number_of_Vertices_Cosmic",";Number of Vertices;",11,0,10);
   TH1F* hTotalZ_cs        = new TH1F("Sum_of_Length_in_Z_Cosmic",";Total Z [wires];",100,0,2000);
   TH1F* hMaxLen_cs        = new TH1F("Max_of_Length_in_Z_Cosmic",";Max Length [wires];",100,0,2000);
   TH1F* hDiffLen_cs       = new TH1F("Diff_of_Length_in_Z_Cosmic",";Length Diff. [wires';",100,0,2000);
   TH1F* hTotalHits_cs     = new TH1F("Total_Number_of_Hits_Cosmic",";Number of Hits;",100,0,500);
   TH2F* hZLenVAngle_cs    = new TH2F("Z_Length_vs_Angle_Cosmic",";Z Length [wires];",100,0,2000,100,0.,2.);
   TH2F* hZLenVCharge_cs   = new TH2F("Z_Length_vs_Charge_Cosmic",";Z Length [wires];",100,0,2000,100,0.,2048);
   
   TChain *tree = new TChain("analysistree/anatree");
   tree -> Add("/uboone/data/users/aschu/CosmicTagging/anahist.root");
   Int_t  nentries = tree->GetEntries();

   const int kMaxClusters_nu = 500;

   Int_t   run;
   Int_t   subrun;
   Int_t   event;
   Short_t nclusters;                                 //number of clusters in a  given event
   Short_t cluncosmictags_tagger[kMaxClusters_nu];    //No. of cosmic tags associated to this cluster
   Float_t clucosmicscore_tagger[kMaxClusters_nu];    //Cosmic score associated to this cluster. In the case of more than one tag, the first one is associated.
   Short_t clucosmictype_tagger[kMaxClusters_nu];     //Cosmic tag type for this cluster.
   Short_t cluster_StartWire[kMaxClusters_nu];        //wire coordinate of the start of the cluster 
   Short_t cluster_StartTick[kMaxClusters_nu];        //tick coordinate of the start of the cluster in time ticks
   Float_t cluster_StartAngle[kMaxClusters_nu];       //angle of the start of the cluster in radians
   Float_t cluster_StartCharge[kMaxClusters_nu];      //charge on the first wire of the cluster in ADC
   Short_t cluster_EndWire[kMaxClusters_nu];          //wire coordinate of the end of the cluster
   Short_t cluster_EndTick[kMaxClusters_nu];          //tick coordinate of the end of the cluster in time ticks
   Float_t cluster_EndAngle[kMaxClusters_nu];         //angle of the end of the cluster in radians
   Float_t cluster_EndCharge[kMaxClusters_nu];
   Short_t cluster_NHits[kMaxClusters_nu];            //Number of hits in the cluster
   Short_t clusterView[kMaxClusters_nu];              //which plane this cluster belongs to        
   Float_t cluster_Integral[kMaxClusters_nu];            //returns the total charge of the cluster from hit shape in ADC
  

 
   std::cout << "Running over neutrino Set. " << nentries << " entries." << std::endl;
   //for(Int_t ientry=0; ientry < nentries; ientry++)
   for(Int_t ientry=0; ientry < 5000; ientry++)
   {
      tree -> SetBranchAddress("run", &run);
      tree -> SetBranchAddress("subrun", &subrun);
      tree -> SetBranchAddress("event", &event);
      tree -> SetBranchAddress("cluncosmictags_tagger", cluncosmictags_tagger);
      tree -> SetBranchAddress("clucosmicscore_tagger", clucosmicscore_tagger);
      tree -> SetBranchAddress("clucosmictype_tagger", clucosmictype_tagger);
      tree -> SetBranchAddress("nclusters", &nclusters);
      tree -> SetBranchAddress("cluster_StartWire", cluster_StartWire);
      tree -> SetBranchAddress("cluster_EndWire", cluster_EndWire);
      tree -> SetBranchAddress("cluster_StartTick", cluster_StartTick);
      tree -> SetBranchAddress("cluster_EndTick", cluster_EndTick);
      tree -> SetBranchAddress("cluster_StartAngle", cluster_StartAngle);
      tree -> SetBranchAddress("cluster_EndAngle", cluster_EndAngle);
      tree -> SetBranchAddress("cluster_StartCharge", cluster_StartCharge);
      tree -> SetBranchAddress("cluster_EndCharge", cluster_EndCharge);
      tree -> SetBranchAddress("cluster_NHits", cluster_NHits);
      tree -> SetBranchAddress("clusterView", clusterView);
      tree -> SetBranchAddress("cluster_Integral", cluster_Integral);

      tree -> GetEntry(ientry);
      if(ientry % 100 == 0) std::cout << "this is event\t" << ientry << endl;
      if(nclusters < minClusters) continue; 

      // Check that there are at least 2 good clusters
      // Add index to vector
      std::vector<int> goodClusts;
      for(int j = 0; j < nclusters; j++)
      {
        // std::cout<<"Cluster Start Angle: "<<cluster_StartAngle[j]<<std::endl;
	if(clusterView[j] != plane) continue;
        if(cluster_NHits[j] < minHits) continue;
        if(clucosmicscore_tagger[j] >= maxScore) continue;
        if((cluster_StartTick[j] < minTick) || (cluster_StartTick[j] > maxTick) || (cluster_EndTick[j] < minTick) || (cluster_EndTick[j] > maxTick)) continue;
        if((cluster_StartWire[j] < minWire) || (cluster_StartWire[j] > maxWire) || (cluster_EndWire[j] < minWire) || (cluster_EndWire[j] > maxWire)) continue;
        if(fabs(cluster_StartAngle[j]) > maxAngle && fabs(cluster_StartWire[j]-cluster_EndWire[j]) > maxLengthcut) continue;
        if(fabs(cluster_StartWire[j] - cluster_EndWire[j]) < minLen && fabs(cluster_StartTick[j] - cluster_EndTick[j]) < minTime) continue;
        goodClusts.push_back(j);
      }
      if(goodClusts.size() < minClusters) continue;

      Int_t numVtcs = 0;
      for(unsigned int k1 = 0; k1 < goodClusts.size(); ++k1)
      {
        int cid1 = goodClusts.at(k1);
        if(fabs(cluster_EndWire[cid1] - cluster_StartWire[cid1]) > maxLength && cluster_StartWire[cid1] < cluster_EndWire[cid1])
        {
          // Starts with StartWire
          for(unsigned int k2 = 0; k2 < goodClusts.size(); ++k2)
          {
            if(k2 == k1) continue;
            int cid2 = goodClusts.at(k2);
            if(fabs(cluster_StartWire[cid2] - cluster_EndWire[cid2]) > fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1])) continue;
            if((fabs(cluster_StartWire[cid1] - cluster_StartWire[cid2]) <= maxDistance) && (fabs(cluster_StartTick[cid1] - cluster_StartTick[cid2]) <= maxTime))
            {
              if(cluster_StartCharge[cid2]>cluster_StartCharge[cid1]) 
              {
 		Float_t openangle=fabs(cluster_StartAngle[cid1] - cluster_StartAngle[cid2]);
            	if(openangle>0.1 && openangle<1.35)
		{
		  numVtcs++;
	      	  Int_t   maxzlen  = fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]);
                  Float_t difflen  = fabs(fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]) - fabs(cluster_StartWire[cid2] - cluster_EndWire[cid2]));
              	  Float_t minangle = std::min(fabs(cluster_StartAngle[cid1]),fabs(cluster_StartAngle[cid2]));
              	  Float_t charge   = cluster_StartCharge[cid1] + cluster_StartCharge[cid2];
		  Float_t totalZ   = fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]) + fabs(cluster_StartWire[cid2] - cluster_EndWire[cid2]);
              	  Float_t HitIntensity= cluster_StartCharge[cid2];
              	  Float_t EvtHitIntensity= cluster_Integral[cid2]+cluster_Integral[cid2];
              	  hOpenAngle_nu   -> Fill(openangle);
              	  hStartAngle_nu  -> Fill(minangle);
              	  hStartCharge_nu -> Fill(charge);
              	  hHitIntensity_nu-> Fill(HitIntensity);
              	  hHitEvtIntensity_nu-> Fill(EvtHitIntensity);
              	  hTotalZ_nu      -> Fill(totalZ);
              	  hMaxLen_nu      -> Fill(maxzlen);
                  hDiffLen_nu     -> Fill(difflen);
              	  hTotalHits_nu   -> Fill(cluster_NHits[cid1] + cluster_NHits[cid2]);
              	  hZLenVAngle_nu  -> Fill(maxzlen,minangle);
              	  hZLenVCharge_nu -> Fill(maxzlen,charge);
                  myFile << "nu," << ientry << "," << minangle << "," << openangle << "," << charge << "," << maxzlen << "," << totalZ << "," <<difflen <<","<<EvtHitIntensity<<","<<HitIntensity<<std::endl;
                }
              }
            }
            else if((fabs(cluster_StartWire[cid1] - cluster_EndWire[cid2]) <= maxDistance) && (fabs(cluster_StartTick[cid1] - cluster_EndTick[cid2]) <= maxTime))
            {
              if(cluster_EndCharge[cid2]>cluster_StartCharge[cid1]) 
              {
 		Float_t openangle=fabs(cluster_StartAngle[cid1] - cluster_EndAngle[cid2]);
            	if(openangle>0.1 && openangle<1.35)
		{
		  numVtcs++;
	      	  Int_t   maxzlen  = fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]);
                  Float_t difflen  = fabs(fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]) - fabs(cluster_StartWire[cid2] - cluster_EndWire[cid2]));
              	  Float_t minangle = std::min(fabs(cluster_StartAngle[cid1]),fabs(cluster_EndAngle[cid2]));
              	  Float_t charge   = cluster_StartCharge[cid1] + cluster_EndCharge[cid2];
		  Float_t totalZ   = fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]) + fabs(cluster_StartWire[cid2] - cluster_EndWire[cid2]);
              	  Float_t HitIntensity= cluster_EndCharge[cid2];
              	  Float_t EvtHitIntensity= cluster_Integral[cid2]+cluster_Integral[cid2];
              	  hOpenAngle_nu   -> Fill(openangle);
              	  hStartAngle_nu  -> Fill(minangle);
              	  hStartCharge_nu -> Fill(charge);
              	  hHitIntensity_nu-> Fill(HitIntensity);
              	  hHitEvtIntensity_nu-> Fill(EvtHitIntensity);
              	  hTotalZ_nu      -> Fill(totalZ);
              	  hMaxLen_nu      -> Fill(maxzlen);
                  hDiffLen_nu     -> Fill(difflen);
              	  hTotalHits_nu   -> Fill(cluster_NHits[cid1] + cluster_NHits[cid2]);
              	  hZLenVAngle_nu  -> Fill(maxzlen,minangle);
              	  hZLenVCharge_nu -> Fill(maxzlen,charge);
                  myFile << "nu," << ientry << "," << minangle << "," << openangle << "," << charge << "," << maxzlen << "," << totalZ << "," <<difflen <<","<<EvtHitIntensity<<","<<HitIntensity<<std::endl;
                }
              }
            }
          }
        } //end startwire
        else if(fabs(cluster_EndWire[cid1] - cluster_StartWire[cid1]) > maxLength && cluster_StartWire[cid1] > cluster_EndWire[cid1])
        {
          // Starts with EndWire
          for(unsigned int k2 = 0; k2 < goodClusts.size(); ++k2)
          {
            if(k2 == k1) continue;
            int cid2 = goodClusts.at(k2);
            if(fabs(cluster_StartWire[cid2] - cluster_EndWire[cid2]) > fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1])) continue;
            if((fabs(cluster_EndWire[cid1] - cluster_StartWire[cid2]) <= maxDistance) && (fabs(cluster_EndTick[cid1] - cluster_StartTick[cid2]) <= maxTime))
            {
              if(cluster_StartCharge[cid2]>cluster_EndCharge[cid1]) 
              {
 		Float_t openangle=fabs(cluster_EndAngle[cid1] - cluster_StartAngle[cid2]);
            	if(openangle>0.1 && openangle<1.35)
		{
		  numVtcs++;
	      	  Int_t   maxzlen  = fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]);
                  Float_t difflen  = fabs(fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]) - fabs(cluster_StartWire[cid2] - cluster_EndWire[cid2]));
              	  Float_t minangle = std::min(fabs(cluster_EndAngle[cid1]),fabs(cluster_StartAngle[cid2]));
              	  Float_t charge   = cluster_EndCharge[cid1] + cluster_StartCharge[cid2];
		  Float_t totalZ   = fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]) + fabs(cluster_StartWire[cid2] - cluster_EndWire[cid2]);
		  Float_t HitIntensity  =  cluster_StartCharge[cid2];
		  Float_t EvtHitIntensity  =  cluster_Integral[cid2]+cluster_Integral[cid2];
             	  hOpenAngle_nu   -> Fill(openangle);
              	  hStartAngle_nu  -> Fill(minangle);
              	  hStartCharge_nu -> Fill(charge);
              	  hHitIntensity_nu-> Fill(HitIntensity);
              	  hHitEvtIntensity_nu-> Fill(EvtHitIntensity);
              	  hTotalZ_nu      -> Fill(totalZ);
              	  hMaxLen_nu      -> Fill(maxzlen);
                  hDiffLen_nu     -> Fill(difflen);
              	  hTotalHits_nu   -> Fill(cluster_NHits[cid1] + cluster_NHits[cid2]);
              	  hZLenVAngle_nu  -> Fill(maxzlen,minangle);
              	  hZLenVCharge_nu -> Fill(maxzlen,charge);
                  myFile << "nu," << ientry << "," << minangle << "," << openangle << "," << charge << "," << maxzlen << "," << totalZ << "," <<difflen <<","<<EvtHitIntensity<<","<<HitIntensity<<std::endl;
                }
              }
            }
            else if((fabs(cluster_EndWire[cid1] - cluster_EndWire[cid2]) <= maxDistance) && (fabs(cluster_EndTick[cid1] - cluster_EndTick[cid2]) <= maxTime))
            {
              if(cluster_EndCharge[cid2]>cluster_EndCharge[cid1]) 
              {
 		Float_t openangle=fabs(cluster_EndAngle[cid1] - cluster_EndAngle[cid2]);
            	if(openangle>0.1 && openangle<1.35)
		{
		  numVtcs++;
	      	  Int_t   maxzlen  = fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]);
                  Float_t difflen  = fabs(fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]) - fabs(cluster_StartWire[cid2] - cluster_EndWire[cid2]));
              	  Float_t minangle = std::min(fabs(cluster_EndAngle[cid1]),fabs(cluster_EndAngle[cid2]));
              	  Float_t charge   = cluster_EndCharge[cid1] + cluster_EndCharge[cid2];
		  Float_t totalZ   = fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]) + fabs(cluster_StartWire[cid2] - cluster_EndWire[cid2]);
              	  Float_t HitIntensity= cluster_EndCharge[cid2];
              	  Float_t EvtHitIntensity= cluster_Integral[cid2]+cluster_Integral[cid2];
	          hOpenAngle_nu   -> Fill(openangle);
              	  hStartAngle_nu  -> Fill(minangle);
              	  hStartCharge_nu -> Fill(charge);
              	  hHitIntensity_nu-> Fill(HitIntensity);
              	  hHitEvtIntensity_nu-> Fill(EvtHitIntensity);
              	  hTotalZ_nu      -> Fill(totalZ);
              	  hMaxLen_nu      -> Fill(maxzlen);
                  hDiffLen_nu     -> Fill(difflen);
              	  hTotalHits_nu   -> Fill(cluster_NHits[cid1] + cluster_NHits[cid2]);
              	  hZLenVAngle_nu  -> Fill(maxzlen,minangle);
              	  hZLenVCharge_nu -> Fill(maxzlen,charge);
                  myFile << "nu," << ientry << "," << minangle << "," << openangle << "," << charge << "," << maxzlen << "," << totalZ <<","<<difflen << ","<<EvtHitIntensity<<","<<HitIntensity<<std::endl;
                }
              }
            }
          }
        } //end endwire
      } // end loop over clusters
      hNVertices_nu->Fill(numVtcs);
   } // end loop over events
   std::cout << "Num. evt w/ good clusters (1 day): " << hNVertices_nu->GetEntries()*.124 << std::endl;
   std::cout << "Num. evt w/ matched clusters(1 day): " <<(hNVertices_nu->GetEntries()-hNVertices_nu->GetBinContent(1))*.124<< std::endl;
   std::cout << "Num. matching clusters (1 day): " << hOpenAngle_nu->GetEntries()*.124 << std::endl;

   TChain *tree = new TChain("analysistree/anatree");
   tree -> Add("/uboone/data/users/aschu/CosmicTagging/anahist-cosmic.root");
   Int_t  nentries = tree->GetEntries();

   const int kMaxClusters_cs = 1500;

   Int_t   run;
   Int_t   subrun;
   Int_t   event;
   Short_t nclusters_cs;                                    //number of clusters in a  given event
   Short_t cluncosmictags_tagger_cs[kMaxClusters_cs];    //No. of cosmic tags associated to this cluster
   Float_t clucosmicscore_tagger_cs[kMaxClusters_cs];    //Cosmic score associated to this cluster. In the case of more than one tag, the first one is associated.
   Short_t clucosmictype_tagger_cs[kMaxClusters_cs];     //Cosmic tag type for this cluster.
   Short_t cluster_StartWire_cs[kMaxClusters_cs];        //wire coordinate of the start of the cluster 
   Short_t cluster_StartTick_cs[kMaxClusters_cs];        //tick coordinate of the start of the cluster in time ticks
   Float_t cluster_StartAngle_cs[kMaxClusters_cs];       //angle of the start of the cluster in radians
   Float_t cluster_StartCharge_cs[kMaxClusters_cs];      //charge on the first wire of the cluster in ADC
   Short_t cluster_EndWire_cs[kMaxClusters_cs];          //wire coordinate of the end of the cluster
   Short_t cluster_EndTick_cs[kMaxClusters_cs];          //tick coordinate of the end of the cluster in time ticks
   Float_t cluster_EndAngle_cs[kMaxClusters_cs];         //angle of the end of the cluster in radians
   Float_t cluster_EndCharge_cs[kMaxClusters_cs];      //charge on the first wire of the cluster in ADC
   Short_t cluster_NHits_cs[kMaxClusters_cs];            //Number of hits in the cluster
   Short_t clusterView_cs[kMaxClusters_cs];              //which plane this cluster belongs to        
   Float_t cluster_Integral_cs[kMaxClusters_cs];            //returns the total charge of the cluster from hit shape in ADC

   std::cout << "Running over cosmic Set. " << nentries << " entries." << std::endl;
   //for(Int_t ientry=0; ientry < nentries; ientry++)
   for(Int_t ientry=0; ientry < 5000; ientry++)
   {
      tree -> SetBranchAddress("run", &run);
      tree -> SetBranchAddress("subrun", &subrun);
      tree -> SetBranchAddress("event", &event);
      tree -> SetBranchAddress("nclusters", &nclusters_cs);
      tree -> SetBranchAddress("cluncosmictags_tagger", cluncosmictags_tagger_cs);
      tree -> SetBranchAddress("clucosmicscore_tagger", clucosmicscore_tagger_cs);
      tree -> SetBranchAddress("clucosmictype_tagger", clucosmictype_tagger_cs);
      tree -> SetBranchAddress("cluster_StartWire", cluster_StartWire_cs);
      tree -> SetBranchAddress("cluster_EndWire", cluster_EndWire_cs);
      tree -> SetBranchAddress("cluster_StartTick", cluster_StartTick_cs);
      tree -> SetBranchAddress("cluster_EndTick", cluster_EndTick_cs);
      tree -> SetBranchAddress("cluster_StartAngle", cluster_StartAngle_cs);
      tree -> SetBranchAddress("cluster_EndAngle", cluster_EndAngle_cs);
      tree -> SetBranchAddress("cluster_StartCharge", cluster_StartCharge_cs);
      tree -> SetBranchAddress("cluster_EndCharge", cluster_EndCharge_cs);
      tree -> SetBranchAddress("cluster_NHits", cluster_NHits_cs);
      tree -> SetBranchAddress("clusterView", clusterView_cs);
      tree -> SetBranchAddress("cluster_Integral", cluster_Integral_cs);

      tree -> GetEntry(ientry);
      if(ientry % 100 == 0) std::cout << "this is event\t" << ientry << endl;
      if(nclusters_cs < minClusters) continue;

      // Check that there are at least 2 good clusters
      // Add index to vector
      std::vector<int> goodClusts;
      for(int j = 0; j < nclusters_cs; j++)
      {
         if(clusterView_cs[j] != plane) continue;
         if(cluster_NHits_cs[j] < minHits) continue;
         if(clucosmicscore_tagger_cs[j] >= maxScore) continue;
         if((cluster_StartTick_cs[j] < minTick) || (cluster_StartTick_cs[j] > maxTick) || (cluster_EndTick_cs[j] < minTick) || (cluster_EndTick_cs[j] > maxTick)) continue;
         if((cluster_StartWire_cs[j] < minWire) || (cluster_StartWire_cs[j] > maxWire) || (cluster_EndWire_cs[j] < minWire) || (cluster_EndWire_cs[j] > maxWire)) continue;
         if(fabs(cluster_StartAngle_cs[j]) > maxAngle && fabs(cluster_StartWire_cs[j]-cluster_EndWire_cs[j]) > maxLengthcut) continue;
         if(fabs(cluster_StartWire_cs[j] - cluster_EndWire_cs[j]) < minLen && fabs(cluster_StartTick_cs[j] - cluster_EndTick_cs[j]) < minTime) continue;
         goodClusts.push_back(j);
      }
      if(goodClusts.size() < minClusters) continue;

      Int_t numVtcs = 0;
      for(unsigned int k1 = 0; k1 < goodClusts.size(); ++k1)
      {
        int cid1 = goodClusts.at(k1);
        if(fabs(cluster_EndWire_cs[cid1] - cluster_StartWire_cs[cid1]) > maxLength && cluster_StartWire_cs[cid1] < cluster_EndWire_cs[cid1])
        {
          // Starts with StartWire
          for(unsigned int k2 = 0; k2 < goodClusts.size(); ++k2)
          {
            if(k2 == k1) continue;
            int cid2 = goodClusts.at(k2);
            if(fabs(cluster_StartWire_cs[cid2] - cluster_EndWire_cs[cid2]) > fabs(cluster_StartWire_cs[cid1] - cluster_EndWire_cs[cid1])) continue;
            if((fabs(cluster_StartWire_cs[cid1] - cluster_StartWire_cs[cid2]) <= maxDistance) && (fabs(cluster_StartTick_cs[cid1] - cluster_StartTick_cs[cid2]) <= maxTime))
            {
              if(cluster_StartCharge_cs[cid2]>cluster_StartCharge_cs[cid1]) 
              {
 		Float_t openangle=fabs(cluster_StartAngle_cs[cid1] - cluster_StartAngle_cs[cid2]);
                if(openangle>0.1 && openangle<1.35)
		{
		  numVtcs++;
	      	  Int_t   maxzlen  = fabs(cluster_StartWire_cs[cid1] - cluster_EndWire_cs[cid1]);
                  Float_t difflen  = fabs(fabs(cluster_StartWire_cs[cid1] - cluster_EndWire_cs[cid1]) - fabs(cluster_StartWire_cs[cid2] - cluster_EndWire_cs[cid2]));
              	  Float_t minangle = std::min(fabs(cluster_StartAngle_cs[cid1]),fabs(cluster_StartAngle_cs[cid2]));
              	  Float_t charge   = cluster_StartCharge_cs[cid1] + cluster_StartCharge_cs[cid2];
		  Float_t totalZ   = fabs(cluster_StartWire_cs[cid1] - cluster_EndWire_cs[cid1]) + fabs(cluster_StartWire_cs[cid2] - cluster_EndWire_cs[cid2]);
		  Float_t HitIntensity   = cluster_StartCharge_cs[cid2];
		  Float_t EvtHitIntensity   = cluster_Integral_cs[cid2]+cluster_Integral_cs[cid1];
              	  hOpenAngle_cs   -> Fill(openangle);
              	  hStartAngle_cs  -> Fill(minangle);
              	  hStartCharge_cs -> Fill(charge);
              	  hHitIntensity_cs-> Fill(HitIntensity);
              	  hHitEvtIntensity_cs-> Fill(EvtHitIntensity);
              	  hTotalZ_cs      -> Fill(totalZ);
              	  hMaxLen_cs      -> Fill(maxzlen);
                  hDiffLen_cs     -> Fill(difflen);
              	  hTotalHits_cs   -> Fill(cluster_NHits_cs[cid1] + cluster_NHits_cs[cid2]);
              	  hZLenVAngle_cs  -> Fill(maxzlen,minangle);
              	  hZLenVCharge_cs -> Fill(maxzlen,charge);
                  myFile<< "cos," << ientry << "," << minangle << "," << openangle << "," << charge << "," << maxzlen << "," << totalZ << "," <<difflen <<","<<EvtHitIntensity<<","<<HitIntensity<<std::endl;
                }
              }
            }
            else if((fabs(cluster_StartWire_cs[cid1] - cluster_EndWire_cs[cid2]) <= maxDistance) && (fabs(cluster_StartTick_cs[cid1] - cluster_EndTick_cs[cid2]) <= maxTime))
            {
              if(cluster_EndCharge_cs[cid2]>cluster_StartCharge_cs[cid1]) 
              {
 		Float_t openangle=fabs(cluster_StartAngle_cs[cid1] - cluster_EndAngle_cs[cid2]);
            	if(openangle>0.1 && openangle<1.35)
		{
		  numVtcs++;
	      	  Int_t   maxzlen  = fabs(cluster_StartWire_cs[cid1] - cluster_EndWire_cs[cid1]);
                  Float_t difflen  = fabs(fabs(cluster_StartWire_cs[cid1] - cluster_EndWire_cs[cid1]) - fabs(cluster_StartWire_cs[cid2] - cluster_EndWire_cs[cid2]));
              	  Float_t minangle = std::min(fabs(cluster_StartAngle_cs[cid1]),fabs(cluster_EndAngle_cs[cid2]));
              	  Float_t charge   = cluster_StartCharge_cs[cid1] + cluster_EndCharge_cs[cid2];
		  Float_t totalZ   = fabs(cluster_StartWire_cs[cid1] - cluster_EndWire_cs[cid1]) + fabs(cluster_StartWire_cs[cid2] - cluster_EndWire_cs[cid2]);
		  Float_t HitIntensity   = cluster_EndCharge_cs[cid2];
		  Float_t EvtHitIntensity   = cluster_Integral_cs[cid2]+cluster_Integral_cs[cid1];
              	  hOpenAngle_cs   -> Fill(openangle);
              	  hStartAngle_cs  -> Fill(minangle);
              	  hStartCharge_cs -> Fill(charge);
              	  hHitIntensity_cs-> Fill(HitIntensity);
              	  hHitEvtIntensity_cs-> Fill(EvtHitIntensity);
              	  hTotalZ_cs      -> Fill(totalZ);
              	  hMaxLen_cs      -> Fill(maxzlen);
                  hDiffLen_cs     -> Fill(difflen);
              	  hTotalHits_cs   -> Fill(cluster_NHits_cs[cid1] + cluster_NHits_cs[cid2]);
              	  hZLenVAngle_cs  -> Fill(maxzlen,minangle);
              	  hZLenVCharge_cs -> Fill(maxzlen,charge);
                  myFile<< "cos," << ientry << "," << minangle << "," << openangle << "," << charge << "," << maxzlen << "," << totalZ << "," <<difflen <<","<<EvtHitIntensity<<","<<HitIntensity<<std::endl;
                }
              }
            }
          }
        } //end startwire
        else if(fabs(cluster_EndWire_cs[cid1] - cluster_StartWire_cs[cid1]) > maxLength && cluster_StartWire_cs[cid1] > cluster_EndWire_cs[cid1])
        {
          // Starts with EndWire
          for(unsigned int k2 = 0; k2 < goodClusts.size(); ++k2)
          {
            if(k2 == k1) continue;
            int cid2 = goodClusts.at(k2);
            if(fabs(cluster_StartWire_cs[cid2] - cluster_EndWire_cs[cid2]) > fabs(cluster_StartWire_cs[cid1] - cluster_EndWire_cs[cid1])) continue;
            if((fabs(cluster_EndWire_cs[cid1] - cluster_StartWire_cs[cid2]) <= maxDistance) && (fabs(cluster_EndTick_cs[cid1] - cluster_StartTick_cs[cid2]) <= maxTime))
            {
              if(cluster_StartCharge_cs[cid2]>cluster_EndCharge_cs[cid1]) 
              {
 	        Float_t openangle=fabs(cluster_EndAngle_cs[cid1] - cluster_StartAngle_cs[cid2]);
            	if(openangle>0.1 && openangle<1.35)
		{
		  numVtcs++;
	      	  Int_t   maxzlen  = fabs(cluster_StartWire_cs[cid1] - cluster_EndWire_cs[cid1]);
                  Float_t difflen  = fabs(fabs(cluster_StartWire_cs[cid1] - cluster_EndWire_cs[cid1]) - fabs(cluster_StartWire_cs[cid2] - cluster_EndWire_cs[cid2]));
              	  Float_t minangle = std::min(fabs(cluster_EndAngle_cs[cid1]),fabs(cluster_StartAngle_cs[cid2]));
              	  Float_t charge   = cluster_EndCharge_cs[cid1] + cluster_StartCharge_cs[cid2];
		  Float_t totalZ   = fabs(cluster_StartWire_cs[cid1] - cluster_EndWire_cs[cid1]) + fabs(cluster_StartWire_cs[cid2] - cluster_EndWire_cs[cid2]);
		  Float_t HitIntensity   = cluster_StartCharge_cs[cid2];
		  Float_t EvtHitIntensity   = cluster_Integral_cs[cid2]+cluster_Integral_cs[cid1];
              	  hOpenAngle_cs   -> Fill(openangle);
              	  hStartAngle_cs  -> Fill(minangle);
              	  hStartCharge_cs -> Fill(charge);
              	  hHitIntensity_cs-> Fill(HitIntensity);
              	  hHitEvtIntensity_cs-> Fill(EvtHitIntensity);
              	  hTotalZ_cs      -> Fill(totalZ);
              	  hMaxLen_cs      -> Fill(maxzlen);
                  hDiffLen_cs     -> Fill(difflen);
              	  hTotalHits_cs   -> Fill(cluster_NHits_cs[cid1] + cluster_NHits_cs[cid2]);
              	  hZLenVAngle_cs  -> Fill(maxzlen,minangle);
              	  hZLenVCharge_cs -> Fill(maxzlen,charge);
                  myFile<< "cos," << ientry << "," << minangle << "," << openangle << "," << charge << "," << maxzlen << "," << totalZ << "," <<difflen <<","<<EvtHitIntensity<<","<<HitIntensity<<std::endl;
                }
              }
            }
            else if((fabs(cluster_EndWire_cs[cid1] - cluster_EndWire_cs[cid2]) <= maxDistance) && (fabs(cluster_EndTick_cs[cid1] - cluster_EndTick_cs[cid2]) <= maxTime))
            {
              if(cluster_EndCharge_cs[cid2]>cluster_EndCharge_cs[cid1]) 
              {
 		Float_t openangle=fabs(cluster_EndAngle_cs[cid1] - cluster_EndAngle_cs[cid2]);
                if(openangle>0.1 && openangle<1.35)
		{
		  numVtcs++;
	      	  Int_t   maxzlen  = fabs(cluster_StartWire_cs[cid1] - cluster_EndWire_cs[cid1]);
                  Float_t difflen  = fabs(fabs(cluster_StartWire_cs[cid1] - cluster_EndWire_cs[cid1]) - fabs(cluster_StartWire_cs[cid2] - cluster_EndWire_cs[cid2]));
              	  Float_t minangle = std::min(fabs(cluster_EndAngle_cs[cid1]),fabs(cluster_EndAngle_cs[cid2]));
              	  Float_t charge   = cluster_EndCharge_cs[cid1] + cluster_EndCharge_cs[cid2];
		  Float_t totalZ   = fabs(cluster_StartWire_cs[cid1] - cluster_EndWire_cs[cid1]) + fabs(cluster_StartWire_cs[cid2] - cluster_EndWire_cs[cid2]);
		  Float_t HitIntensity   = cluster_EndCharge_cs[cid2];
		  Float_t EvtHitIntensity   = cluster_Integral_cs[cid2]+cluster_Integral_cs[cid1];
              	  hOpenAngle_cs   -> Fill(openangle);
              	  hStartAngle_cs  -> Fill(minangle);
              	  hStartCharge_cs -> Fill(charge);
              	  hHitIntensity_cs-> Fill(HitIntensity);
              	  hHitEvtIntensity_cs-> Fill(EvtHitIntensity);
              	  hTotalZ_cs      -> Fill(totalZ);
              	  hMaxLen_cs      -> Fill(maxzlen);
                  hDiffLen_cs     -> Fill(difflen);
              	  hTotalHits_cs   -> Fill(cluster_NHits_cs[cid1] + cluster_NHits_cs[cid2]);
              	  hZLenVAngle_cs  -> Fill(maxzlen,minangle);
              	  hZLenVCharge_cs -> Fill(maxzlen,charge);
                  myFile<< "cos," << ientry << "," << minangle << "," << openangle << "," << charge << "," << maxzlen << "," << totalZ << "," <<difflen <<","<<EvtHitIntensity<<","<<HitIntensity<<std::endl;
                }
              }
            }
          }
        } //end endwire
      } // end loop over clusters
      hNVertices_cs->Fill(numVtcs);
   } // end event loop

   myFile.close();
   //std::cout << "Num. evt w/ good clusters (1 day): " << hNVertices_cs->GetEntries()*18400./9200.*.031*540. << std::endl;
   //std::cout << "Num. matching clusters (1 day): " << hOpenAngle_cs->GetEntries()*18400./9200.*.031*540. << std::endl;
   std::cout << "Num. evt w/ good clusters (1 day): " << hNVertices_cs->GetEntries()*.124*540. << std::endl;
   std::cout << "Num. evt w/ matched clusters (1 day): " << (hNVertices_cs->GetEntries()-hNVertices_cs->GetBinContent(1))*.124*540. << std::endl;
   std::cout << "Num. matching clusters (1 day): " << hOpenAngle_cs->GetEntries()*.124*540. << std::endl;

   Float_t norm_nu = .124;
   //Float_t norm_nu = .031;
   hNVertices_nu     -> Scale(norm_nu);
   hOpenAngle_nu     -> Scale(norm_nu);
   hStartAngle_nu    -> Scale(norm_nu);
   hStartCharge_nu   -> Scale(norm_nu);
   hTotalZ_nu        -> Scale(norm_nu);
   hMaxLen_nu        -> Scale(norm_nu);
   hDiffLen_nu       -> Scale(norm_nu);
   hHitIntensity_nu  -> Scale(norm_nu);
   hHitEvtIntensity_nu  -> Scale(norm_nu);
   hTotalHits_nu     -> Scale(norm_nu);
   hZLenVAngle_nu    -> Scale(norm_nu);
   hZLenVCharge_nu   -> Scale(norm_nu);

   //Float_t norm_cs = 18400./9200. * 0.031 * 540;
   Float_t norm_cs =  .124* 540;
  // Float_t norm_cs =  18400./9200.*.031* 540;
   hNVertices_cs     -> Scale(norm_cs);
   hOpenAngle_cs     -> Scale(norm_cs);
   hStartAngle_cs    -> Scale(norm_cs);
   hStartCharge_cs   -> Scale(norm_cs);
   hTotalZ_cs        -> Scale(norm_cs);
   hHitIntensity_cs  -> Scale(norm_cs);
   hHitEvtIntensity_cs  -> Scale(norm_cs);
   hMaxLen_cs        -> Scale(norm_cs);
   hDiffLen_cs       -> Scale(norm_cs);
   hTotalHits_cs     -> Scale(norm_cs);
   hZLenVAngle_cs    -> Scale(norm_cs);
   hZLenVCharge_cs   -> Scale(norm_cs);

   hNVertices_nu     -> Write();
   hOpenAngle_nu     -> Write();
   hStartAngle_nu    -> Write();
   hStartCharge_nu   -> Write();
   hTotalZ_nu        -> Write();
   hMaxLen_nu        -> Write();
   hDiffLen_nu       -> Write();
   hTotalHits_nu     -> Write();
   hHitIntensity_nu  -> Write();
   hHitEvtIntensity_nu  -> Write();
   hZLenVAngle_nu    -> Write();
   hZLenVCharge_nu   -> Write();
   
   hNVertices_cs     -> Write();
   hOpenAngle_cs     -> Write();
   hStartAngle_cs    -> Write();
   hStartCharge_cs   -> Write();
   hTotalZ_cs        -> Write();
   hMaxLen_cs        -> Write();
   hDiffLen_cs       -> Write();
   hTotalHits_cs     -> Write();
   hHitIntensity_cs  -> Write();
   hHitEvtIntensity_cs  -> Write();
   hZLenVAngle_cs    -> Write();
   hZLenVCharge_cs   -> Write();

   return 0;

}
