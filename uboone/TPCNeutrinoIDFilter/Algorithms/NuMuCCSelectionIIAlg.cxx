/**
 *  @file   NuMuCCSelectionIIAlg.cxx
 * 
 *  @brief  Implementation of the Selection II Neutrino ID alg
 *          This module outputs associations between vertices
 *          and tracks that are found to be within the cut value
 *
 *  @authors xiao.luo@yale.edu, tjyang@fnal.gov
 *
 */

// The main include
#include "uboone/TPCNeutrinoIDFilter/Algorithms/NuMuCCSelectionIIAlg.h"

// Framework Includes
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

//------------------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace neutrinoid {

NuMuCCSelectionIIAlg::NuMuCCSelectionIIAlg(fhicl::ParameterSet const &pset) :
    fMyProducerModule(0),
    fGeometry(lar::providerFrom<geo::Geometry>()),
    fDetector(lar::providerFrom<detinfo::DetectorPropertiesService>())
{
    this->reconfigure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

NuMuCCSelectionIIAlg::~NuMuCCSelectionIIAlg()
{
}
    
void NuMuCCSelectionIIAlg::reconfigure(fhicl::ParameterSet const &inputPset)
{
    // Assume we could be called externally with the top level module's complete parameter set
    const fhicl::ParameterSet& pset = inputPset.get<fhicl::ParameterSet>("NuMuCCSelectionIIAlg");
    
    fTrackModuleLabel        = pset.get<std::string> ("TrackModuleLabel");
    fVertexModuleLabel       = pset.get<std::string> ("VertexModuleLabel");
    fOpFlashModuleLabel      = pset.get<std::string> ("OpFlashModuleLabel");
    fCalorimetryModuleLabel  = pset.get<std::string> ("CalorimetryModuleLabel");
    
    fDistToEdgeX             = fGeometry->DetHalfWidth()   - pset.get<double>("DistToEdgeX",   20.);
    fDistToEdgeY             = fGeometry->DetHalfHeight()  - pset.get<double>("DistToEdgeY",   20.);
    fDistToEdgeZ             = fGeometry->DetLength() / 2. - pset.get<double>("DistToEdgeZ",   10.);
    
    fBeamMin                 = pset.get<double>      ("BeamMin",                              3.3);
    fBeamMax                 = pset.get<double>      ("BeamMax",                              4.9);
    fPEThresh                = pset.get<double>      ("PEThresh",                              50.);
    fTrk2FlashDist           = pset.get<double>      ("Trk2FlashDist",                         70.);
    fMinTrk2VtxDist          = pset.get<double>      ("MinTrk2VtxDist",                         3.);
    fMinTrackLen             = pset.get<double>      ("MinTrackLen",                           15.);
    fMaxCosineAngle          = pset.get<double>      ("MaxCosineAngle",                        0.9);
    fMaxCosy1stTrk           = pset.get<double>      ("MaxCosy1stTrk",                         0.6);
    fMinTrackLen2ndTrk       = pset.get<double>      ("MinTrackLen2ndTrk",                     30.);
    fMaxCosySingle           = pset.get<double>      ("MaxCosySingle",                         0.7);
    fMinTrackLenSingle       = pset.get<double>      ("MinTrackLenSingle",                     40.);
    fMindEdxRatioSingle      = pset.get<double>      ("MindEdxRatioSingle",                    1.5);
    fMaxTrkLengthySingle     = pset.get<double>      ("MaxTrkLengthySingle",                   25.);
    fMinStartdEdx1stTrk      = pset.get<double>      ("MinStartdEdx1stTrk",                    2.5);
    fMaxEnddEdx1stTrk        = pset.get<double>      ("MaxEnddEdx1stTrk",                      4.0);
    fDoHists                 = pset.get<bool>        ("FillHistograms",                      false);
    fDebug                   = pset.get<int>         ("Debug",                                   0);
}
    
void NuMuCCSelectionIIAlg::beginJob(art::ServiceHandle<art::TFileService>& tfs)
{
    
    return;
}
    
void NuMuCCSelectionIIAlg::produces(art::EDProducer* owner)
{
    fMyProducerModule = owner;
    fMyProducerModule->produces< art::Assns<recob::Vertex, recob::Track> >();
    fMyProducerModule->produces< art::Assns<recob::Vertex, recob::PFParticle> >();
}

    
bool NuMuCCSelectionIIAlg::findNeutrinoCandidates(art::Event & evt) const
{
    // Agreed convention is to ALWAYS output to the event store so get a pointer to our collection
    std::unique_ptr<art::Assns<recob::Vertex, recob::Track>>      vertexTrackAssociations(new art::Assns<recob::Vertex, recob::Track>);
    std::unique_ptr<art::Assns<recob::Vertex, recob::PFParticle>> vertexPFParticleAssociations(new art::Assns<recob::Vertex, recob::PFParticle>);

    // tracks
    art::Handle< std::vector<recob::Track> > trackListHandle;
    std::vector<art::Ptr<recob::Track> > tracklist;
    if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
      art::fill_ptr_vector(tracklist, trackListHandle);
    
    // vertices
    art::Handle< std::vector<recob::Vertex> > vtxListHandle;
    std::vector<art::Ptr<recob::Vertex> > vtxlist;
    if (evt.getByLabel(fVertexModuleLabel,vtxListHandle))
      art::fill_ptr_vector(vtxlist, vtxListHandle);
    
    // flashes
    art::Handle< std::vector<recob::OpFlash> > flashListHandle;
    std::vector<art::Ptr<recob::OpFlash> > flashlist;
    if (evt.getByLabel(fOpFlashModuleLabel, flashListHandle))
      art::fill_ptr_vector(flashlist, flashListHandle);

    // associations
    art::FindMany<anab::Calorimetry>  fmcal(trackListHandle, evt, fCalorimetryModuleLabel);

    //check the flash info
    double FlashPEmax=0;
    int NuFlashID=-1;
   
    for (size_t i = 0; i<flashlist.size(); ++i){
      if (flashlist[i]->TotalPE()>fPEThresh && flashlist[i]->Time()>fBeamMin && flashlist[i]->Time()<fBeamMax){
        if (flashlist[i]->TotalPE()>FlashPEmax){
          FlashPEmax = flashlist[i]->TotalPE();
          NuFlashID = i;
        }
      }
    }

    //Did not find the desired flash, return
    if (NuFlashID == -1){
      evt.put(std::move(vertexTrackAssociations));
      evt.put(std::move(vertexPFParticleAssociations));
      return false;
    }

    //Save basic track information in vectors
    std::vector<double> trkstartx(tracklist.size());
    std::vector<double> trkstarty(tracklist.size());
    std::vector<double> trkstartz(tracklist.size());
    std::vector<double> trkendx(tracklist.size());
    std::vector<double> trkendy(tracklist.size());
    std::vector<double> trkendz(tracklist.size());
    std::vector<double> trkstartdcosx(tracklist.size());
    std::vector<double> trkstartdcosy(tracklist.size());
    std::vector<double> trkstartdcosz(tracklist.size());
    std::vector<double> trkenddcosx(tracklist.size());
    std::vector<double> trkenddcosy(tracklist.size());
    std::vector<double> trkenddcosz(tracklist.size());
    std::vector<double> trklen(tracklist.size());
    double larStart[3];
    double larEnd[3];
    std::vector<double> trackStart;
    std::vector<double> trackEnd;
    for (size_t i = 0; i<tracklist.size(); ++i){
      trackStart.clear();
      trackEnd.clear();
      memset(larStart, 0, 3);
      memset(larEnd, 0, 3);
      tracklist[i]->Extent(trackStart,trackEnd); 
      tracklist[i]->Direction(larStart,larEnd);
      trkstartx[i]      = trackStart[0];
      trkstarty[i]      = trackStart[1];
      trkstartz[i]      = trackStart[2];
      trkendx[i]        = trackEnd[0];
      trkendy[i]        = trackEnd[1];
      trkendz[i]        = trackEnd[2];
      trkstartdcosx[i]  = larStart[0];
      trkstartdcosy[i]  = larStart[1];
      trkstartdcosz[i]  = larStart[2];
      trkenddcosx[i]    = larEnd[0];
      trkenddcosy[i]    = larEnd[1];
      trkenddcosz[i]    = larEnd[2];
      trklen[i]         = tracklist[i]->Length();
    }

    //Match each track with the selected flash
    std::vector<bool> trackflashmatch(tracklist.size());
    bool foundtrackflashmatch = false;

    //double TaggedFlashYCenter = flashlist[NuFlashID]->YCenter();
    double TaggedFlashZCenter = flashlist[NuFlashID]->ZCenter();
    for (size_t i = 0; i<tracklist.size(); ++i){
      double FlashTrackDis = 1e10;
      if ((trkstartz[i]<TaggedFlashZCenter && trkendz[i]>TaggedFlashZCenter)||
          (trkstartz[i]>TaggedFlashZCenter && trkendz[i]<TaggedFlashZCenter)){
        FlashTrackDis = 0;
      }
      else{
        FlashTrackDis = std::min(std::abs(trkstartz[i] - TaggedFlashZCenter),
                                 std::abs(trkendz[i] - TaggedFlashZCenter));
      }
      if (FlashTrackDis<fTrk2FlashDist){
        trackflashmatch[i] = true;
        foundtrackflashmatch = true;
      }
      else{
        trackflashmatch[i] = false;
      }
    }
    if (!foundtrackflashmatch) {
      if (fDebug) std::cout<<"Did not find any tracks matching flash."<<std::endl;
      evt.put(std::move(vertexTrackAssociations));
      evt.put(std::move(vertexPFParticleAssociations));    
      return false;
    }

    //Match tracks with vertices
    if (!vtxlist.size()) {
      if (fDebug) std::cout<<"No vertex found"<<std::endl;
      evt.put(std::move(vertexTrackAssociations));
      evt.put(std::move(vertexPFParticleAssociations));    
      return false;
    }

    std::vector<std::vector<int>> trkindex(vtxlist.size());
    std::vector<std::vector<bool>> fliptrack(vtxlist.size());
    for (size_t i = 0; i<vtxlist.size(); ++i){
      double xyz[3];
      vtxlist[i]->XYZ(xyz);
      if (!inFV(xyz[0], xyz[1], xyz[2])) continue;
        for (size_t j = 0; j<tracklist.size(); ++j){
          double vtxtrkStartDis = sqrt(pow(trkstartx[j]-xyz[0],2)+
                                       pow(trkstarty[j]-xyz[1],2)+
                                       pow(trkstartz[j]-xyz[2],2));
          double vtxtrkEndDis = sqrt(pow(trkendx[j]-xyz[0],2)+
                                     pow(trkendy[j]-xyz[1],2)+
                                     pow(trkendz[j]-xyz[2],2));
          double vtxtrkDis = std::min(vtxtrkStartDis, vtxtrkEndDis);
          if (vtxtrkDis<fMinTrk2VtxDist){
            trkindex[i].push_back(j);
            if (vtxtrkEndDis<vtxtrkStartDis){
              fliptrack[i].push_back(true);
            }
            else{
              fliptrack[i].push_back(false);
            }
          }
        }//Loope over all tracks
    }//Loop over all vertices
    
    //calculate average dE/dx near the track start and track end
    std::vector<double> trkStartdEdx(tracklist.size());
    std::vector<double> trkEnddEdx(tracklist.size());
    for (size_t i = 0; i<tracklist.size(); ++i){
      if (fmcal.isValid()){
        std::vector<const anab::Calorimetry*> calos = fmcal.at(i);
        int icalo = -1;
        int totalnhits = 0;
        for (size_t j = 0; j<calos.size(); ++j){
          if (int(calos[j]->dEdx().size())>totalnhits){
            icalo = j;
            totalnhits = calos[j]->dEdx().size();
          }
        }
        
        double sumdEdxStart=0;
        double sumdEdxEnd=0;
        
        int MaxHits=0;
        if(totalnhits>=20){
          MaxHits=10;
        }
        else if(totalnhits>0){
          MaxHits=totalnhits/2;
        }
        for(int ihit=0;ihit<MaxHits;ihit++){
          sumdEdxStart += calos[icalo]->dEdx()[ihit]*scaledEdx(calos[icalo]->XYZ()[ihit].X(), calos[icalo]->PlaneID().Plane, evt.isRealData());
          sumdEdxEnd += calos[icalo]->dEdx()[totalnhits-ihit-1]*scaledEdx(calos[icalo]->XYZ()[totalnhits-ihit-1].X(), calos[icalo]->PlaneID().Plane, evt.isRealData());
        }
        //        if (fDebug) std::cout<<trkxyz[i*3*2000*3+iplane*2000*3+0*3+0]<<" "
        //            <<trkxyz[i*3*2000*3+iplane*2000*3+0*3+1]<<" "
        //            <<trkxyz[i*3*2000*3+iplane*2000*3+0*3+2]<<" "
        //            <<trkxyz[i*3*2000*3+iplane*2000*3+(totalnhits-1)*3+0]<<" "
        //            <<trkxyz[i*3*2000*3+iplane*2000*3+(totalnhits-1)*3+1]<<" "
        //            <<trkxyz[i*3*2000*3+iplane*2000*3+(totalnhits-1)*3+2]<<" "
        //            <<trkstartx[i]<<" "<<trkstarty[i]<<" "<<trkstartz[i]<<" "
        //            <<trkendx[i]<<" "<<trkendy[i]<<" "<<trkendz[i]<<endl;
        if (icalo!=-1&&
            sqrt(pow(calos[icalo]->XYZ()[0].X()-trkstartx[i],2)+
                 pow(calos[icalo]->XYZ()[0].Y()-trkstarty[i],2)+
                 pow(calos[icalo]->XYZ()[0].Z()-trkstartz[i],2))>
            sqrt(pow(calos[icalo]->XYZ()[0].X()-trkendx[i],2)+
                 pow(calos[icalo]->XYZ()[0].Y()-trkendy[i],2)+
                 pow(calos[icalo]->XYZ()[0].Z()-trkendz[i],2))){
          std::swap(sumdEdxEnd, sumdEdxStart);
        }
        if(MaxHits>0 && sumdEdxEnd>0 && sumdEdxStart>0){
          trkStartdEdx[i] = sumdEdxStart/MaxHits;
          trkEnddEdx[i] = sumdEdxEnd/MaxHits;
        }
        if (fDebug) std::cout<<"Trkid = "<<i<<" MaxHits = "<<MaxHits<<" sumdEdxStart "<<sumdEdxStart<<" sumdEdxEnd "<<sumdEdxEnd<<" icalo "<<icalo<<std::endl;
      }
    }//Loop over tracks
    
    //Examine tracks around each vertex and select neutrino candidates
    std::vector<bool> nuvtx(vtxlist.size());
    std::vector<double> cosangle(vtxlist.size()); //angle between two longest tracks
    std::vector<double> dcosylong(vtxlist.size()); //dcosy of the longest track
    std::vector<double> trklen2nd(vtxlist.size()); //track length of the second longest track
    for (size_t i = 0; i<vtxlist.size(); ++i){
      
      //Check if there are tracks associated with this vertex
      if (!trkindex[i].size()) {
          //no tracks associated with vertex
        nuvtx[i] = false;
        if (fDebug) std::cout<<"ivtx = "<<i<<" no tracks associated with this vertex."<<std::endl;
        continue;
      }
      
      //Check if at least one track matches the flash
      bool flashmatch = false;
      for (size_t j = 0; j<trkindex[i].size(); ++j){
        if (trackflashmatch[trkindex[i][j]]){
          flashmatch = true;
        }
      }
      if (!flashmatch){
        //no tracks matched to the flash around this vertex
          nuvtx[i] = false;
          if (fDebug) std::cout<<"ivtx = "<<i<<" no tracks around the vertex matched to the flash."<<std::endl;
          continue;
      }
      
      //study if mult>=2
      if (trkindex[i].size()>1){
        //find two longest tracks
        int j0 = -1;
        int j1 = -1;
        double trklen0 = -1;
        double trklen1 = -1;
        //find the highest track
        int jhigh = -1;
        double highy= -1000;
        for (size_t j = 0; j<trkindex[i].size(); ++j){//Loop over all tracks around vertex
          if (trklen[trkindex[i][j]]>trklen0){
            j1 = j0;
            trklen1 = trklen0;
            j0 = j;
            trklen0 = trklen[trkindex[i][j]];
          }
          else if (trklen[trkindex[i][j]]>trklen1){
            j1 = j;
            trklen1 = trklen[trkindex[i][j]];
          }
          
          if (trkstarty[trkindex[i][j]]>highy){
            highy = trkstarty[trkindex[i][j]];
            jhigh = j;
          }
          if (trkendy[trkindex[i][j]]>highy){
            highy = trkendy[trkindex[i][j]];
            jhigh = j;
          }
        }//Loop over all tracks around vertex
        float dcosx0 = trkstartdcosx[trkindex[i][j0]];
        float dcosy0 = trkstartdcosy[trkindex[i][j0]];
        float dcosz0 = trkstartdcosz[trkindex[i][j0]];
        float dcosx1 = trkstartdcosx[trkindex[i][j1]];
        float dcosy1 = trkstartdcosy[trkindex[i][j1]];
        float dcosz1 = trkstartdcosz[trkindex[i][j1]];
        if (fliptrack[i][j0]){
          dcosx0 = trkenddcosx[trkindex[i][j0]];
          dcosy0 = trkenddcosy[trkindex[i][j0]];
          dcosz0 = trkenddcosz[trkindex[i][j0]];
        }
        if (fliptrack[i][j1]){
          dcosx1 = trkenddcosx[trkindex[i][j1]];
          dcosy1 = trkenddcosy[trkindex[i][j1]];
          dcosz1 = trkenddcosz[trkindex[i][j1]];
        }
        //cosine angle between two longest tracks, 1 indicates broken track
        cosangle[i] = std::abs(dcosx0*dcosx1+dcosy0*dcosy1+dcosz0*dcosz1);
        if (j0 == jhigh){
          dcosylong[i] = fliptrack[i][j0]?std::abs(trkstartdcosy[trkindex[i][j0]]):std::abs(trkenddcosy[trkindex[i][j0]]);
          trklen2nd[i] = trklen[trkindex[i][j1]];
        }
        if (fDebug) std::cout<<i<<" "<<trkindex[i].size()<<" "<<j0<<" "<<jhigh<<" "<<cosangle[i]<<" "<<dcosylong[i]<<" "<<trklen2nd[i]<<std::endl;
        if (cosangle[i]>fMaxCosineAngle) {
          nuvtx[i] = false;
          continue;
        }
        if (j0 == jhigh){
          if (dcosylong[i]>fMaxCosy1stTrk&&trklen2nd[i]<fMinTrackLen2ndTrk){
            nuvtx[i] = false;
            continue;
          }
        }
      }//Multi>1
      
      //If there are more than 2 tracks, accept the vertex 
      if (trkindex[i].size()>2){
        nuvtx[i] = true;
        if (fDebug) std::cout<<"ivtx = "<<i<<" track multiplicity = "<<trkindex[i].size()<<std::endl;
        continue;
      }
      
      //Single track
      if (trkindex[i].size()==1){
        nuvtx[i] = false; //will update later
        //Only select contained track
        int itrk = trkindex[i][0];
        //Track is fully contained
        if (inFV(trkstartx[itrk],
                 trkstarty[itrk],
                 trkstartz[itrk])&&
            inFV(trkendx[itrk],
                 trkendy[itrk],
                 trkendz[itrk])){
          if (std::abs(trkstartdcosy[itrk])>fMaxCosySingle) {
            nuvtx[i] = false;
            continue;
          }
          //At least 40 cm
          if (trklen[itrk]>fMinTrackLenSingle){
            double TrackLengthYNu = trklen[itrk]*std::abs(trkstartdcosy[itrk]);
            double LongSingleTrackdEdxRatio = -999;              
            if (trkstarty[itrk]>trkendy[itrk]){
              LongSingleTrackdEdxRatio = trkStartdEdx[itrk]/trkEnddEdx[itrk];
            }
            else{
              LongSingleTrackdEdxRatio = trkEnddEdx[itrk]/trkStartdEdx[itrk];
            }
            if (fDebug) std::cout<<"LongSingleTrackdEdxRatio = "<<LongSingleTrackdEdxRatio<<" TrackLengthYNu "<<TrackLengthYNu<<" trkStartdEdx "<<trkStartdEdx[itrk]<<" trkEnddEdx "<<trkEnddEdx[itrk]<<std::endl;
            if(LongSingleTrackdEdxRatio>fMindEdxRatioSingle || (TrackLengthYNu<=fMaxTrkLengthySingle && LongSingleTrackdEdxRatio<=fMindEdxRatioSingle)){
              nuvtx[i] = true;
            }
          }//At least 40 cm
        }//Track is contained
        if (fDebug) std::cout<<"ivtx = "<<i<<" single track "<<nuvtx[i]<<std::endl;
        continue;
        }//Single track
      
      //Multiplicity = 2
      if (trkindex[i].size()==2){
        bool isMichel = false;
        //if (fDebug) std::cout<<trklen[trkindex[i][0]]<<" "<<fliptrack[i][0]<<" "<<trkStartdEdx[trkindex[i][0]]<<" "<<trkEnddEdx[trkindex[i][0]]<<" "<<trkstarty[trkindex[i][0]]<<" "<<trkendy[trkindex[i][0]]<<" "<<trklen[trkindex[i][1]]<<" "<<fliptrack[i][1]<<" "<<trkStartdEdx[trkindex[i][1]]<<" "<<trkEnddEdx[trkindex[i][1]]<<" "<<trkstarty[trkindex[i][1]]<<" "<<trkendy[trkindex[i][1]]<<std::endl;
        double trkstartdedx0 = 0;
        double trkenddedx0 = 0;
        double trkendy0 = 0;
        double trklen1 = 0;
        if (trklen[trkindex[i][0]]>trklen[trkindex[i][1]]){//first track is longer
          if (!fliptrack[i][0]){
            trkstartdedx0 = trkStartdEdx[trkindex[i][0]];
            trkenddedx0 = trkEnddEdx[trkindex[i][0]];
            trkendy0 = trkendy[trkindex[i][0]];
            trklen1 = trklen[trkindex[i][1]];
          }
          else{
            trkstartdedx0 = trkEnddEdx[trkindex[i][0]];
            trkenddedx0 = trkStartdEdx[trkindex[i][0]];
            trkendy0 = trkstarty[trkindex[i][0]];
            trklen1 = trklen[trkindex[i][1]];
          }
        }//first track is longer
        else{//second track is longer
          if (!fliptrack[i][1]){
            trkstartdedx0 = trkStartdEdx[trkindex[i][1]];
            trkenddedx0 = trkEnddEdx[trkindex[i][1]];
            trkendy0 = trkendy[trkindex[i][1]];
            trklen1 = trklen[trkindex[i][0]];
          }
          else{
            trkstartdedx0 = trkEnddEdx[trkindex[i][1]];
            trkenddedx0 = trkStartdEdx[trkindex[i][1]];
            trkendy0 = trkstarty[trkindex[i][1]];
            trklen1 = trklen[trkindex[i][0]];
          }
        }//second track is longer
        if (((trkstartdedx0>trkenddedx0&&
              trkstartdedx0>fMinStartdEdx1stTrk&&trkenddedx0<fMaxEnddEdx1stTrk)||
               trkendy0>fDistToEdgeY)&&trklen1<fMinTrackLen2ndTrk)
          isMichel = true;
        
        if (isMichel){
          nuvtx[i] = false;
        }
        else{
          nuvtx[i] = true;
        }
        if (fDebug) std::cout<<"ivtx = "<<i<<" mul = 2 "<<nuvtx[i]<<std::endl;
        continue;
        }//Multiplicity = 2
      
    }//Loop over all vertices
    
    //Find the longest track
    int ivtx = -1;
    int itrk = -1;
    double longesttracklength = -1;
    for (size_t i = 0; i<vtxlist.size(); ++i){
      if (!nuvtx[i]) continue;
      for (size_t j = 0; j<trkindex[i].size(); ++j){
        if (fDebug) std::cout<<"ivtx = "<<i<<" trkid = "<<trkindex[i][j]<<" tracklen "<<trklen[trkindex[i][j]]<<" trackflashmatch "<<trackflashmatch[trkindex[i][j]]<<std::endl;
        if (trklen[trkindex[i][j]]>longesttracklength&&
            trackflashmatch[trkindex[i][j]]&&
            trklen[trkindex[i][j]]>fMinTrackLen){
          longesttracklength = trklen[trkindex[i][j]];
          ivtx = i;
          itrk = j;
        }
      }
    }//Loop over all vertices
    if (ivtx!=-1 && itrk!=-1){
      //outputfile[isample]<<run<<" "<<subrun<<" "<<event<<" "<<ivtx<<" "<<trkindex[ivtx][itrk]<<" "<<trkindex[ivtx].size()<<std::endl;
      util::CreateAssn(*fMyProducerModule, evt, tracklist[itrk], vtxlist[ivtx], *vertexTrackAssociations);
    }
    
    // Add associations to event.
    evt.put(std::move(vertexTrackAssociations));
    evt.put(std::move(vertexPFParticleAssociations));    
    
    return true;
}
    
bool NuMuCCSelectionIIAlg::inFV(double x, double y, double z) const
{
    double distInX = x - fGeometry->DetHalfWidth();
    double distInY = y;
    double distInZ = z - 0.5 * fGeometry->DetLength();
    
    if (std::abs(distInX) < fDistToEdgeX && std::abs(distInY) < fDistToEdgeY && std::abs(distInZ) < fDistToEdgeZ) return true;
    
    return false;
}

double NuMuCCSelectionIIAlg::scaledEdx(double x, int plane, bool isdata) const{
  double dEdx = 1.63;
  double p0_data[3] = {1.927, 2.215, 1.936};
  double p1_data[3] = {0.001495, 0.0001655, 0.001169};
  double p0_mc[3] = {1.843, 1.904, 1.918};
  double p1_mc[3] = {-0.0008329, -0.001357, -0.0007563};
  if (isdata){
    return dEdx/(p0_data[plane]+x*p1_data[plane]);
  }
  else{
    return  dEdx/(p0_mc[plane]+x*p1_mc[plane]);
  }
}


} // namespace
