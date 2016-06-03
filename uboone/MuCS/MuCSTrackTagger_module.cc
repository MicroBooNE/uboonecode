////////////////////////////////////////////////////////////////////////
/// \file  MuCSTagger_module.cc
/// \brief EDProducer for tagging tracks as MuCS.
///
/// \version $Id: MuCSTagger_module.cxx
/// \author  Matthew.Bass@physics.ox.ac.uk && 
////////////////////////////////////////////////////////////////////////
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iterator>

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/geo.h"

#include "lardataobj/RecoBase/Track.h"

#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "TVector3.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"


class MuCSTrackTagger : public art::EDProducer {
public:
  explicit MuCSTrackTagger(fhicl::ParameterSet const & p);
  virtual ~MuCSTrackTagger();

  void produce(art::Event & e) override;

  void beginJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;


private:

  
  double length(const art::Ptr<recob::Track> track); //< Length of reconstructed track, trajectory by trajectory.
  bool intersectsBoxes(const TVector3 & start, const TVector3& dir);

  std::string fTrackModuleLabel; //< Track label to find MuCS tags in
  std::vector<float> fMuCSTopBox, fMuCSBottomBox; //< Box Edge Positions (x1,x2,y1,y2,z1,z2)
  float fBoxExtension; //< Amount to extend acceptance for box interception [cm]
  unsigned int fDirFromNPoints; //< Number of points to use to determine track direction (0=use track end direction)
  float fMinTrackLength; //< Minimum length of track to consider [cm]
  
  //hists
  TH2F* fTopBoxPosHist;
  TH2F* fBottomBoxPosHist;
};

MuCSTrackTagger::MuCSTrackTagger(fhicl::ParameterSet const & p){
  this->reconfigure(p);
  // Call appropriate Produces<>() functions here.
  produces< std::vector<anab::CosmicTag> >();
  produces< art::Assns<recob::Track, anab::CosmicTag> >();
}

MuCSTrackTagger::~MuCSTrackTagger() {}

void MuCSTrackTagger::beginJob(){
  art::ServiceHandle<art::TFileService> tfs;
  
  fTopBoxPosHist     = tfs->make<TH2F>("topboxpos","TopBoxPositions;X;Z",100, 
                                  0.9*(fMuCSTopBox[0]-fBoxExtension), 1.1*(fMuCSTopBox[1]+fBoxExtension)
                                  ,100, 0.9*(fMuCSTopBox[4]-fBoxExtension), 1.1*(fMuCSTopBox[5]+fBoxExtension));

  fBottomBoxPosHist  = tfs->make<TH2F>("bottomboxpos","BottomBoxPositions;X;Z",100, 
                                  0.9*(fMuCSBottomBox[0]-fBoxExtension), 1.1*(fMuCSBottomBox[1]+fBoxExtension)
                                  ,100, 0.9*(fMuCSBottomBox[4]-fBoxExtension), 1.1*(fMuCSBottomBox[5]+fBoxExtension));
                                                                    
}


bool MuCSTrackTagger::intersectsBoxes(const TVector3 & start, const TVector3& dir){
  //return true if this trajector will intersect both boxes
  
  TVector3 newpTop, newpBottom;
  newpTop.SetXYZ(start.X() + (fMuCSTopBox[2]-start.Y())*dir.X()/dir.Y(),
                 fMuCSTopBox[2],
                 start.Z() + (fMuCSTopBox[2]-start.Y())*dir.Z()/dir.Y());
  newpBottom.SetXYZ(start.X() + (fMuCSBottomBox[2]-start.Y())*dir.X()/dir.Y(),
                 fMuCSBottomBox[2],
                 start.Z() + (fMuCSBottomBox[2]-start.Y())*dir.Z()/dir.Y());  
                 
  //populate hists regardless of intersection
  fTopBoxPosHist->Fill(newpTop.X(),newpTop.Z());
  fBottomBoxPosHist->Fill(newpBottom.X(),newpBottom.Z());

  
  if(newpTop.X() > fMuCSTopBox[0]-fBoxExtension && newpTop.X() < fMuCSTopBox[1]+fBoxExtension
     && newpTop.Z() > fMuCSTopBox[4]-fBoxExtension && newpTop.Z() < fMuCSTopBox[5]+fBoxExtension
     && newpBottom.X() > fMuCSBottomBox[0]-fBoxExtension && newpBottom.X() < fMuCSBottomBox[1]+fBoxExtension
     && newpBottom.Z() > fMuCSBottomBox[4]-fBoxExtension && newpBottom.Z() < fMuCSBottomBox[5]+fBoxExtension)
   return true;
   else
    return false;
}

void MuCSTrackTagger::produce(art::Event & e) {
  // Implementation of required member function here.

  std::unique_ptr< std::vector< anab::CosmicTag > >              cosmicTagTrackVector( new std::vector<anab::CosmicTag> );
  std::unique_ptr< art::Assns<recob::Track, anab::CosmicTag > >  assnOutCosmicTagTrack( new art::Assns<recob::Track, anab::CosmicTag>);

  art::Handle<std::vector<recob::Track> > Trk_h;
  e.getByLabel( fTrackModuleLabel, Trk_h );
  std::vector<art::Ptr<recob::Track> > TrkVec;
  art::fill_ptr_vector(TrkVec, Trk_h);

  for (auto trk: TrkVec){
    if(length(trk)<fMinTrackLength) continue;
    
    //choose highest edge as track start
    TVector3 start, end, startDir, endDir;
    if(trk->Vertex()[1]>trk->End()[1]){
      start=trk->Vertex();
      end=trk->End();
      startDir=trk->VertexDirection();
      endDir=trk->EndDirection();
    }else{
      start=trk->End();
      end=trk->Vertex();
      startDir=trk->EndDirection();
      endDir=trk->VertexDirection();
    }
    
    //find which end of the trajectory to use to get direction
    unsigned int pStart;
    TVector3 dir;
    int pSign;
    if(trk->LocationAtPoint(0)==start){
      pStart=0;
      pSign=1; //go forward for track direction
    }else if(trk->LocationAtPoint(trk->NumberTrajectoryPoints()-1)==start){
      pStart=trk->NumberTrajectoryPoints()-1;
      pSign=-1; //go backward for track direction
    }else{
      throw cet::exception("MuCSTrackTagger") << "Start seems to be in wrong position!\n";
    }
    
    if(fDirFromNPoints==0){
      //use reversed track start direction
      dir=-startDir;
    }else{
      //use diff between pstart and pstart+psign*(fDirFromNPoints-1)
      if(fDirFromNPoints>trk->NumberTrajectoryPoints())
        mf::LogInfo("MuCSTrackTagger") << "Track has too few trajectory points ("<<trk->NumberTrajectoryPoints()<<"), skipping it.\n";
      dir=(trk->LocationAtPoint(pStart) - trk->LocationAtPoint(pStart+pSign*(fDirFromNPoints-1))).Unit();
    }
    
    //find interesections and generate tags if appropriate
    bool btag=intersectsBoxes(start,dir);
    
    if (btag){
      cosmicTagTrackVector->emplace_back(-999.);
      util::CreateAssn(*this, e, *cosmicTagTrackVector, trk, *assnOutCosmicTagTrack );
    }
    
  }

 /*std::cout<<"\n"<<Trk_h->size()<<"\t"<<(*cosmicTagTrackVector).size();
   for(unsigned int f=0;f<Trk_h->size();f++){
   	std::cout<<"\n\t"<<f<<"\t"<<(*cosmicTagTrackVector)[f].CosmicScore()<<"\t"<<(*cosmicTagTrackVector)[f].CosmicType();
   }*/
 
  // e.put( std::move(outTracksForTags) );
  e.put( std::move(cosmicTagTrackVector) );
  e.put( std::move(assnOutCosmicTagTrack) );


} // end of produce

// Length of reconstructed track, trajectory by trajectory.
double MuCSTrackTagger::length(art::Ptr<recob::Track> track){
  double result = 0.;
  TVector3 disp = track->LocationAtPoint(0);
  int n = track->NumberTrajectoryPoints();
  for(int i = 1; i < n; ++i) {
    const TVector3& pos = track->LocationAtPoint(i);
    disp -= pos;
    result += disp.Mag();
    disp = pos;
  }
  return result;
}



void MuCSTrackTagger::reconfigure(fhicl::ParameterSet const & p) {

  fTrackModuleLabel = p.get< std::string >("TrackModuleLabel", "track");

  fMuCSBottomBox=p.get< std::vector< float > >("MuCSBottomBox");
  fMuCSTopBox=p.get< std::  vector< float > >("MuCSTopBox");
  if(fMuCSBottomBox.size()!=6 || fMuCSTopBox.size()!=6)
    throw cet::exception("MuCSTrackTagger") << "MuCSBottomBox or MuCSTopBox has wrong size!\n";
  
  fBoxExtension=p.get<float>("BoxExtension",0.);
  fDirFromNPoints=p.get<unsigned int>("DirFromNPoints",0);
  fMinTrackLength=p.get<float>("MinTrackLength",0.);
  
}

DEFINE_ART_MODULE(MuCSTrackTagger)
