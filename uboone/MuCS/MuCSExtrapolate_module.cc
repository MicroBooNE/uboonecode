
////////////////////////////////////////////////////////////////////////
// Class:       MuCSReco
// Module Type: producer
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//
//    trivia : Reco stage of the MuCSMerger process. 
//             Adds angle and position info to hit patterns based on hits.
//    author : Matt Bass
//    e-mail : Matthew.Bass@physics.ox.ac.uk
//
////////////////////////////////////////////////////////////////////////

#ifndef MUCSRECO_H
#define MUCSRECO_H 

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TText.h"
#include "TTimeStamp.h"

#include <memory>
#include <iostream>
#include "vector"
#include "lardata/RawData/TriggerData.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include <sqlite3.h> 
#include "MuCSData.h"
#include "MuCSRecoData.h"

using namespace std;

struct dbfields {
 Int_t matches;
 Int_t nentries;
 Float_t p;
 Float_t q;
 Float_t p_rms;
 Float_t q_rms;
};
  
class MuCSReco;

class MuCSReco : public art::EDProducer {
public:
  explicit MuCSReco( fhicl::ParameterSet const &pset );
  //virtual ~MuCSReco();
  
  void reconfigure( fhicl::ParameterSet const &pset );
  void produce( art::Event &evt ) override;
      
private:
  void getStripNumberFromPMTHit(int pmt, int hit, int &layer, int &strip);
  double getStripDims(int layer,int dim, int strip);
  double getAvgHitPosition(const int pmt, const std::vector<int> hits, const int dim);
  void execFill(const int pmta, const std::vector<int> hitsa, const int pmtb, const std::vector<int> hitsb, dbfields &ldbfields, const int dim);
  unsigned int getChannelFromHit(const int hit,const std::vector<int> lvec);
  double getLayerMiddle(int layer,int dim);

  //fcl parameters
  std::vector<int> fHitMap1, fHitMap2, fHitMap3, fHitMap7; //< Maps from PMT channels to strip positions
  unsigned int fBaseLayer1, fBaseLayer2, fBaseLayer3, fBaseLayer7; //<Index (top to bottom) of the topmost position for this PMT
  std::vector<double> fLayerDims; //<Dimensions of each layer (linearized multi-dim array layers vs dims (x1,x2,y1,y2,z1,z2))
  std::vector<int> fLayerDirections; //<Directions (1 or -1) of increasing strip number for each layer
  double fStripWidth; //<Width of the strips used to map strip number to position along the layer
  
};

void MuCSReco::reconfigure( fhicl::ParameterSet const &p ){
    fHitMap1 = p.get< std::vector<int> >( "HitMap1" );
    fHitMap2 = p.get< std::vector<int> >( "HitMap2" );
    fHitMap3 = p.get< std::vector<int> >( "HitMap3" );
    fHitMap7 = p.get< std::vector<int> >( "HitMap7" );
    fBaseLayer1 = p.get< unsigned int >( "BaseLayer1" );
    fBaseLayer2 = p.get< unsigned int >( "BaseLayer2" );
    fBaseLayer3 = p.get< unsigned int >( "BaseLayer3" );
    fBaseLayer7 = p.get< unsigned int >( "BaseLayer7" );
    fLayerDims = p.get< std::vector<double> >( "LayerDims" );
    fLayerDirections = p.get< std::vector<int> >( "LayerDirections" );
    fStripWidth = p.get< double >("StripWidth");
    return;
}

MuCSReco::MuCSReco( fhicl::ParameterSet const &pset ){
  this->reconfigure( pset );
  
  produces< std::vector<MuCS::MuCSRecoData> >();  
  
}

unsigned int MuCSReco::getChannelFromHit(const int hit,const std::vector<int> lvec){
  unsigned int pos = std::find(lvec.begin(), lvec.end(), hit) - lvec.begin();
  
  if(pos >= lvec.size())
    throw cet::exception("MuCSExtrapolate") << "Invalid hit to channel map!\n";
    
  return pos;
}

void MuCSReco::getStripNumberFromPMTHit(int pmt, int hit, int &layer, int &strip){
  //use pmt hit map vectors to return layer number and a hit
  int baselayer=0;
  std::vector<int> lvec;
  if (pmt==7){
    baselayer=fBaseLayer7;
    lvec=fHitMap7;
  }else if (pmt==3){
    baselayer=fBaseLayer3;
    lvec=fHitMap3;
  }else if (pmt==2){
    baselayer=fBaseLayer2;
    lvec=fHitMap2;
  }else if (pmt==1){
    baselayer=fBaseLayer1;
    lvec=fHitMap1;
  }else
    throw cet::exception("MuCSExtrapolate") << "Invalid pmt number!\n";

  unsigned int chan=getChannelFromHit(hit,lvec);
  layer=baselayer+(chan>11 ? 1 : 0);
  strip=chan-(chan>11 ? 12: 0);
  //std::cout<<"getStripNumberFromPMTHit baselayer, pmt, hit, layer, chan, strip, lvec0="<<baselayer<<","<<pmt<<","<<hit<<","<<layer<<","<<chan<<","<<strip<<","<<lvec[0]<<std::endl;
  
  if (layer<0 || layer>7 || strip<0 ||strip>12)
    throw cet::exception("MuCSExtrapolate") << "Strip/layer dimension error!\n";
}


double MuCSReco::getStripDims(int layer,int dim, int strip){
  //use layer geometry vectors to get avg strip dimensions
  unsigned int i=layer*6+dim*2;
  if (i+1 > fLayerDims.size())
    throw cet::exception("MuCSExtrapolate") << "Invalid layer dimensions index!\n";

  unsigned int offset=(fLayerDirections[layer]==1?0:1);
  
  double low = fLayerDims[i+offset]+ fLayerDirections[layer]*fStripWidth*strip;
  double high = fLayerDims[i+offset]+ fLayerDirections[layer]*fStripWidth*(strip+1);
  //std::cout<<"got strip dims i, layer, dim, strip, low, high="<<i<<","<<layer<<","<<dim<<","<<strip<<","<<low<<","<<high<<std::endl;
  return (low+high)/2;
  
}

double MuCSReco::getLayerMiddle(int layer,int dim){
  //use layer geometry vectors to get avg strip dimensions
  unsigned int i=layer*6+dim*2;
  if (i+1 > fLayerDims.size())
    throw cet::exception("MuCSExtrapolate") << "Invalid layer dimensions index!\n";

  double low = fLayerDims[i];
  double high = fLayerDims[i+1];
  return (low+high)/2;
  
}

double MuCSReco::getAvgHitPosition(const int pmt, const std::vector<int> hits, const int dim){
  //compute top position
  double avgdsum=0.;
  unsigned int avgdcount=0;
  for (int hit: hits){
    int layer=0,strip=0;
    getStripNumberFromPMTHit(pmt,hit,layer,strip);
    avgdsum+=getStripDims(layer,dim,strip);
    avgdcount++;
    //std::cout<<"found pmt,hit,layer,strip="<<pmt<<","<<hit<<","<<layer<<","<<strip<<std::endl;
  }
  
  if (avgdcount==0)
    throw cet::exception("MuCSExtrapolate") << "ERROR! No hits found!\n";
    
  return avgdsum/avgdcount;
}

void MuCSReco::execFill(const int pmta, const std::vector<int> hitsa, const int pmtb, const std::vector<int> hitsb, dbfields &ldbfields, const int dim){
  //map hits to an average position in each layer in the specified dim(ension)
  //hitsa is the top layer set, hitsb is the bottom layer set
  //compute start and angle and put into ldbfields
  

  double p1=getAvgHitPosition(pmta, hitsa, dim);
  double p2=getAvgHitPosition(pmtb, hitsb, dim);
  
  double deltaY=getLayerMiddle(4,1)-getLayerMiddle(0,1); //distance between boxes
  
  ldbfields.p=p1;
  ldbfields.q=atan2(deltaY,(p2-p1)); 

  if (hitsa.size()>hitsb.size())
    ldbfields.matches=hitsa.size();
  else
    ldbfields.matches=hitsb.size();
  
  //std::cout<<"found p1,p2,deltaY="<<p1<<","<<p2<<","<<deltaY<<std::endl;
  //std::cout<<"found q="<<ldbfields.q<<std::endl;
}

void MuCSReco::produce( art::Event &evt ){

  //get MuCSData object
  std::vector< art::Handle< std::vector<MuCS::MuCSData> > > mucslist;
  evt.getManyByType( mucslist );
  art::Handle< std::vector<MuCS::MuCSData> > mucs = mucslist[0]; 
  
  dbfields xfields, zfields;
  Float_t xq=0.,xq_rms=0.,x=0.,x_rms=0.,zq=0.,zq_rms=0.,z=0.,z_rms=0.,y=0.;
  Int_t xmatches=0,zmatches=0;
  //only fill fields if there are hits
  if (mucs->at(0).Hits1().size()>0 && mucs->at(0).Hits2().size()>0 && mucs->at(0).Hits3().size()>0 && mucs->at(0).Hits7().size()>0){
    execFill(3,mucs->at(0).Hits3(),1,mucs->at(0).Hits1(),xfields,0);
    execFill(7,mucs->at(0).Hits7(),2,mucs->at(0).Hits2(),zfields,2);
    xq=xfields.q;
    x=xfields.p;
    zq=zfields.q;
    z=zfields.p;
    xmatches=xfields.matches;
    zmatches=zfields.matches;
    y=(xmatches>0 && zmatches>0) ? getLayerMiddle(0,1) : 0.; //only populate if there were database entries for this field

  }

  //std::cout<<"found x,y,z="<<x<<","<<y<<","<<z<<std::endl;
  //std::cout<<"found qxz,qyz="<<xq<<","<<zq<<std::endl;
  //now create and populate the MuCSReco object
  std::unique_ptr< std::vector<MuCS::MuCSRecoData> > mucsrecocol(new std::vector<MuCS::MuCSRecoData>);
  MuCS::MuCSRecoData mucsrecoevt( xq, xq_rms, x, x_rms, zq, zq_rms, z, z_rms, y, xmatches, zmatches );
  mucsrecocol->push_back( mucsrecoevt );
  evt.put( std::move( mucsrecocol ) );
  
}

DEFINE_ART_MODULE( MuCSReco )

#endif

////////////////////////////////////////////////////////////////////////
//
//    The end !
//
////////////////////////////////////////////////////////////////////////
