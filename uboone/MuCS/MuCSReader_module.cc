
////////////////////////////////////////////////////////////////////////
//
//    trivia : A simple analyser to read merged data, September 2015
//    author : Leonidas N. Kalousis
//    e-mail : kalousis@vt.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef MuCSReader_Module

#define MuCSReader_Module


#include "larsim/Simulation/SimChannel.h"

#include "larsim/Simulation/LArG4Parameters.h"


#include "lardata/Utilities/LArProperties.h"

#include "lardata/Utilities/DetectorProperties.h"


#include "larcore/Geometry/Geometry.h"

#include "larcore/Geometry/OpDetGeo.h"

#include "SimulationBase/MCParticle.h"

#include "SimulationBase/MCTruth.h"

#include "larcore/SimpleTypesAndConstants/geo_types.h"


#include "lardata/RecoBase/Hit.h"

#include "larreco/RecoAlg/SpacePointAlg.h"

#include "lardata/RecoBase/Cluster.h"

#include "lardata/RecoBase/Track.h"

#include "lardata/RecoBase/PFParticle.h"  

#include "lardata/RecoBase/SpacePoint.h"  

#include "lardata/RecoBase/OpHit.h"  

#include "lardata/RecoBase/OpFlash.h"  

#include "lardata/RecoObjects/BezierTrack.h"

#include "larreco/RecoAlg/TrackMomentumCalculator.h"

#include "lardata/RawData/TriggerData.h"


#include "art/Framework/Core/EDAnalyzer.h"

#include "art/Framework/Principal/Event.h"

#include "art/Framework/Principal/Handle.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
 
#include "art/Framework/Services/Optional/TFileService.h"

#include "art/Framework/Core/ModuleMacros.h"

#include "art/Framework/Core/FindManyP.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "fhiclcpp/ParameterSet.h"

// #include "EventDisplay/HeaderDrawer.h"


#include "TMath.h"

#include "TH1.h"

#include "TAxis.h"

#include "TH2.h"

#include "TTree.h"

#include "TLorentzVector.h"

#include "TVector3.h"


#include <map>

#include <vector>

#include <algorithm>

#include <iostream>

#include <string>

#include <cmath>


#include "MuCSData.h"

using namespace std;


namespace MuCSReader 
{
  class MuCSReader : public art::EDAnalyzer 
  {
  public:
    
    explicit MuCSReader( fhicl::ParameterSet const& pset );
    
    virtual ~MuCSReader();
    
    void beginJob();
        
    void beginRun( const art::Run& run );
        
    void reconfigure( fhicl::ParameterSet const& pset );
        
    void analyze ( const art::Event& evt ); 
    
    void endJob();
    
  private:
    
    std::string fMergerProducerLabel; 
    Int_t group = 0;
    Int_t trigID = 0;
    Int_t run0;
    Int_t srun0;
    
  }; 
  
  MuCSReader::MuCSReader( fhicl::ParameterSet const& parameterSet )
    : EDAnalyzer(parameterSet)
  {
    this->reconfigure(parameterSet);
  
  }
  
  MuCSReader::~MuCSReader() 
  {}
    
  void MuCSReader::beginJob()
  {}
  
  void MuCSReader::beginRun( const art::Run& run )
  {}
  
  void MuCSReader::reconfigure( fhicl::ParameterSet const& p )
  {
    fMergerProducerLabel = p.get< std::string >( "MergerProducerLabel" );
        
    return;
    
  }
    
  void MuCSReader::analyze( const art::Event& evt ) 
  {
    if ( trigID==0 )
      {
	cout << "" << endl;
	cout << " starting ... " << endl;
	cout << "" << endl;
		
	run0 = evt.run();
	srun0 = evt.subRun();
		
      }
    
    Int_t event = evt.id().event();
    cout << "" << endl;
    cout << " - run : " << run0 << ", sub. run : " << srun0 << ", event : " << event << endl;
    cout << "" << endl;
    /*
    art::Handle< std::vector<MuCS::MuCSData> > mucsHandle;
    evt.getByLabel( fMergerProducerLabel, mucsHandle );
    std::vector< art::Ptr<MuCS::MuCSData> > mucs;
    art::fill_ptr_vector( mucs, mucsHandle );
    */
    std::vector< art::Handle< std::vector<MuCS::MuCSData> > > mucslist;
    evt.getManyByType( mucslist );
    art::Handle< std::vector<MuCS::MuCSData> > mucs = mucslist[0]; 
    cout << mucs->size() << endl;
    
    Float_t time0 = mucs->at(0).T0();
    cout << time0 << endl;
    
    Int_t tot1 = mucs->at(0).Hits1().size();
    cout <<  tot1 << endl; cout << " hits 1 : ";
    for ( Int_t j=0; j<tot1; j++ ) cout << mucs->at(0).Hits1().at(j) << ", ";
    cout << "" << endl;
    
    Int_t tot2 = mucs->at(0).Hits2().size();
    cout <<  tot2 << endl;  cout << " hits 2 : ";
    for ( Int_t j=0; j<tot2; j++ ) cout << mucs->at(0).Hits2().at(j) << ", ";
    cout << "" << endl;
    
    Int_t tot3 = mucs->at(0).Hits3().size();
    cout <<  tot3 << endl;  cout << " hits 3 : ";
    for ( Int_t j=0; j<tot3; j++ ) cout << mucs->at(0).Hits3().at(j) << ", ";
    cout << "" << endl;
    
    Int_t tot7 = mucs->at(0).Hits7().size();
    cout <<  tot7 << endl;  cout << " hits 7 : ";
    for ( Int_t j=0; j<tot7; j++ ) cout << mucs->at(0).Hits7().at(j) << ", ";
    cout << "" << endl;
    
    cout << "*" << endl;
    cout << "" << endl;
    cout << "" << endl;
    
    cout << " ADC1 : ";
    for ( Int_t j=0; j<24; j++ ) cout << mucs->at(0).ADC1().at(j) << ", ";
    cout << "" << endl;
    
    cout << " ADC2 : ";
    for ( Int_t j=0; j<24; j++ ) cout << mucs->at(0).ADC2().at(j) << ", ";
    cout << "" << endl;
    
    cout << " ADC3 : ";
    for ( Int_t j=0; j<24; j++ ) cout << mucs->at(0).ADC3().at(j) << ", ";
    cout << "" << endl;
    
    cout << " ADC7 : ";
    for ( Int_t j=0; j<24; j++ ) cout << mucs->at(0).ADC7().at(j) << ", ";
    cout << "" << endl;
    
    getchar();
    cout << "*" << endl;
    cout << "" << endl;
    cout << "" << endl;
    
    
    
    trigID++;
    return;
    
  }
  
  void MuCSReader::endJob()
  {
    cout << "" << endl; 
    cout << " - events processed : " << trigID << endl;
    cout << "" << endl;
    cout << " ... ending ! " << endl;
    cout << "" << endl;
    
  }
  
  DEFINE_ART_MODULE( MuCSReader )
  
} 

#endif 

////////////////////////////////////////////////////////////////////////
//
//    The end !
//
////////////////////////////////////////////////////////////////////////
