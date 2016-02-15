
////////////////////////////////////////////////////////////////////////
//
//    trivia : Finding the time offset between the TPC 
//             and the Muon Counter System (MuCS), September 2015
//    author : Leonidas N. Kalousis
//    e-mail : kalousis@vt.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef MuCSDT_Module

#define MuCSDT_Module


#include "larsim/Simulation/SimChannel.h"

#include "larsim/Simulation/LArG4Parameters.h"


#include "lardata/DetectorInfoServices/LArPropertiesService.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"


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

#include "TFile.h"


#include <map>

#include <vector>

#include <algorithm>

#include <iostream>

#include <string>

#include <cmath>


using namespace std;


namespace MuCSDT 
{
  class MuCSDT : public art::EDAnalyzer 
  {
  public:
    
    explicit MuCSDT( fhicl::ParameterSet const& pset );
    
    virtual ~MuCSDT();
    
    void beginJob();
        
    void beginRun( const art::Run& run );
        
    void reconfigure( fhicl::ParameterSet const& pset );
        
    void analyze ( const art::Event& evt ); 
    
    void endJob();
    
  private:
    
    std::string fSwizzlerProducerLabel; 
    
    Int_t group = 0;
        
    Int_t trigID = 0;
    TH1F *hDT;
    Int_t run0;
    Int_t srun0;
    TFile *f1;
        
    Int_t seq;
    Float_t time_sec_low;
    Float_t time_sec_high;
    Float_t time_16ns_low;
    Float_t time_16ns_high;
    Double_t t0;
    
    TTree *my_tree;
    Int_t my_entries;
    
    Double_t previous_trigtime;
    Double_t t_start; 
    
    Double_t TOLER = 20.0;
    Double_t offset = -666.0;
        
  }; 
  
  MuCSDT::MuCSDT( fhicl::ParameterSet const& parameterSet )
    : EDAnalyzer(parameterSet)
  {
    this->reconfigure(parameterSet);
  
  }
  
  MuCSDT::~MuCSDT() 
  {}
    
  void MuCSDT::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    // hDT = tfs->make<TH1F>( "hDT", "", 10800, 0, 10800 );
    hDT = tfs->make<TH1F>( "hDT", "", 1000*10800, 0, 10800 );
    
    f1 = new TFile( Form( "/uboone/data/users/kalousis/MuCS/muons/mega_micro_ana_%d_0.333_0.root", group ), "read" );  
    
    if ( f1->IsZombie() ) 
      {
	cout << " - mucs file not existing ! " << endl;
	return;
	
      }
    
    else
      {
	my_tree = (TTree*)f1->Get( "preselected" );
	
	my_tree->SetBranchStatus( "*", 0 ); 
	my_tree->SetBranchStatus( "seq", 1 );
	my_tree->SetBranchStatus( "time_sec_low", 1 );
	my_tree->SetBranchStatus( "time_sec_high", 1 );
	my_tree->SetBranchStatus( "time_16ns_low", 1 );
	my_tree->SetBranchStatus( "time_16ns_high", 1 );
	my_tree->SetBranchStatus( "t0", 1 );
	
	my_tree->SetBranchAddress( "seq", &seq );
	my_tree->SetBranchAddress( "time_sec_low", &time_sec_low );
	my_tree->SetBranchAddress( "time_sec_high", &time_sec_high );
	my_tree->SetBranchAddress( "time_16ns_low", &time_16ns_low );
	my_tree->SetBranchAddress( "time_16ns_high", &time_16ns_high );
	my_tree->SetBranchAddress( "t0", &t0 );
	
	my_entries = my_tree->GetEntries();
	cout << " - events in mucs : " << my_entries << endl;
	cout << "" << endl;
	cout << "" << endl;
		
      }
    
  }
  
  void MuCSDT::beginRun( const art::Run& run )
  {}
  
  void MuCSDT::reconfigure( fhicl::ParameterSet const& p )
  {
    fSwizzlerProducerLabel = p.get< std::string >( "SwizzlerProducerLabel" );
    
    group = p.get< int >( "group" );
    
    return;
    
  }
    
  void MuCSDT::analyze( const art::Event& evt ) 
  {
    if ( trigID==0 )
      {
	cout << "" << endl;
	cout << " starting ... " << endl;
	cout << "" << endl;
	
	cout << " - group : " << group << endl;
	cout << "" << endl;
	cout << "" << endl;
	
	run0 = evt.run();
	srun0 = evt.subRun();

	previous_trigtime = 0.0;
	
      }
    
    Int_t event = evt.id().event();
    cout << "" << endl;
    cout << " - run : " << run0 << ", sub. run : " << srun0 << ", event : " << event << endl;
    cout << "" << endl;
    
    unsigned long long int tsval = evt.time().value();
    const unsigned long int mask32 = 0xFFFFFFFFUL;
    unsigned long int unix_time_stamp = ( tsval >> 32 ) & mask32;
    // unsigned long int llo = tsval & mask32;
    // TTimeStamp ts(unix_time_stamp, (int)llo);
    cout << " - unix timestamp : " << unix_time_stamp << endl;
    cout << "" << endl; 
    
    art::Handle< std::vector<raw::Trigger> > trigHandle;
    evt.getByLabel( fSwizzlerProducerLabel, trigHandle );
    std::vector< art::Ptr<raw::Trigger> > trigs;
    art::fill_ptr_vector( trigs, trigHandle );
    Double_t trigtime = ( trigs.at(0)->TriggerTime() )*1.0e-6;
    cout << " - trig. time : " << trigtime << ", diff : " << trigtime-previous_trigtime <<  ", " << trigs.size() << endl;
    cout << "" << endl; 
        
    if ( trigID==0 ) t_start = trigtime;
    Double_t t_rel = trigtime-t_start;
    cout << " - relative trig. time : " << t_rel << endl;
    cout << "" << endl; 
        
    for ( Int_t i=0; i<my_entries; i++ )
      {
	my_tree->GetEntry( i );
	
	Float_t DTunix = TMath::Abs( time_sec_high*65536.0+time_sec_low-unix_time_stamp );
		
	if ( DTunix<=TOLER ) 
	  {
	    // cout << " Gotcha !!! " << endl;
	    // cout << "" << endl;
	    // cout << " i : " << i << ", mucs unix timestamp : " << Form( "%.1f", time_sec_high*65536.0+time_sec_low ) << ", diff : " << DTunix << endl; 
	    // cout << "" << endl;
	    // getchar();
	    	    
	    Double_t tmucs = t0*1.0e-9;
	    hDT->Fill( tmucs-t_rel );
	    
	  }
		
      }
    
    previous_trigtime = trigtime;
    
    trigID++;
    return;
    
  }
  
  void MuCSDT::endJob()
  {
    cout << "" << endl; 
    cout << " - events processed : " << trigID << endl;
    cout << "" << endl;
    offset = hDT->GetXaxis()->GetBinCenter( hDT->GetMaximumBin() );
    // offset = hDT->GetXaxis()->GetBinUpEdge( hDT->GetMaximumBin() );
    cout << " - DT : " << Form( "%.6f", offset ) << endl;
    cout << "" << endl;
    cout << " ... ending ! " << endl;
    cout << "" << endl;
    
  }
  
  DEFINE_ART_MODULE( MuCSDT )
  
} 

#endif 

////////////////////////////////////////////////////////////////////////
//
//    The end !
//
////////////////////////////////////////////////////////////////////////
