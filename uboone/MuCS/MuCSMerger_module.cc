
////////////////////////////////////////////////////////////////////////
// Class:       MuCSMerger
// Module Type: producer
// File:        MuCSMerger_module.cc
//
// Generated at Wed May 20 14:08:15 2015 by Leonidas N. Kalousis using artmod
// from cetpkgsupport v1_08_05.
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//
//    trivia : The main routine developed to merge data from the TPC 
//             and the Muon Counter System (MuCS), September 2015
//    author : Leonidas N. Kalousis
//    e-mail : kalousis@vt.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef MUCSMERGER_H
#define MUCSMERGER_H 

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardataobj/RawData/TriggerData.h"


#include <memory>
#include <iostream>
#include "vector"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include "MuCSData.h"
#include "MuCSDTOffset.h"
#include "ifdh.h"  //to handle flux files

using namespace std;

class MuCSMerger;

class MuCSMerger : public art::EDProducer {
public:
  explicit MuCSMerger( fhicl::ParameterSet const &pset );
  virtual ~MuCSMerger();
  
  void reconfigure( fhicl::ParameterSet const &pset ); // override;
  void produce( art::Event &evt ) override;
  void beginRun(art::Run& run);
  
private:
  
  std::string fSwizzlerProducerLabel; 
  std::string fMuCSFile;
  
  Int_t trigID = 0;
  TH1F *hDT;
  Int_t run0;
  Int_t srun0;
  TFile *fTFMuCSData;
  
  Int_t seq;
  Float_t time_sec_low;
  Float_t time_sec_high;
  Float_t time_16ns_low;
  Float_t time_16ns_high;
  Double_t t0;
  
  Double_t fadc1[24];
  Double_t fadc2[24];
  Double_t fadc3[24];
  Double_t fadc7[24];

  std::vector<Int_t> *fhits1 = new std::vector<Int_t>;
  std::vector<Int_t> *fhits2 = new std::vector<Int_t>;
  std::vector<Int_t> *fhits3 = new std::vector<Int_t>;
  std::vector<Int_t> *fhits7 = new std::vector<Int_t>;

  TTree *my_tree;
  Int_t my_entries;
  
  Double_t previous_trigtime;
  Double_t t_start; 
  
  Double_t TOLER = 20.0;
  Double_t offset = -666.0;
  Double_t TOLER2 = 2.5*0.001;

  ifdh_ns::ifdh* fIFDH=0; ///< For MuCS data file retrieval
};

void MuCSMerger::reconfigure( fhicl::ParameterSet const &p ){
  fSwizzlerProducerLabel = p.get< std::string >( "SwizzlerProducerLabel" );
  fMuCSFile = p.get< std::string >( "MuCSFile" );
  return;
}

MuCSMerger::MuCSMerger( fhicl::ParameterSet const &pset ){
  this->reconfigure( pset );
  produces< std::vector<MuCS::MuCSData> >();  
}

MuCSMerger::~MuCSMerger() {
  //close MuCS data root file
  fTFMuCSData->Close();
  //cleanup temp files
  fIFDH->cleanup();  
}

void MuCSMerger::beginRun(art::Run& run){
  //get offset from MuCSDTOffset stored in run
  std::vector< art::Handle< std::vector<MuCS::MuCSDTOffset> > > dtcol;
  run.getManyByType( dtcol );
  art::Handle< std::vector<MuCS::MuCSDTOffset> > tMuCSDTOffset = dtcol[0];   
  
  offset = tMuCSDTOffset->at(0).getoffset();
  cout << " - DT : " << Form( "%.6f", offset ) << endl;
  cout << "" << endl;
  
  //fetch MuCS data file using ifdh
  if ( ! fIFDH ) fIFDH = new ifdh_ns::ifdh;
  mf::LogInfo("MuCSMerger") << "Fetching: "<<fMuCSFile<<"\n";
  std::string fetchedfile(fIFDH->fetchInput(fMuCSFile));
  mf::LogInfo("MuCSMerger") << "Fetched; local path: "<<fetchedfile<<"\n";    

  fTFMuCSData = new TFile( Form( fetchedfile.c_str() ), "read" );  
  
  if ( fTFMuCSData->IsZombie() ) {
    cout << " - mucs file not existing ! " << endl;
    return;
  }else{
    my_tree = (TTree*)fTFMuCSData->Get( "preselected" );
    
    my_tree->SetBranchStatus( "*", 0 ); 
    my_tree->SetBranchStatus( "seq", 1 );
    my_tree->SetBranchStatus( "time_sec_low", 1 );
    my_tree->SetBranchStatus( "time_sec_high", 1 );
    my_tree->SetBranchStatus( "time_16ns_low", 1 );
    my_tree->SetBranchStatus( "time_16ns_high", 1 );
    my_tree->SetBranchStatus( "t0", 1 );
    
    my_tree->SetBranchStatus( "ADC1", 1 );
    my_tree->SetBranchStatus( "ADC2", 1 );
    my_tree->SetBranchStatus( "ADC3", 1 );
    my_tree->SetBranchStatus( "ADC7", 1 );
    
    my_tree->SetBranchStatus( "hits1", 1 );
    my_tree->SetBranchStatus( "hits2", 1 );
    my_tree->SetBranchStatus( "hits3", 1 );
    my_tree->SetBranchStatus( "hits7", 1 );
    
    my_tree->SetBranchAddress( "seq", &seq );
    my_tree->SetBranchAddress( "time_sec_low", &time_sec_low );
    my_tree->SetBranchAddress( "time_sec_high", &time_sec_high );
    my_tree->SetBranchAddress( "time_16ns_low", &time_16ns_low );
    my_tree->SetBranchAddress( "time_16ns_high", &time_16ns_high );
    my_tree->SetBranchAddress( "t0", &t0 );
    
    my_tree->SetBranchAddress( "ADC1", &fadc1 );
    my_tree->SetBranchAddress( "ADC2", &fadc2 );
    my_tree->SetBranchAddress( "ADC3", &fadc3 );
    my_tree->SetBranchAddress( "ADC7", &fadc7 );
    
    my_tree->SetBranchAddress( "hits1", &fhits1 );
    my_tree->SetBranchAddress( "hits2", &fhits2 );
    my_tree->SetBranchAddress( "hits3", &fhits3 );
    my_tree->SetBranchAddress( "hits7", &fhits7 );
    
    my_entries = my_tree->GetEntries();
    cout << " - events in mucs : " << my_entries << endl;
    cout << "" << endl;
    cout << "" << endl;
  }
}

void MuCSMerger::produce( art::Event &evt ){
  if ( trigID==0 ){
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
  
  std::unique_ptr< std::vector<MuCS::MuCSData> > mucsdatacol(new std::vector<MuCS::MuCSData>);
  Float_t time0 = -1.0;
  Float_t adc1[24], adc2[24], adc3[24], adc7[24];
  for ( Int_t j=0; j<24; j++ ) { adc1[j]=0; adc2[j]=0; adc3[j]=0; adc7[j]=0; }
  std::vector<Int_t> hits1, hits2, hits3, hits7;
  hits1.clear(); hits2.clear(); hits3.clear(); hits7.clear();
  
  Int_t ntimes=0;
  
  Int_t dtmin=10000;
  
  for ( Int_t i=0; i<my_entries; i++ ){
    my_tree->GetEntry( i );
    Float_t DTunix = TMath::Abs( time_sec_high*65536.0+time_sec_low-unix_time_stamp );
    if ( DTunix<=TOLER ){
      Double_t tmucs = t0*1.0e-9;
      Double_t dt0 =  tmucs-t_rel;
      Double_t dt = dt0-offset; 
      
      if ( TMath::Abs(dt)<TOLER2 ){
        cout << " Gotcha !!! " << endl;
        cout << "" << endl;
        cout << " i : " << i << ", mucs unix timestamp : " << Form( "%.1f", time_sec_high*65536.0+time_sec_low ) << ", diff : " << DTunix << endl; 
        cout << "" << endl;
        cout << " - mucs t0 : " << tmucs << ", " << "diff : " << dt << endl;
        cout << "" << endl;
        
        if ( TMath::Abs( dt )<TMath::Abs( dtmin ) ){
            time0 = tmucs;
            
            for ( Int_t j=0; j<24; j++ ) { 
              adc1[j]=fadc1[j]; adc2[j]=fadc2[j]; // cout << fadc1[j] << ", " << fadc2[j] << endl;
              adc3[j]=fadc3[j]; adc7[j]=fadc7[j]; // cout << fadc3[j] << ", " << fadc7[j] << endl;
            }
            
            Int_t tot1 = fhits1->size(); // cout << tot1 << endl;
            Int_t tot2 = fhits2->size(); // cout << tot2 << endl;
            Int_t tot3 = fhits3->size(); // cout << tot3 << endl;
            Int_t tot7 = fhits7->size(); // cout << tot7 << endl;
                    
            hits1.clear();
            for ( Int_t j=0; j<tot1; j++ ) hits1.push_back( fhits1->at(j) );
            
            hits2.clear();
            for ( Int_t j=0; j<tot2; j++ ) hits2.push_back( fhits2->at(j) );
            
            hits3.clear();
            for ( Int_t j=0; j<tot3; j++ ) hits3.push_back( fhits3->at(j) );
            
            hits7.clear();
            for ( Int_t j=0; j<tot7; j++ ) hits7.push_back( fhits7->at(j) );
                
            ntimes++;
            if ( ntimes>=2 ) { cout << " - MULTIPLE PAIRS, PROPERLY TREATED !!! " << endl; } // getchar(); }
            
            dtmin = dt;
        }
      }
    }
  }
  
  MuCS::MuCSData mucsevt( time0, adc1, adc2, adc3, adc7, hits1, hits2, hits3, hits7 ); 
  mucsdatacol->push_back( mucsevt );
  evt.put( std::move( mucsdatacol ) );
  
  previous_trigtime = trigtime;
    
  trigID++;
  return;
    
}

DEFINE_ART_MODULE( MuCSMerger )

#endif
