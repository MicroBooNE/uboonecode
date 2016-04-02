
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

//larsoft includes
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/SimpleTypesAndConstants/geo_types.h"
#include "lardata/RawData/TriggerData.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

//ROOT includes
#include "TMath.h"
#include "TH1.h"
#include "TAxis.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"

//C includes
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>

#include "MuCSDTOffset.h"
#include "ifdh.h"  //to handle flux files

using namespace std;
namespace MuCSDT 
{
  class MuCSDT : public art::EDProducer
  {
  public:
    explicit MuCSDT( fhicl::ParameterSet const& pset );
    virtual ~MuCSDT();
    void beginJob();
    void beginRun( art::Run& run );
    void endRun( art::Run& run );
    void reconfigure( fhicl::ParameterSet const& pset );
    void produce ( art::Event& evt ); 
    void endJob();
    
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
    
    TTree *my_tree;
    Int_t my_entries;
    
    Double_t previous_trigtime;
    Double_t t_start; 
    
    Double_t TOLER = 20.0;
    Double_t offset = -666.0;
    
    ifdh_ns::ifdh* fIFDH=0; ///< For MuCS data file retrieval
  }; 
  
  MuCSDT::MuCSDT( fhicl::ParameterSet const& pset ){
    this->reconfigure(pset);
    produces< std::vector<MuCS::MuCSDTOffset>, art::InRun >();
  }
  
  MuCSDT::~MuCSDT() {
    //close MuCS data root file
    fTFMuCSData->Close();
    //cleanup temp files
    fIFDH->cleanup();  
  }
    
  void MuCSDT::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    // hDT = tfs->make<TH1F>( "hDT", "", 10800, 0, 10800 );
    hDT = tfs->make<TH1F>( "hDT", "", 1000*10800, 0, 10800 );
    
    //fetch MuCS data file using ifdh
    if ( ! fIFDH ) fIFDH = new ifdh_ns::ifdh;
    mf::LogInfo("CorsikaGen") << "Fetching: "<<fMuCSFile<<"\n";
    std::string fetchedfile(fIFDH->fetchInput(fMuCSFile));
    mf::LogInfo("CorsikaGen") << "Fetched; local path: "<<fetchedfile<<"\n";    
    
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
  
  void MuCSDT::beginRun( art::Run& run ){}
  void MuCSDT::endRun( art::Run& run ){
    std::unique_ptr< std::vector<MuCS::MuCSDTOffset> > dtcol(new std::vector<MuCS::MuCSDTOffset>);
    MuCS::MuCSDTOffset dt( hDT->GetXaxis()->GetBinCenter( hDT->GetMaximumBin()  ));
    dtcol->push_back(dt);
    run.put( std::move( dtcol ) );
    //store found offset onto run
    //std::unique_ptr<std::vector<double>> MuCSDTOffsets(new std::vector<double>);
    //MuCSDTOffsets->push_back(hDT->GetXaxis()->GetBinCenter( hDT->GetMaximumBin() ));
    //*MuCSDTOffset=hDT->GetXaxis()->GetBinCenter( hDT->GetMaximumBin() );
    //run.put(std::move(MuCSDTOffsets),"MuCSDTOffsets");
  }
   
  void MuCSDT::reconfigure( fhicl::ParameterSet const& p ){
    fSwizzlerProducerLabel = p.get< std::string >( "SwizzlerProducerLabel" );
    fMuCSFile = p.get< std::string >( "MuCSFile" );
    return;
  }
    
  void MuCSDT::produce( art::Event& evt ){
    if ( trigID==0 ){
      cout << "" << endl;
      cout << " starting ... " << endl;
      cout << "" << endl;

      cout << " - MuCSFile : " << fMuCSFile << endl;
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
        
    for ( Int_t i=0; i<my_entries; i++ ){
      my_tree->GetEntry( i );
  
      Float_t DTunix = TMath::Abs( time_sec_high*65536.0+time_sec_low-unix_time_stamp );
    
      if ( DTunix<=TOLER ){
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
  
  void MuCSDT::endJob(){
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
