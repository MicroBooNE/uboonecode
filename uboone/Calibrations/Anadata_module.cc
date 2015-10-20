////////////////////////////////////////////////////////////////////////
//
// module to create a TTree for TPCnoise analysis
//
//
////////////////////////////////////////////////////////////////////////
#ifndef ANADATA_H
#define ANADATA_H

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/FindMany.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
<<<<<<< HEAD
#include "SimpleTypesAndConstants/RawTypes.h"                       
=======
#include "SimpleTypesAndConstants/RawTypes.h"                                                                                             
>>>>>>> aaed7c028546a8718d7a68dd62af7449d1723ba2
#include "Filters/ChannelFilter.h"
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "RawData/BeamInfo.h"
#include "RawData/DAQHeader.h"
#include "RawData/OpDetWaveform.h"
#include "RawData/OpDetPulse.h"
#include "Geometry/Geometry.h"
<<<<<<< HEAD
#include "Utilities/LArFFT.h"

#include <iostream>
#include <cstring> 
#include <vector>
#include <map>
#include <iterator>
=======

#include <cstring> // std::memcpy()
#include <vector>
#include <map>
#include <iterator> // std::begin(), std::end()
>>>>>>> aaed7c028546a8718d7a68dd62af7449d1723ba2
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
<<<<<<< HEAD
#include <functional>
=======
#include <functional> // std::mem_fun_ref
>>>>>>> aaed7c028546a8718d7a68dd62af7449d1723ba2
#include <typeinfo>

#include "TTree.h"
#include "TTimeStamp.h"
#include "TH2.h"
#include "TFile.h"
<<<<<<< HEAD
#include "TString.h"


const int kMaxPrimaries = 20000; 
=======


const int kMaxPrimaries = 20000; //maximum number of primary particles
>>>>>>> aaed7c028546a8718d7a68dd62af7449d1723ba2
namespace microboone {
  
  class Anadata : public art::EDAnalyzer {
    
  public:
    
    explicit Anadata(fhicl::ParameterSet const& pset); 
    virtual ~Anadata();
    
    void analyze(const art::Event& evt);
    void beginJob();
    void beginSubRun(const art::SubRun& sr);
<<<<<<< HEAD
    void cal_ped_rms( std::vector<short> const& rawData,float& pedestal,
		      float&                    noise,
		      float&                    maxADC,
		      float&                    minADC,
		      int&                      maxADCbin,
		      int&                      minADCbin,
		      int const&                startBin,
		      int const&                endBin);
    
=======
>>>>>>> aaed7c028546a8718d7a68dd62af7449d1723ba2
    
  private:
    
    
    void   ResetVars();
    
    TTree* fTree;
    //run information
    Int_t    run;                  
    Int_t    subrun;               
    Int_t    event;                
    Double_t evttime;              
    Double_t beamtime;             
    Char_t   isdata;//flag, 0=MC 1=data
<<<<<<< HEAD
=======
    uint32_t channel[kMaxPrimaries];
>>>>>>> aaed7c028546a8718d7a68dd62af7449d1723ba2
    Double_t adc_rms[kMaxPrimaries];
    Double_t adc_rms_wire_plane_0[kMaxPrimaries];
    Double_t adc_rms_wire_plane_1[kMaxPrimaries];
    Double_t adc_rms_wire_plane_2[kMaxPrimaries];
    Double_t pedstal_wire_plane_0[kMaxPrimaries];
    Double_t pedstal_wire_plane_1[kMaxPrimaries];
    Double_t pedstal_wire_plane_2[kMaxPrimaries];
<<<<<<< HEAD
    Int_t tot_channel;
    uint32_t channel[kMaxPrimaries];
    Float_t val_adc[kMaxPrimaries];
    Int_t    tot_waveform;
    Float_t pedstal_forchannel[kMaxPrimaries];
    Int_t fMaxTicks; Int_t fchannel;
    Float_t fMax, fMin, fMean, fRms;
    Int_t fMaxtime, fMintime;
    
    
    TH1F *fHistMod;
    TH1F** waveNoiseHists;
    TH1F *currentHist;
    TH1F *currentFFTHist;
    TH1F *ped_histo;
    TH1F *rms_histo;
    TH2F *ped_histo_2d;
    TH2F *rms_histo_2d;
    TH2F *ped_histo_nevt;
    TH2F *rms_histo_nevt;
    TH2F *rms_histo_wire_plane_0;
    TH2F *rms_histo_wire_plane_1;
    TH2F *rms_histo_wire_plane_2;
    
    Int_t N_evt_processed = 1;
    TH1F *h_adc_Time_[8256];
    TH1F *h_amp_Freq_[8256];
    TH2F *h_plane_Wire_[8256];
=======
    Int_t    tot_channel;
    Int_t    tot_waveform;
    Double_t pedstal_forchannel[kMaxPrimaries];
    Int_t fDataSize;
    std::vector<double>holder;
    std::vector<double>diff_val_holder;
    Int_t event_n;
    std::vector<short> rawadc;
    std::vector<double> channels_pedestal; 
    
    TH2D *ped_histo_nevt;
    TH2D *rms_histo_nevt;
    Int_t N_evt_processed = 10;
    
    
    
    
>>>>>>> aaed7c028546a8718d7a68dd62af7449d1723ba2
    std::string processname[kMaxPrimaries];
    std::string fRawDigitModuleLabel;
    std::string fRawOpDetWaveModuleLabel;
    
<<<<<<< HEAD
    
  };
}//name space
microboone::Anadata::Anadata(fhicl::ParameterSet const& pset) : 
  EDAnalyzer(pset), fMaxTicks(9595), fchannel(8256),
=======
  };
}//name space

microboone::Anadata::Anadata(fhicl::ParameterSet const& pset) : 
  EDAnalyzer(pset),
>>>>>>> aaed7c028546a8718d7a68dd62af7449d1723ba2
  fRawDigitModuleLabel (pset.get<std::string>("RawDigitModuleLabel"))
{
}

//-------------------------------------------------
microboone::Anadata::~Anadata()
{
}

void microboone::Anadata::beginJob(){
  
  art::ServiceHandle<art::TFileService> tfs;
  fTree= tfs->make<TTree>("noiseonrawtree","noiseonraw");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("evttime",&evttime,"evttime/D");
  fTree->Branch("beamtime",&beamtime,"beamtime/D");
  fTree->Branch("isdata",&isdata,"isdata/B"); 
  fTree->Branch("tot_channel",&tot_channel,"tot_channel/I");
<<<<<<< HEAD
  fTree->Branch("channel", channel, "channel[tot_channel]/I");
  fTree->Branch("pedstal_forchannel", pedstal_forchannel, "pedstal_forchannel[tot_channel]/D");
  fTree->Branch("adc_rms", adc_rms, "adc_rms[tot_channel]/D");
=======
  fTree->Branch("pedstal_forchannel", pedstal_forchannel, "pedstal_forchannel[tot_channel]/D");
  fTree->Branch("channel", channel, "channel[tot_channel]/I");
  fTree->Branch("adc_rms", adc_rms, "adc_rms[tot_channel]/D");
  fTree->Branch("event_n", &event_n,  "event_n/I");
>>>>>>> aaed7c028546a8718d7a68dd62af7449d1723ba2
  fTree->Branch("adc_rms_wire_plane_0", adc_rms_wire_plane_0, "adc_rms_wire_plane_0[tot_channel]/D");
  fTree->Branch("adc_rms_wire_plane_1", adc_rms_wire_plane_1, "adc_rms_wire_plane_1[tot_channel]/D");
  fTree->Branch("adc_rms_wire_plane_2", adc_rms_wire_plane_2, "adc_rms_wire_plane_2[tot_channel]/D");
  fTree->Branch("pedstal_wire_plane_0", pedstal_wire_plane_0, "pedstal_wire_plane_0[tot_channel]/D");
  fTree->Branch("pedstal_wire_plane_1", pedstal_wire_plane_1, "pedstal_wire_plane_1[tot_channel]/D");
  fTree->Branch("pedstal_wire_plane_2", pedstal_wire_plane_2, "pedstal_wire_plane_2[tot_channel]/D");
<<<<<<< HEAD
  
  ped_histo = tfs->make<TH1F> ("ped_histo", "pedestal_forchannel", fchannel, 0 ,fchannel);
  rms_histo = tfs->make<TH1F> ("rms_histo", "rms_forchannel", fchannel, 0 ,fchannel);
  ped_histo_2d = tfs->make<TH2F> ("ped_histo_2d", "pedestal_forchannel_2d", fchannel, 0 ,fchannel,100,0,2500);
  rms_histo_2d = tfs->make<TH2F> ("rms_histo_2d", "rms_forchannel_2d", fchannel, 0 ,fchannel,100,0,50);

  ped_histo_nevt = tfs->make<TH2F> ("ped_histo_nevt", "pedestal_forchannel_nevents_2d", fchannel, 0 ,fchannel,100,0,2500);
  rms_histo_nevt = tfs->make<TH2F> ("rms_histo_nevt", "rms_forchannel_nevents_2d", fchannel, 0 ,fchannel, 100,0,100);
  
  rms_histo_wire_plane_0 = tfs->make<TH2F>("rms_histo_wire_plane_0","rms_forchannel_wireplane_0",2400, 0, 2400,100, 0, 100);
  rms_histo_wire_plane_1 = tfs->make<TH2F>("rms_histo_wire_plane_1","rms_forchannel_wireplane_1",2400, 0, 2400,100, 0, 100);
  rms_histo_wire_plane_2 = tfs->make<TH2F>("rms_histo_wire_plane_2","rms_forchannel_wireplane_2",3456, 0,3456,100, 0, 100);
  
  
  for(int i=0 ; i<fchannel; i++)
    { std::ostringstream  start;
      start<<i;
      TString start1 = start.str();
      TString numstr = start1;
      TString histoname1 = "h_adc_time_"+numstr;
      TString histoname2 = "h_amp_freq_"+numstr;
      TString histoname3 = "h_plane_wire_"+numstr;
      h_adc_Time_[i]     = tfs->make <TH1F>(histoname1,"h_adc_time_plot",fMaxTicks,-0.5,fMaxTicks-0.5);
      h_amp_Freq_[i]     = tfs->make <TH1F>(histoname2, "h_FFT_plot" , fMaxTicks,-0.5,fMaxTicks-0.5);
      h_plane_Wire_[i]   = tfs->make <TH2F>(histoname3, "plane and wire no", 3,0,3,4000,0,4000); 
    }
  
  
  fHistMod         = new TH1F("fHistMod","",fMaxTicks,-0.5,fMaxTicks-0.5);
  currentHist      = new TH1F("currentHist","",fMaxTicks,-0.5,fMaxTicks-0.5);
  currentFFTHist   = new TH1F("currentFFTHist","",fMaxTicks,-0.5,fMaxTicks-0.5);
  waveNoiseHists   = new TH1F * [9595];
  for(unsigned i = 0; i < 9595; i++){
    waveNoiseHists[i] = new TH1F(Form("waveNoiseHist_number%d",i),";Tick;ADC Value",fMaxTicks,-0.5,fMaxTicks-0.5);}
=======

  ped_histo_nevt = tfs->make<TH2D> ("ped_histo_nevt", "pedestal_forchannel_nevents", 8256, 0 ,8256, 100,0,2500);
  rms_histo_nevt = tfs->make<TH2D> ("rms_histo_nevt", "rms_forchannel_nevents", 8256, 0 ,8256, 100,0,50);
>>>>>>> aaed7c028546a8718d7a68dd62af7449d1723ba2
  
}

void microboone::Anadata::beginSubRun(const art::SubRun& sr)
{
}

void microboone::Anadata::analyze(const art::Event& evt)
{  
  ResetVars();
<<<<<<< HEAD
=======
  //services
>>>>>>> aaed7c028546a8718d7a68dd62af7449d1723ba2
  
  run = evt.run();
  subrun = evt.subRun();
  event = evt.event();
<<<<<<< HEAD
  
=======

>>>>>>> aaed7c028546a8718d7a68dd62af7449d1723ba2
  art::Timestamp ts = evt.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  evttime = tts.AsDouble();
  
<<<<<<< HEAD
=======
  //copied from MergeDataPaddles.cxx
>>>>>>> aaed7c028546a8718d7a68dd62af7449d1723ba2
  art::Handle< raw::BeamInfo > beam;
  if (evt.getByLabel("beam",beam)){
    beamtime = (double)beam->get_t_ms(); /*millisecond*/
    beamtime/=1000.; //in second
  }
<<<<<<< HEAD
  
=======

>>>>>>> aaed7c028546a8718d7a68dd62af7449d1723ba2
  art::ServiceHandle<geo::Geometry> fGeometry;
  unsigned int maxChannels  = fGeometry->Nchannels();
  unsigned int wireMaxNum[] = {fGeometry->Nwires(0),fGeometry->Nwires(1),fGeometry->Nwires(2)};
  /* 3456 Y wires arrayed vertically and the 2400 U and 2400 V*/
<<<<<<< HEAD
  std::cout<< " /maxChannels:"<< maxChannels << " /wireMaxNum[0]:"<< wireMaxNum[0]<<
    " /wireMaxNum[1]:"<< wireMaxNum[1]<< " /wireMaxNum[2]:"<< wireMaxNum[2]<<std::endl;
=======
  std::cout<< " /maxChannels:"<< maxChannels << " /wireMaxNum[0]:"<< wireMaxNum[0]<< " /wireMaxNum[1]:"<< wireMaxNum[1]<< " /wireMaxNum[2]:"<< wireMaxNum[2]<<std::endl;
  
>>>>>>> aaed7c028546a8718d7a68dd62af7449d1723ba2
  
  if (evt.isRealData()){
    isdata = 1;
  }
  else isdata = 0;
  
<<<<<<< HEAD
  //======================Start======================  
  
  art::Handle< std::vector<raw::RawDigit> > rawDigitHandle;
  evt.getByLabel(fRawDigitModuleLabel,rawDigitHandle);
  tot_channel =rawDigitHandle->size();
  std::vector< short > waveform;
  std::vector<raw::RawDigit> const& rawDigitVector(*rawDigitHandle);
  
  for(int ich = 0; ich < tot_channel; ich++)
    {
      fHistMod->Reset();
      currentFFTHist->Reset();
      waveform.clear();
      const size_t n_samp = rawDigitVector.at(ich).NADC();
      channel[ich] = rawDigitVector.at(ich).Channel();
      
      for(unsigned int i_samp = 0 ; i_samp < n_samp ; i_samp++ )
	{
	  fHistMod->SetBinContent(i_samp+1, rawDigitVector.at(ich).ADCs().at(i_samp) );
	  val_adc[i_samp] = fHistMod->GetBinContent(i_samp+1);
	}
      
      for( int s = 0 ; s < fHistMod->GetNbinsX() ; s++ )
	{
	  waveNoiseHists[s]->SetBinContent(s+1, fHistMod->GetBinContent(s+1) );
	  waveform.push_back(waveNoiseHists[s]->GetBinContent( s+1 ));
	}	
      
      cal_ped_rms( waveform, fMean , fRms, fMax, fMin, fMaxtime, fMintime, 0, (int) n_samp);
      
      pedstal_forchannel[ich] = (double)fMean;
      adc_rms[ich] = (double)fRms;

      
      
      std::vector<geo::WireID> wire_ids;
      wire_ids = fGeometry->ChannelToWire(rawDigitVector.at(ich).Channel());
      int plane = wire_ids[0].Plane;
      int wire = wire_ids[0].Wire;
      
      h_plane_Wire_[ich]->Fill(plane, wire);
      TH1 *hm = 0;
      hm = fHistMod->FFT(fHistMod, "MAG");
      Int_t n_Bins = fHistMod->GetNbinsX();
      for(Int_t i = 0; i<n_Bins; i++)
        {      
	  h_adc_Time_[ich]->Fill(i,val_adc[i]);
	  h_amp_Freq_[ich]->Fill(i+1,hm->GetBinContent(i+2));
	  
        }
      
      ped_histo->Fill(rawDigitVector.at(ich).Channel(),fMean) ;
      rms_histo->Fill(rawDigitVector.at(ich).Channel(),fRms) ;
      ped_histo_2d->Fill(rawDigitVector.at(ich).Channel(),fMean) ;
      rms_histo_2d->Fill(rawDigitVector.at(ich).Channel(),fRms) ;
      ped_histo_nevt->Fill(rawDigitVector.at(ich).Channel(),fMean, 1/(float)N_evt_processed) ;
      rms_histo_nevt->Fill(rawDigitVector.at(ich).Channel(),fRms, 1/(float)N_evt_processed) ;
      
      if(wire_ids[0].Plane==0) {adc_rms_wire_plane_0[ich]=fRms; pedstal_wire_plane_0[ich]=fMean; rms_histo_wire_plane_0->Fill(wire_ids[0].Wire,fRms);}
      if(wire_ids[0].Plane==1) {adc_rms_wire_plane_1[ich]=fRms;  pedstal_wire_plane_1[ich]=fMean; rms_histo_wire_plane_1->Fill(wire_ids[0].Wire,fRms);}
      if(wire_ids[0].Plane==2) {adc_rms_wire_plane_2[ich]=fRms;  pedstal_wire_plane_2[ich]=fMean; rms_histo_wire_plane_2->Fill(wire_ids[0].Wire,fRms);}
      
      
    }//tot_channel
  
  
  
=======
  
  
  art::Handle< std::vector<raw::RawDigit> > rawDigitHandle;
  evt.getByLabel(fRawDigitModuleLabel,rawDigitHandle);
  tot_channel = rawDigitHandle->size();
  
  
  event_n = 0;
  for(size_t it = 0 ; it< rawDigitHandle->size(); ++it)
    {
      
      holder.clear();
      diff_val_holder.clear();
      rawadc.clear();
      
      art::Ptr<raw::RawDigit> digitVec(rawDigitHandle, it);
      fDataSize= digitVec->Samples();
      holder.resize(fDataSize);
      diff_val_holder.resize(fDataSize);
      rawadc.resize(fDataSize);
      
      double adc_pedestal = 0.;
      double adc_stdev    = 0.;
      double adc_diff_sq  = 0.;
      
      raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
      
      for(int bin = 0; bin < fDataSize; ++bin){
	holder[bin] = rawadc[bin];   // Rawhitfinder_module.
	adc_pedestal += holder[bin]/fDataSize;
	
      }
      
      pedstal_forchannel[it] = adc_pedestal;
      channel[it]            = digitVec->Channel();
      
      ped_histo_nevt->Fill(digitVec->Channel(),adc_pedestal,1/(float)N_evt_processed);
      for (int bin = 0 ; bin < fDataSize; ++bin)
	{ 
	  diff_val_holder[bin] = (rawadc[bin]-adc_pedestal);
	  adc_diff_sq += (diff_val_holder[bin]*diff_val_holder[bin])/fDataSize;
	}
      adc_stdev = sqrt(adc_diff_sq);
      adc_rms[it] = adc_stdev;
      rms_histo_nevt->Fill(digitVec->Channel(),adc_stdev,1/(float)N_evt_processed);
      
      //get the wire plane... stuffs...
      std::vector<geo::WireID> wire_ids;
      wire_ids = fGeometry->ChannelToWire(digitVec->Channel());
      
      /*three wire planes*/ 
      /*std::cout<< "/wire_ID: "<< wire_ids[0].Plane<< " /wire:"<<wire_ids[0].Wire<<  "/pedstal_forchannel[it]:"<<pedstal_forchannel[it]<<std::endl;*/
      
      if(wire_ids[0].Plane==0) {adc_rms_wire_plane_0[it]=adc_stdev; pedstal_wire_plane_0[it]=adc_pedestal;}
      
      if(wire_ids[0].Plane==1) {adc_rms_wire_plane_1[it]=adc_stdev;  pedstal_wire_plane_1[it]=adc_pedestal;}
      
      if(wire_ids[0].Plane==2) {adc_rms_wire_plane_2[it]=adc_stdev;  pedstal_wire_plane_2[it]=adc_pedestal;}    
      
    }
  
  
  event_n++;
>>>>>>> aaed7c028546a8718d7a68dd62af7449d1723ba2
  fTree->Fill();
  
}

<<<<<<< HEAD
void microboone::Anadata::cal_ped_rms( std::vector<short> const& rawData, float & pedestal, float & noise, float & maxADC,
				       float & minADC, int & maxADCbin, int & minADCbin, int const& startBin, int const& endBin)
{
  int n_samples = (int) rawData.size();
  
  if(startBin >= endBin || startBin < 0 || endBin > n_samples) return;
  pedestal = 0;
  noise = 0; 
  maxADC = 0; 
  minADC = 4095;
  maxADCbin =0;
  minADCbin =0;
  int count = 0;
  for(unsigned int it= (unsigned int)startBin; it< (unsigned int)endBin; it++){
    short thisADC = rawData.at(it);
    count++;
    pedestal += thisADC;
    if ( thisADC > maxADC ){
      maxADC = thisADC;
      maxADCbin = it;
    }
    if ( thisADC < minADC ){
      minADC = thisADC;
      minADCbin = it;
    }
  }
  if( count > 0 )
    pedestal /= (float)count;
  
  count = 0;
  for(unsigned int it= (unsigned int) startBin; it< (unsigned int) endBin; it++){
    count++;
    short thisADC = rawData.at(it);
    noise += (thisADC-pedestal)*(thisADC-pedestal);
  }
  if( count -1 > 0 )
    noise = sqrt( noise / (double)(count -1) );
  
  return;
}






=======
>>>>>>> aaed7c028546a8718d7a68dd62af7449d1723ba2
void microboone::Anadata::ResetVars(){
  
  run      = -99999;
  subrun   = -99999;
  event    = -99999;
  evttime  = -99999;
  beamtime = -99999;
  isdata   = -99;
<<<<<<< HEAD
  //event_n = -999999;
=======
  event_n = -999999;
>>>>>>> aaed7c028546a8718d7a68dd62af7449d1723ba2
  
  for (int i = 0; i<kMaxPrimaries; ++i)
    {
      pedstal_forchannel[i] = -999999;
      channel[i]            = -999999;
      adc_rms[i]            = -999999.;
      adc_rms_wire_plane_0[i]=-100.;
      adc_rms_wire_plane_1[i]=-100.;
      adc_rms_wire_plane_2[i]=-100.;
      pedstal_wire_plane_0[i]=-100.;
      pedstal_wire_plane_1[i]=-100.;
      pedstal_wire_plane_2[i]=-100.;
<<<<<<< HEAD
      val_adc[i] = -100.;
=======
>>>>>>> aaed7c028546a8718d7a68dd62af7449d1723ba2
    }
}
namespace microboone{
  
  DEFINE_ART_MODULE(Anadata)
  
} 

#endif



