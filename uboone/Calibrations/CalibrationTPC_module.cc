#ifndef CALIBRATIONTPC_H
#define CALIBRATIONTPC_H
/*!
 * Title:   CalibrationTPC class
 * Authors:  wketchum@lanl.gov
 *           dcaratelli@nevis.columbia.edu
 *           ansmith01@email.wm.edu
 * Inputs:  raw::RawDigit
 * Outputs: Histograms and other nice data
 *
 * Description:
 * This analyzer is intended to look at RawData from calibration 
 * and calibration-like runs. It will include a number of calibration 
 * algorithms that create a number of monitoring histograms and produce results
 * of calibration tests.
 */

#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TTree.h>


#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Handle.h" 
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

#include "RawData/RawDigit.h"
#include "Geometry/Geometry.h"
#include "Utilities/DetectorProperties.h"

#include "CalibrationTPC_Algs.cc"

namespace calibration {

  class CalibrationTPC : public art::EDAnalyzer {

  public:
    explicit CalibrationTPC(fhicl::ParameterSet const& pset);
    virtual ~CalibrationTPC();

    void analyze(art::Event const& evt);
    void reconfigure(fhicl::ParameterSet const& pset);
    
    void beginJob();
    void endJob();

    //likely we will need begin/end run and subrun functions
    void beginRun(art::Run const& run);
    void endRun(art::Run const& run);
    void beginSubRun(art::SubRun const& subrun);
    void endSubRun(art::SubRun const& subrun);
    void MakeHisto( TH1D* & histo,
		    std::vector< std::vector<float> > Data);
    void GetWaveforms( std::vector<raw::RawDigit> const& rawDigit,
		       std::vector<TH1I*> & Waveforms,
		       std::map< unsigned int, uint32_t > chanmap );
    
  private:

    //Keep Track of event number
    int fEvtNum; //Number of current event
    
    std::string       fRawDigitModuleLabel;   //label for rawdigit module
    unsigned int      fNFFTBins;              //number of bins in noise FFT
    uint32_t          fNChanMax;              //Max Channel Number used
    int               fChanBegin;             //Chan Num from where to start
    int               fNChannels;             //Number of Channels from Geometry
    uint32_t          fsubRunNum;             //Subrun Number taken from event
    int               fEmptyRun;
    int               fGainRun;
    int               fShortRun;
    int               fprePulseTicks;         //N of ticks to use for baseline when pulse is present
    std::vector<float> fVmin;
    std::vector<float> fVmax;
    std::vector<float> fVstep;
    std::vector<float> fLinFitVmaxColl;
    std::vector<float> fLinFitVmaxInd;
    int               fLinBegin;
    std::vector<int>  fSubRunGain;
    std::vector<int>  fSubRunShape;
    // there will be a lot of other things here ...

    // these are containers for the calibration results
    // Intended design: each of these is reinitialized at subrun begin
    // Thus, fPedestalData[ie][ich] = mean pedestal for event ie, channel ich
    std::vector< std::vector<float> > fPedestalData;
    std::vector< std::vector<float> > fNoiseData;
    std::vector< std::vector<float> > fGainData;
    std::vector< std::vector<float> > fTimeData;
    //Create a Map to link ordering of channel [ich] in above vector
    //to channel number.
    std::vector< std::map< unsigned int, uint32_t> > fChanMap;
    std::vector< std::map< unsigned int, uint32_t> > fChanMapEndRun;
    

    //Vector for Noise Spectrum (a vector of floats) per event per channel
    std::vector< std::vector< std::vector<float> > > fNoiseSpectra;

    //Output Histograms etc..
    //----------------------
    //TTree for baseline and noise values
    TTree* fDistributions;
    //Variables for TTree
    Double_t fBaselineVal;
    Double_t fNoiseVal;
    //Histograms for Empty run - generally just one per run
    std::vector<TH1I*> fAllChan;
    std::vector<TH1D*> fChannelBaseline;
    std::vector<TH1D*> fChannelNoise;
    //Histograms one per run: one per gain setting, per rise-time, per input voltage
    std::vector<TH1D*> fChannelGain;
    std::vector<TH1D*> fChannelTimeMax;
    std::vector<TH1D*> fChannelGainVal;
    std::vector<TH1I*> fWaveforms; 
    //make histograms per channel, per gain, per shaping time
    std::vector<TH1D*>  fLinearity;
    std::vector<TF1*>   fFit;

  }; //end class CalibrationTPC


  //-------------------------------------------------------------------
  CalibrationTPC::CalibrationTPC(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset){ 
    this->reconfigure(pset); 
  }


  //-------------------------------------------------------------------
  CalibrationTPC::~CalibrationTPC(){}


  //-------------------------------------------------------------------
  void CalibrationTPC::reconfigure(fhicl::ParameterSet const& pset){
    fRawDigitModuleLabel = pset.get<std::string>("RawDigitModuleLabel");
    fNFFTBins            = pset.get<unsigned int>("NFFTBins");
    fNChanMax            = pset.get<int>("NChanMax");
    fGainRun             = pset.get<int>("GainRun");
    fShortRun            = pset.get<int>("ShortRun");
    fEmptyRun            = pset.get<int>("EmptyRun");
    fprePulseTicks       = pset.get<int>("prePulseTicks");
    fVmin                = pset.get<std::vector <float> >("Vmin"); 
    fVmax                = pset.get<std::vector <float> >("Vmax");
    fVstep               = pset.get<std::vector <float> >("Vstep");
    fLinFitVmaxColl      = pset.get<std::vector <float> >("LinFitVmaxColl");
    fLinFitVmaxInd       = pset.get<std::vector <float> >("LinFitVmaxInd");
    fLinBegin            = pset.get<int>("LinBegin");
    fSubRunGain          = pset.get< std::vector <int> >("SubRunGain");
    fSubRunShape         = pset.get< std::vector <int> >("SubRunShape");
  }


  //-------------------------------------------------------------------
  void CalibrationTPC::beginJob(){
    
    art::ServiceHandle<geo::Geometry> geo;
    fChanBegin = 0;
    //fNChannels = 20000;//geo->Nchannels();
    //fNChanMax  = fNChannels;
    fNChannels = fNChanMax;
    fsubRunNum = 0;

  }

  //-------------------------------------------------------------------
  void CalibrationTPC::endJob(){
  }

  //-------------------------------------------------------------------
  void CalibrationTPC::beginRun(art::Run const& run){
  }

  //-------------------------------------------------------------------
  void CalibrationTPC::endRun(art::Run const& run){

    art::ServiceHandle<art::TFileService> tfs;

    if ( fGainRun ){
      //Here take histos produced at the end of each subrun and produce
      //linearity-test plots etc...
      
      //fLinBegin gives number of subrun from start=0 where
      //linearity tests begin
      
      unsigned int lincounter = 0;

      for ( uint32_t sr=0; sr < fsubRunNum; sr++ ){ // loop through subruns
	//subruns where linearity start
	//2, 13, 24, ... (starts from 1)
	if ( lincounter == fSubRunShape.size() ) {
	  return;
	}
	else {
	  if ( ((sr-1)%11) == 0 ){ //means we have a new linearity start (Vin=0)
	    //Prepare Gain Histograms: fill with linear fit value (slope)
	    fChannelGainVal.push_back( tfs->make<TH1D>(Form("GainVal_Gain_%d_ShT_%d",fSubRunGain.at(lincounter), fSubRunShape.at(lincounter)),
						       "Value of Linear Fit",
						       fNChanMax+1, fChanBegin, fChanBegin+fNChanMax+1) );
	    //Number of steps in voltage:
	    int nsteps = 1+(int)((fVmax.at(lincounter)-fVmin.at(lincounter))/fVstep.at(lincounter));
	    for ( int ch = 0; ch < fNChannels; ch++){
	      //check and see if there are events for this chan number
	      if ( fAllChan.at(0)->GetBinContent(ch) > 0 ){ //means there are events
		fLinearity.push_back( tfs->make<TH1D>(Form("Linearity_Ch_%d_Gain_%d_ShT_%d", ch, fSubRunGain.at(lincounter), fSubRunShape.at(lincounter)),
						      Form("Linearity Test Ch %d Gain %d mV/fC Shaping Time %d usec; Pulse Voltage [mV]; Amplitude [ADCs]", ch, fSubRunGain.at(lincounter), fSubRunShape.at(lincounter) ), 
						      nsteps,
						      0-(fVstep.at(lincounter)/2.),
						      fVmax.at(lincounter)+(fVstep.at(lincounter)/2.) ) );
		//now loop through next 11 subruns and fill histogram
		for (int i = 0; i < nsteps; i++){
		  //Need to avoid points with Vin < 50 mV from being plotted
		  //That is because the pulser used on MRT sends huge pulses if input set to < 50 mV.
		  //This will distort linearity plot and fit therefore affecting the gain measurement
		  //Vin = i*fVstep.at(lincounter)
		  if ( (i == 0) or (i*fVstep.at(lincounter) >= 0.05) ) {
		    fLinearity.back()->SetBinContent( i+1, fChannelGain.at(sr+i)->GetBinContent(ch) );
		    fLinearity.back()->SetBinError( i+1,  fChannelGain.at(sr+i)->GetBinError(ch) / sqrt( fAllChan.at(sr+i)->GetBinContent(ch) ) );
		  }
		}//for all voltages
		
		//fit histogram
		//determine where max is (from LinFitVmax)
		//also depends on induction vs collection wire
		double fitmax;
		if ( (ch%100) < 32 )
		  fitmax = fLinFitVmaxInd.at(lincounter);
		else
		  fitmax = fLinFitVmaxColl.at(lincounter);
		  fFit.push_back( tfs->make<TF1>(Form("Linearity_Ch_%d_Gain_%d_ShT_%d", ch, fSubRunGain.at(lincounter), fSubRunShape.at(lincounter))
						 , "[0]+x*[1]", 0, fitmax ) );
		  fFit.back()->SetParameter(0,0);
		  fFit.back()->SetParameter(1,1);
		  fFit.back()->SetParName(0,"b");
		  fFit.back()->SetParName(1,"a");

		gStyle->SetOptStat(0);
		fLinearity.back()->Fit( fFit.back(), "R" );

		//Add Fit Slope to fChannelGainVal histogram
		fChannelGainVal.back()->SetBinContent( ch+1, fFit.back()->GetParameter(1) );
		
	      }//for channel is not empty
	    }//for all channels
	    //begin of new linearity test
	    lincounter += 1;
	  }//if the beginning of a new linearity test
	}
      }//for all subruns
      
    }//if gain run

    return;
  }


  //-------------------------------------------------------------------
  void CalibrationTPC::beginSubRun(art::SubRun const& subrun){

    art::ServiceHandle<art::TFileService> tfs;

    gStyle->SetOptStat(0);

    fWaveforms.clear();

    //First Make a channel map: number of events per channel
    fAllChan.push_back( tfs->make<TH1I>(Form("AllChan_SubRun_%d",fsubRunNum),
					"Channels ;Chan Num; Num Evts", 
					fNChanMax+1, fChanBegin, fChanBegin+fNChanMax+1) );
    

    //if subrunnumber indicates pulse-type subrun
    if ( fGainRun or fShortRun ){

      fChannelGain.push_back( tfs->make<TH1D>(Form("ChannelGain_%d",fsubRunNum),
					      "ADC Gain by channel ;Chan Num; Max [ADCs]",
					      fNChanMax+1, fChanBegin, fChanBegin+fNChanMax+1) );

      fChannelTimeMax.push_back( tfs->make<TH1D>(Form("ChannelTimeMax_%d",fsubRunNum),
					      "Time of Max ADC count ;Chan Num; Max Time [Tick Number]",
					      fNChanMax+1, fChanBegin, fChanBegin+fNChanMax+1) );

      fChannelBaseline.push_back( tfs->make<TH1D>(Form("Channelbaseline_%d",fsubRunNum),
						  "Baseline By Channel ;Chan Num; Baseline [ADCs]",
						  fNChanMax+1, fChanBegin, fChanBegin+fNChanMax+1) );

      fChannelNoise.push_back( tfs->make<TH1D>(Form("ChannelNoise_%d",fsubRunNum),
					       "RMS Noise by Channel ;Chan Num; RMS Noise [ADCs]",
					       fNChanMax+1, fChanBegin, fChanBegin+fNChanMax+1) );

    }

    if ( fEmptyRun ){

      fDistributions = tfs->make<TTree>(Form("Distributions_SubRun_%d",fsubRunNum),
					Form("Distributions_SubRun_%d",fsubRunNum) );

      fChannelBaseline.push_back( tfs->make<TH1D>(Form("Channelbaseline_%d",fsubRunNum),
						  "Baseline By Channel ;Chan Num; Baseline [ADCs]",
						  fNChanMax+1, fChanBegin, fChanBegin+fNChanMax+1) );

      fDistributions->Branch("BaselineVal", &fBaselineVal, "BaselineVal/D");

      fChannelNoise.push_back( tfs->make<TH1D>(Form("ChannelNoise_%d",fsubRunNum),
					       "RMS Noise by Channel ;Chan Num; RMS Noise [ADCs]",
					       fNChanMax+1, fChanBegin, fChanBegin+fNChanMax+1) );

      fDistributions->Branch("NoiseVal", &fNoiseVal, "BaselineVal/D");
            
    }
     
    fsubRunNum += 1;

  }

  //-------------------------------------------------------------------
  void CalibrationTPC::endSubRun(art::SubRun const& subrun){

    //Take containers with all info and insert that info in Histograms
    art::ServiceHandle<art::TFileService> tfs;

    std::vector<int> EvtsPerChan(fNChanMax, 0);
    
    int currentChan = 0;
    for ( size_t evtN=0; evtN < fChanMap.size(); evtN++ ){
      // std::cout << "Event Number: " << evtN << " Number of Channels: " << fChanMap[evtN].size() << std::endl;
      for ( size_t chanIndex=0; chanIndex < fChanMap[evtN].size(); chanIndex++ ){
	currentChan = fChanMap[evtN].at(chanIndex);
	//	std::cout << "Current Channel: " << currentChan << std::endl;
	EvtsPerChan.at(currentChan)     += 1;
      }
    }

    for ( uint32_t chanN=0; chanN < fNChanMax; chanN++){
      if ( EvtsPerChan[chanN] > 0 ){
	fAllChan.at(fsubRunNum-1)->SetBinContent( chanN+1, EvtsPerChan.at(chanN) );
      }
    }


    //Determine type of subrun and fill corresponding histogram
    if ( fEmptyRun ){
      MakeHisto( fChannelBaseline.at(fsubRunNum-1), fPedestalData);
      MakeHisto( fChannelNoise.at(fsubRunNum-1), fNoiseData);

      //Loop through vectors and fill TTree:
      for ( size_t ev = 0; ev < fPedestalData.size(); ev++ ){
	for ( size_t ch = 0; ch < fPedestalData.at(ev).size(); ch++) {
	  fBaselineVal = fPedestalData.at(ev).at(ch);
	  fNoiseVal = fNoiseData.at(ev).at(ch);
	  fDistributions->Fill();
	}
      }

    }//end emptyrun
    
    if ( fGainRun or fShortRun ){
      MakeHisto( fChannelGain.at(fsubRunNum-1), fGainData);
      MakeHisto( fChannelBaseline.at(fsubRunNum-1), fPedestalData);
      MakeHisto( fChannelNoise.at(fsubRunNum-1), fNoiseData);
      MakeHisto( fChannelTimeMax.at(fsubRunNum-1), fTimeData);
    }//end gainrun and shortrun

    //Clear per-subrun vectors
    fChanMapEndRun = fChanMap;

    fPedestalData.clear();
    fNoiseData.clear();
    fNoiseSpectra.clear();
    fChanMap.clear();
    fGainData.clear();

  }


  //-------------------------------------------------------------------
  void CalibrationTPC::MakeHisto( TH1D* & histo,
				  std::vector< std::vector<float> > Data){

    std::vector<int> EvtsPerChan(fNChanMax, 0);
    std::vector<float> DataPerChan(fNChanMax, 0);
    std::vector<float> DataErrPerChan(fNChanMax, 0);

    gStyle->SetOptStat(0);

    int currentChan = 0;

    //fill data and event number vectors
    for ( size_t evtN=0; evtN < fChanMap.size(); evtN++ ){
      for ( size_t chanIndex=0; chanIndex < fChanMap[evtN].size(); chanIndex++ ){
	currentChan = fChanMap[evtN].at(chanIndex);
	EvtsPerChan.at(currentChan)     += 1;
	DataPerChan.at(currentChan) += Data.at(evtN).at(chanIndex);
      }
    }

    //average out baseline and noise values: divide by N events
    for ( uint32_t chanN=0; chanN < fNChanMax; chanN++)
      if ( EvtsPerChan[chanN] != 0 )
	DataPerChan.at(chanN) /= EvtsPerChan[chanN];
      else
	DataPerChan.at(chanN) = 0;

    //Now with Data per channel calculate
    //standard deviation of these values channel-by-channel
    for ( size_t evtN=0; evtN < fChanMap.size(); evtN++ ){
      for ( size_t chanIndex=0; chanIndex < fChanMap[evtN].size(); chanIndex++ ){
	currentChan = fChanMap[evtN].at(chanIndex);
	DataErrPerChan.at(currentChan) += (Data.at(evtN).at(chanIndex)-DataPerChan.at(currentChan))*
	  (Data.at(evtN).at(chanIndex)-DataPerChan.at(currentChan));
      }
    }
    
    //calculate spread in values
    for ( uint32_t chanN=0; chanN < fNChanMax; chanN++)
      DataErrPerChan.at(chanN) = sqrt(DataErrPerChan.at(chanN));

    //fill histogram chan-by-chan
    for ( uint32_t chanN=0; chanN < fNChanMax; chanN++){
      if ( EvtsPerChan[chanN] > 0 ){
	histo->SetBinContent( chanN+1, DataPerChan.at(chanN));
	histo->SetBinError( chanN+1, DataErrPerChan.at(chanN));
      }
    } 

    gStyle->SetOptStat(0);

  }

 
  //-------------------------------------------------------------------------
  void CalibrationTPC::GetWaveforms( std::vector<raw::RawDigit> const& rawDigit,
				     std::vector<TH1I*> & Waveforms,
				     std::map< unsigned int, uint32_t > chanmap ){

    art::ServiceHandle<art::TFileService> tfs;
    //Make a directory specific for waveforms
    //tfs->mkdir(Form("WaveForms_sr_%d_ev_%d",fsubRunNum, fEvtNum) ,"Histograms of ADC Waveforms");

    for (unsigned int ich=0; ich < chanmap.size(); ich++){

      Waveforms.push_back( tfs->make<TH1I>(Form("WF_subRun_%d_ev_%d_chan_%d", fsubRunNum, fEvtNum, chanmap.at(ich) ), "Waveform", 
					   rawDigit.at(ich).fADC.size(),
					   0, rawDigit.at(ich).fADC.size() ) );

      for ( unsigned int tick=0; tick < rawDigit.at(ich).fADC.size(); tick++){
	Waveforms.back()->SetBinContent( tick+1, rawDigit.at(ich).fADC.at(tick) );
      }
    }//for all channels

    return;
    
  }
  


  
  //-------------------------------------------------------------------
  void CalibrationTPC::analyze(art::Event const& evt){

    fEvtNum = evt.event();

    unsigned int run = evt.run();
    unsigned int subrun = evt.subRun();
    
    LOG_INFO ("CalibrationTPC")
      << "Processing Run " << run 
      << ", Subrun " << subrun 
      << ", Event " << evt.event();
    
    art::Handle< std::vector<raw::RawDigit> > rawDigitHandle;
    evt.getByLabel(fRawDigitModuleLabel,rawDigitHandle);
    std::vector<raw::RawDigit> const& rawDigitVector(*rawDigitHandle);


    //initialize per-event vectors
    const size_t n_channels = rawDigitVector.size();
    std::vector<float> pedestals(n_channels);
    std::vector<float> noise(n_channels);
    std::vector< std::vector<float> > noise_spectrum(n_channels, std::vector<float>(fNFFTBins));
    std::vector<float> maxADCs(n_channels);
    std::vector<float> minADCs(n_channels);
    std::vector<float> timeMax(n_channels);
    std::map< unsigned int, uint32_t > chanmap;
    

    //decide if empty or pulse subrun
    if (fEmptyRun){
    
      //now run the code
      analyzeEmptyEvent(rawDigitVector,
			pedestals,
			noise,
			noise_spectrum);
      
      //make map for channel index to channel number
      genChanMap( rawDigitVector, chanmap );
      
      fPedestalData.push_back( pedestals );
      fNoiseData.push_back( noise );
      fNoiseSpectra.push_back( noise_spectrum );
      fChanMap.push_back( chanmap );
    }//empty subrun

    if (fGainRun){ //pulse data

      //now run the code
      analyzeGainEvent(rawDigitVector,
		       pedestals,
		       noise,
		       maxADCs,
		       minADCs,
		       timeMax,
		       fprePulseTicks);
      
      //make map for channel index to channel number
      genChanMap( rawDigitVector, chanmap );

      fPedestalData.push_back( pedestals );
      fNoiseData.push_back( noise );
      fGainData.push_back( maxADCs );
      fChanMap.push_back( chanmap );
      fTimeData.push_back( timeMax );
    }//gain subrun


    //if short run: generally for single-subrun run
    //plot waveform
    if (fShortRun){

      //Save waveforms to histogram
      if ( evt.event() == 1 )
	GetWaveforms( rawDigitVector, fWaveforms, chanmap );

    }//short subrun

  }



  DEFINE_ART_MODULE(CalibrationTPC)

} //end namespace calibration

#endif //CALIBRATIONTPC_H
