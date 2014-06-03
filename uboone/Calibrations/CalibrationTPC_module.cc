#ifndef CALIBRATIONTPC_H
#define CALIBRATIONTPC_H
/*!
 * Title:   CalibrationTPC class
 * Author:  wketchum@lanl.gov
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

  private:
    
    std::string       fRawDigitModuleLabel;   //label for rawdigit module
    unsigned int      fNFFTBins; //number of bins in noise FFT
    uint32_t          fNChanMax; //Max Channel Number used
    int               fNChannels; //Number of Channels from Geometry
    uint32_t          fsubRunNum; //Subrun Number taken from event
    // there will be a lot of other things here ...

    // these are containers for the calibration results
    // Intended design: each of these is reinitialized at subrun begin
    // Thus, pedestal_data[ie][ich] = mean pedestal for event ie, channel ich
    std::vector< std::vector<float> > fPedestalData;
    std::vector< std::vector<float> > fNoiseData;
    //Create a Map to link ordering of channel [ich] in above vector
    //to channel number.
    std::vector< std::map< unsigned int, uint32_t> > fChanMap;

    // noise_spectra[ie][ich][ifbin] 
    //   = fft spectrum amplitude for freq. bin ifbin, channel ich, event ie
    std::vector< std::vector< std::vector<float> > > fNoiseSpectra;

    //Output Histograms etc..
    TH1I* fAllChan;
    TH1D* fChannelBaseline;
    TH1D* fChannelNoise;
    

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
    fNChanMax            = 0;
  }


  //-------------------------------------------------------------------
  void CalibrationTPC::beginJob(){
    
    art::ServiceHandle<geo::Geometry> geo;
    fNChannels = 20000;//geo->Nchannels();
    fNChanMax  = fNChannels;

  }

  //-------------------------------------------------------------------
  void CalibrationTPC::endJob(){
  }

  //-------------------------------------------------------------------
  void CalibrationTPC::beginRun(art::Run const& run){
  }

  //-------------------------------------------------------------------
  void CalibrationTPC::endRun(art::Run const& run){
  }

  //-------------------------------------------------------------------
  void CalibrationTPC::beginSubRun(art::SubRun const& subrun){
  }

  //-------------------------------------------------------------------
  void CalibrationTPC::endSubRun(art::SubRun const& subrun){

    std::cout << "Subrun " << fsubRunNum << " just finished!" << std::endl;
    char subrunname[50];
    sprintf(subrunname,"subRun_%d",fsubRunNum);

    //Take containers with all info and insert that info in Histograms
    art::ServiceHandle<art::TFileService> tfs;


    char AllChanName[50];
    sprintf(AllChanName,"AllChan_%d",fsubRunNum);
    fAllChan            = tfs->make<TH1I>(AllChanName,           "Channels ;Chan Num; Num Evts", 
					  fNChanMax, 0, fNChanMax);
    char ChannelbaselineName[50];
    sprintf(ChannelbaselineName,"Channelbaseline_%d",fsubRunNum);
    fChannelBaseline    = tfs->make<TH1D>(ChannelbaselineName,   "Baseline By Channel ;Chan Num; Baseline [ADCs]",
					  fNChanMax, 0, fNChanMax);

    char ChannelNoiseName[50];
    sprintf(ChannelNoiseName,"ChannelNoise_%d",fsubRunNum);
    fChannelNoise       = tfs->make<TH1D>(ChannelNoiseName,      "RMS Noise by Channel ;Chan Num; RMS Noise [ADCs]",
					  fNChanMax, 0, fNChanMax);

    //loop over events and see how many events each channel has
    //then fill info in fAllChan
    //should be using fNChanMax instead of 1000
    std::vector<int> EvtsPerChan(fNChannels);
    std::vector<float> BaselinePerChan(fNChannels);
    std::vector<float> BaselineErrPerChan(fNChannels);
    std::vector<float> NoisePerChan(fNChannels);
    std::vector<float> NoiseErrPerChan(fNChannels);
    for ( uint32_t chanN=0; chanN < fNChanMax; chanN++){
      EvtsPerChan.at(chanN) = 0;
      BaselinePerChan.at(chanN) = 0;
      BaselineErrPerChan.at(chanN) = 0;
      NoisePerChan.at(chanN) = 0;
      NoiseErrPerChan.at(chanN) = 0;
    }

    int currentChan = 0;
    for ( size_t evtN=0; evtN < fChanMap.size(); evtN++ ){
      for ( size_t chanIndex=0; chanIndex < fChanMap[evtN].size(); chanIndex++ ){
	currentChan = fChanMap[evtN].at(chanIndex);
	EvtsPerChan.at(currentChan)     += 1;
	BaselinePerChan.at(currentChan) += fPedestalData.at(evtN).at(chanIndex);
	NoisePerChan.at(currentChan)    += fNoiseData.at(evtN).at(chanIndex);
      }
    }

    //average out baseline and noise values: divide by N events
    for ( uint32_t chanN=0; chanN < fNChanMax; chanN++){
      BaselinePerChan.at(chanN) /= EvtsPerChan[chanN];
      NoisePerChan.at(chanN) /= EvtsPerChan[chanN];
    }

    //Now with baseline and RMS noise per channel calculate
    //standard deviation of these values channel-by-channel
    for ( size_t evtN=0; evtN < fChanMap.size(); evtN++ ){
      for ( size_t chanIndex=0; chanIndex < fChanMap[evtN].size(); chanIndex++ ){
	currentChan = fChanMap[evtN].at(chanIndex);
	BaselineErrPerChan.at(currentChan) += (fPedestalData.at(evtN).at(chanIndex)-BaselinePerChan.at(currentChan))*
	  (fPedestalData.at(evtN).at(chanIndex)-BaselinePerChan.at(currentChan));
	NoiseErrPerChan.at(currentChan) += (fNoiseData.at(evtN).at(chanIndex)-NoisePerChan.at(currentChan))*
	  (fNoiseData.at(evtN).at(chanIndex)-NoisePerChan.at(currentChan));
      }
    }

    for ( uint32_t chanN=0; chanN < fNChanMax; chanN++){
      BaselineErrPerChan.at(chanN) = sqrt(BaselineErrPerChan.at(chanN));
      NoiseErrPerChan.at(chanN) = sqrt(NoiseErrPerChan.at(chanN));
    }

    for ( uint32_t chanN=0; chanN < fNChanMax; chanN++){
      if ( EvtsPerChan[chanN] > 0 ){
	fAllChan->SetBinContent( chanN+1, EvtsPerChan.at(chanN) );
	fChannelBaseline->SetBinContent( chanN+1, BaselinePerChan.at(chanN));
	fChannelBaseline->SetBinError( chanN+1, BaselineErrPerChan.at(chanN));
	fChannelNoise->SetBinContent( chanN+1, NoisePerChan.at(chanN));
	fChannelNoise->SetBinError( chanN+1, NoiseErrPerChan.at(chanN));
      }
    } 


  }

  
  //-------------------------------------------------------------------
  void CalibrationTPC::analyze(art::Event const& evt){

    //get subrun number to use when saving histograms
    fsubRunNum = evt.subRun();
    
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
    std::vector< std::vector<float> > 
      noise_spectrum(n_channels, 
		     std::vector<float>(fNFFTBins));
    std::map< unsigned int, uint32_t > chanmap;

    //now run the code
    analyzeEmptyEvent(rawDigitVector,
		      pedestals,
		      noise,
		      noise_spectrum);

    //make map for channel index to channel number
    genChanMap( rawDigitVector, chanmap, fNChanMax );

    fPedestalData.push_back( pedestals );
    fNoiseData.push_back( noise );
    fNoiseSpectra.push_back( noise_spectrum );
    fChanMap.push_back( chanmap );

    //    calcPedestal(rawDigitVector,pedestals);
  }



  DEFINE_ART_MODULE(CalibrationTPC)

} //end namespace calibration

#endif //CALIBRATIONTPC_H
