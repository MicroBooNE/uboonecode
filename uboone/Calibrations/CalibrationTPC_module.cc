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

#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Handle.h" 

#include "RawData/RawDigit.h"

#include "CalibrationTPC_Algs.h"

namespace calibration {

  class CalibrationTPC : public art::EDAnalyzer {

  public:
    explicit CalibrationTPC(fhicl::ParameterSet const& pset);
    virtual ~CalibrationTPC();

    void analyze(art::Event& evt);
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
    // there will be a lot of other things here ...

    // these are containers for the calibration results
    // Intended design: each of these is reinitialized at subrun begin
    // Thus, pedestal_data[ie][ich] = mean pedestal for event ie, channel ich
    std::vector< std::vector<float> > fPedestalData;
    std::vector< std::vector<float> > fNoiseData;

    // noise_spectra[ie][ich][ifbin] 
    //   = fft spectrum amplitude for freq. bin ifbin, channel ich, event ie
    std::vector< std::vector< std::vector<float> > > fNoiseSpectra;

  }; //end class CalibrationTPC


  //-------------------------------------------------------------------
  CalibrationTPC::CalibrationTPC(fhicl::ParameterSet const& pset){ 
    this->reconfigure(pset); 
  }


  //-------------------------------------------------------------------
  CalibrationTPC::~CalibrationTPC(){}


  //-------------------------------------------------------------------
  void CalibrationTPC::reconfigure(fhicl::ParameterSet const& pset){
    fRawDigitModuleLabel = p.get<std::string>("RawDigitModuleLabel");
    fNFFTBins            = p.get<unsigned int>("NFFTBins");
  }


  //-------------------------------------------------------------------
  void CalibrationTPC::beginJob(){
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
  }

  
  //-------------------------------------------------------------------
  void CalibrationTPC::analyze(art::Event& evt){

    unsigned int run = evt.run();
    unsigned int subrun = evt.subRun();

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


    //now run the code
    analyzeEmptyEvent(rawDigitVector,
		      pedestals,
		      noise,
		      noise_spectrum);
  }



  DEFINE_ART_MODULE(CalibrationTPC)

} //end namespace calibration

#endif //CALIBRATIONTPC_H
