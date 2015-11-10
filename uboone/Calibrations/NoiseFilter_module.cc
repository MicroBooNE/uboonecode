#ifndef NOISEFILTERMODULE_H
#define NOISEFILTERMODULE_H

#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TTree.h>
#include <THStack.h>
#include <TGraph.h>
#include <TMath.h>
#include <TComplex.h>
#include "TApplication.h"

#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Handle.h" 
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

#include "Utilities/LArFFT.h"
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "Geometry/Geometry.h"
// #include "Utilities/DetectorPropertiesService.h"
#include "CalibrationDBI/Interface/IDetPedestalService.h"
#include "CalibrationDBI/Interface/IDetPedestalProvider.h"

namespace calibration {

  class NoiseFilter : public art::EDProducer {

  public:
    explicit NoiseFilter(fhicl::ParameterSet const& pset);
    virtual ~NoiseFilter();

    void produce(art::Event & evt);
    void reconfigure(fhicl::ParameterSet const& pset);    
    void beginJob();
    void endJob();

    void calcMinMax( std::vector<short> const& rawData,
                     float&                    pedestal,
                     float&                    noise,
                     float&                    maxADC,
                     float&                    minADC,
                     int&                      maxADCbin,
                     int&                      minADCbin,
                     int const&                startBin,
                     int const&                endBin);

    void CalcMeanRMSWithFlags(TH1F* fHist, double &theMean, double &theRMS);
    void ChirpFilterAlg(TH1F* fHist,float & isChirpFrac);
    void simpleZigzagFilterAlg(TH1F* fHist);
    void ZigzagFilterAlg(TH1F* fHist);
    void SignalFilterAlg(TH1F* fHist);
    void NoisyFilterAlg(TH1F* fHist, int planeNum);
    void TransientNoiseFilterAlg(TH1F* fHist, int planeNum);
    void WaveFilterAlg(TH1F **filtHists);
    void RawAdaptiveBaselineAlg(TH1F *filtHist);
    void FinalNoisyFilterAlg(TH1F* fHist, int planeNum);
    void RemoveFilterFlags(TH1F *filtHist);
    
  private:

    //Keep Track of event number
    int fEvtNum; //Number of current event

    //******************************
    //Variables Taken from FHICL File
    std::string       fRawDigitModuleLabel;   //label for rawdigit module
    int fDoChirp, fDoZigZag, fDoSignalFilter, fDoNoisyFilter,fDoTransientNoiseFilter, fSaveFilterWF, fSaveNoiseTree, fSaveNoiseTreeWF, fRunWaveFilterAlg, fRunRawAdaptiveBaselineAlg, fRemoveFilterFlags, fDoFinalNoisyFilter;

    TTree *fNoise;
    int fEvent, fRun,fSubrun,fChan, fPlaneNum;
    float fMax, fMin, fMean, fRms;
    int fMaxtime, fMintime;
    bool fisChirp;
    float fZigZagVal, fChirp;
    std::vector<short> fWf;
      
    int fMaxTicks; ///< maximum number of ticks expected (should get from RawDigit in future)

    TH1F *fHistMod;
    TH1F** waveNoiseHists;
    TH1F *currentHist;
    TH1F *currentFFTHist;

    const lariov::IDetPedestalProvider&  fPedestalRetrievalAlg; ///< Keep track of an instance to the pedestal retrieval alg
  }; //end class Noise


  //-------------------------------------------------------------------
  NoiseFilter::NoiseFilter(fhicl::ParameterSet const& pset)
    : EDProducer(),
    fMaxTicks(9595),
    fPedestalRetrievalAlg(art::ServiceHandle<lariov::IDetPedestalService>()->GetPedestalProvider())
    {
    this->reconfigure(pset);
    produces<std::vector<raw::RawDigit> >();
  }

  //-------------------------------------------------------------------
  NoiseFilter::~NoiseFilter(){}

  //-------------------------------------------------------------------
  void NoiseFilter::reconfigure(fhicl::ParameterSet const& pset){
    fRawDigitModuleLabel       = pset.get<std::string>("RawDigitModuleLabel");
    fDoChirp                   = pset.get<int>("doChirp",                  1);
    fDoZigZag                  = pset.get<int>("doZigZag",                 1);
    fDoSignalFilter            = pset.get<int>("doSignalFilter",           1);
    fDoNoisyFilter             = pset.get<int>("doNoisyFilter",            1);
    fRunWaveFilterAlg 	       = pset.get<int>("doWaveFilter",             1);
    fRunRawAdaptiveBaselineAlg = pset.get<int>("doAdaptiveBaseline",       1);
    fDoFinalNoisyFilter        = pset.get<int>("doFinalNoisyFilter",       1);
    fRemoveFilterFlags         = pset.get<int>("doRemoveFilterFlags",      1);
    fSaveFilterWF              = pset.get<int>("saveFilterWF",             1);
    fSaveNoiseTree 	           = pset.get<int>("saveNoiseTree",            0);
    fSaveNoiseTreeWF 	       = pset.get<int>("saveNoiseTreeWF",          0);
    fDoTransientNoiseFilter    = pset.get<int>("doTransientNoiseFilter",   0);
  }

  //-------------------------------------------------------------------
  void NoiseFilter::beginJob(){
    //art::ServiceHandle<geo::Geometry> geo;
    art::ServiceHandle<art::TFileService> tfs;

    //noise filter variables
    fHistMod = new TH1F("fHistMod","",fMaxTicks,-0.5,fMaxTicks-0.5);
    currentHist = new TH1F("","",fMaxTicks,-0.5,fMaxTicks-0.5);
    waveNoiseHists = new TH1F * [48];
    for(unsigned i = 0; i < 48; i++)
    	waveNoiseHists[i] = new TH1F(Form("waveNoiseHist_number%d",i),";Tick;ADC Value",fMaxTicks,-0.5,fMaxTicks-0.5);

    //noise tree variable
    fNoise = tfs->make<TTree>("NoiseTree","NoiseTree");
    fNoise->Branch("run", &fRun, "run/s");
    fNoise->Branch("subrun", &fSubrun, "subrun/s");
    fNoise->Branch("event", &fEvent, "event/s");
    fNoise->Branch("chan", &fChan, "chan/s");
    fNoise->Branch("max", &fMax, "max/F");
    fNoise->Branch("min", &fMin, "min/F");
    fNoise->Branch("mean", &fMean, "mean/F");
    fNoise->Branch("rms", &fRms, "rms/F");
    fNoise->Branch("chirp", &fChirp, "chirp/F");
    fNoise->Branch("wf", "vector<short>", &fWf);
  }

  //-------------------------------------------------------------------
  void NoiseFilter::endJob(){
    art::ServiceHandle<art::TFileService> tfs;
  }
  
  //-------------------------------------------------------------------
  void NoiseFilter::produce(art::Event & evt){
    LOG_INFO ("NoiseFilter Module") << "Processing Run " << fRun << ", Subrun " << fSubrun << ", Event " << fEvent;

    auto const* fGeometry = lar::providerFrom<geo::Geometry>();
  //  auto const* fDetectorProperties = lar::providerFrom<util::DetectorPropertiesService>();
    
    //filtered raw digits	
    std::unique_ptr<std::vector<raw::RawDigit> > filteredRawDigit(new std::vector<raw::RawDigit>);

    //get raw digit container
    art::Handle< std::vector<raw::RawDigit> > rawDigitHandle;
    evt.getByLabel(fRawDigitModuleLabel,rawDigitHandle);
    std::vector<raw::RawDigit> const& rawDigitVector(*rawDigitHandle);

    //get some event specific variables
    fRun = evt.run();
    fSubrun = evt.subRun();
    fEvent = evt.event();

    //define max # channels, max # channels per plane
    unsigned int maxChannels  = fGeometry->Nchannels();
    unsigned int wireMaxNum[] = {fGeometry->Nwires(0),fGeometry->Nwires(1),fGeometry->Nwires(2)};
    //unsigned int maxTimeSamples = fDetectorProperties->NumberTimeSamples();

    //define array to store channel numbers corresponding to wire plane, number, really terrible implementation
    std::vector<std::vector<int>> wirePlaneNum;
    wirePlaneNum.resize(fGeometry->Nplanes());
      
    for(size_t idx = 0; idx < fGeometry->Nplanes(); idx++)
        wirePlaneNum[idx].resize(fGeometry->Nwires(idx),-1);

    //define variables related to correlated noise hists
    std::vector< short > waveform;
    int waveNoiseGroupNum  = 48;
    Int_t waveNoiseCounter = -1;
    int waveNoiseHistsCh[waveNoiseGroupNum];
    float waveNoiseHistsChirp[waveNoiseGroupNum];

    //loop over channels in raw digit container, get mapping to wire plane and number
    const unsigned int n_channels = rawDigitVector.size();
    for(unsigned int ich=0; ich<n_channels; ich++)
    {
        raw::ChannelID_t channel = rawDigitVector.at(ich).Channel();
        bool goodChan(true);

        //get wire IDs
        std::vector<geo::WireID> wids;
        try {
            wids = fGeometry->ChannelToWire(channel);
        }
        catch(...)
        {
            goodChan = false;
        }
        if (channel >= maxChannels || !goodChan) continue;
        
        // Recover plane and wire in the plane
        unsigned int view = wids[0].Plane;
        unsigned int wire = wids[0].Wire;
        if( view < 0 || view > fGeometry->Nplanes()) continue;
        if( wire < 0 || wire > wireMaxNum[view] ) continue;
        wirePlaneNum[view][wire] = ich;
    }//end loop over raw digit container channels

    //loop over planes, wire #, apply noise filtering
    for(unsigned int pnum = 0 ; pnum < 3 ; pnum++){
        for(unsigned int wnum = 0; wnum < wireMaxNum[pnum] ; wnum++){
            int ich = wirePlaneNum[pnum][wnum];
            
            if( ich < 0 ) continue;
        
            //Setup channel specific variables
            const size_t n_samp = rawDigitVector.at(ich).NADC();
            fPlaneNum = (int) pnum;
            fHistMod->Reset();
            waveNoiseCounter++;
      		waveNoiseHists[waveNoiseCounter]->Reset();
            waveNoiseHistsCh[waveNoiseCounter] = rawDigitVector.at(ich).Channel();
            waveNoiseHistsChirp[waveNoiseCounter] = 0;

            //load current waveform into fHist object for modifications
            for( unsigned int s = 0 ; s < n_samp ; s++ )
                fHistMod->SetBinContent(s+1, rawDigitVector.at(ich).ADCs().at(s) );

            //filter out chirping
            float chirpVal = 0;
            if(fDoChirp == 1) ChirpFilterAlg( fHistMod, chirpVal );
            waveNoiseHistsChirp[waveNoiseCounter] = chirpVal;

            //filter zig zag (note flag applied in method)
			simpleZigzagFilterAlg(fHistMod);
			//ZigzagFilterAlg();
	
            //protect signals
            if(fDoSignalFilter == 1) SignalFilterAlg(fHistMod);

            //filter very noisy channels
            if(fDoNoisyFilter == 1) NoisyFilterAlg(fHistMod, (int) pnum);

            //filter transient noise
            if( fDoTransientNoiseFilter == 1 ) TransientNoiseFilterAlg(fHistMod, (int) pnum );
		
            //load partially filtered waveform into "wave group" placeholder
            for( int s = 0 ; s < fHistMod->GetNbinsX() ; s++ )
                waveNoiseHists[waveNoiseCounter]->SetBinContent(s+1, fHistMod->GetBinContent(s+1) );

            //do correlated noise removal
      		if(waveNoiseCounter == waveNoiseGroupNum-1)
      		{
                if( fRunWaveFilterAlg == 1)
                    WaveFilterAlg(waveNoiseHists);
                for(Int_t k = 0; k < waveNoiseGroupNum; k++)
        		{
                    if( fRunRawAdaptiveBaselineAlg == 1)
                        RawAdaptiveBaselineAlg(waveNoiseHists[k]);

                    //get mean for waveform
                    double meanVal, rmsVal;
                    CalcMeanRMSWithFlags(waveNoiseHists[k], meanVal, rmsVal);
                    for( unsigned int s = 0 ; s < n_samp ; s++ ){
                        if( waveNoiseHists[k]->GetBinContent( s+1 ) != 10000.0 )
                            waveNoiseHists[k]->SetBinContent(s+1, waveNoiseHists[k]->GetBinContent( s+1 ) - meanVal );
                    }

                    if( fDoFinalNoisyFilter == 1 )
                        FinalNoisyFilterAlg( waveNoiseHists[k] , (int) pnum);

                    //remove flags
                    if( fRemoveFilterFlags == 1)
                        RemoveFilterFlags(waveNoiseHists[k]);

                    //output waveform to file
                    waveform.clear();
                    unsigned int n_samp =  waveNoiseHists[k]->GetNbinsX();
                    // Recover the database version of the pedestal
                    double pedVal = fPedestalRetrievalAlg.PedMean(waveNoiseHistsCh[k]);

                    for( unsigned int s = 0 ; s < n_samp ; s++ ){
                        waveform.push_back( waveNoiseHists[k]->GetBinContent( s+1 ) + pedVal - meanVal );
                    }

                    //do basic noise measurement, save in tree
                    fChan = waveNoiseHistsCh[k];
                    fChirp = waveNoiseHistsChirp[k];
                    calcMinMax( waveform, fMean , fRms, fMax, fMin, fMaxtime, fMintime, 0, (int) n_samp);
                    fRms = rmsVal;
                    fMean = meanVal;

                    fWf.clear();
                    if( fSaveNoiseTreeWF == 1 ){
                        for( unsigned int s = 0 ; s < n_samp ; s++ ){
                            fWf.push_back(  waveNoiseHists[k]->GetBinContent( s+1 ) );
                            //fWf.push_back( fHistMod->GetBinContent( s+1 ) );
                        }
                    }
                    if( fSaveNoiseTree == 1)
                        fNoise->Fill();
				
                    if( fSaveFilterWF == 1 ){
                        filteredRawDigit->emplace_back( raw::RawDigit( waveNoiseHistsCh[k] , n_samp, waveform, raw::kNone) );
                        filteredRawDigit->back().SetPedestal(fMean+pedVal,fRms);
                    }
                }
                waveNoiseCounter = -1;
            }

            //place new raw digit in filtered continer, update pedestal and rms measurment
            //filteredRawDigit->emplace_back( raw::RawDigit(fchan, n_samp, rawDigitVector.at(ich).ADCs(), raw::kNone) ); //pipe out output
            //filteredRawDigit->emplace_back( raw::RawDigit(fchan, n_samp, waveform, raw::kNone) );
        }//end wire # loop
    }//end plane # loop

    if( fSaveFilterWF == 1 )
    	evt.put(std::move(filteredRawDigit));

    //if(hasCompressedRawDigit(rawDigitVector))
    //  throw "ERORR! You can't run the CalibrationTPCTest analyzer with compressed rawDigits!";
  }

  //---------------------------------------
  void NoiseFilter::calcMinMax( std::vector<short> const& rawData,
                               float & pedestal,
                               float & noise,
                               float & maxADC,
                               float & minADC,
                               int & maxADCbin,
                               int & minADCbin,
                               int const& startBin,
                               int const& endBin){

    	int n_samples = (int) rawData.size();
    
    	if( startBin >= endBin || startBin < 0 || endBin > n_samples ) return;

    	pedestal  = 0;
    	noise     = 0;
    	maxADC    = 0;
    	minADC    = 4095;
    	maxADCbin = 0;
    	minADCbin = 0;

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

    	//calculate RMS on pre-pulse to make sure it is not too high (what threshold?)
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

void NoiseFilter::CalcMeanRMSWithFlags(TH1F* fHist, double &theMean, double &theRMS)
{
  Double_t ADCval;
  theMean = 0.0;
  theRMS  = 0.0;
  Int_t waveformSize = fHist->GetNbinsX();
  Int_t counter = 0;
  for(Int_t i = 0; i < waveformSize; i++)
  {
    ADCval = fHist->GetBinContent(i+1);

    if(ADCval < 4096.0)
    {
      theMean += ADCval;
      theRMS += TMath::Power(ADCval,2.0);
      counter++;
    }
  }
  
  if(counter == 0)
  {
    theMean = 0.0;
    theRMS = 0.0;
  }
  else
  {
    theMean /= (Double_t)counter;
    theRMS /= (Double_t)counter;
    theRMS = TMath::Sqrt(theRMS-TMath::Power(theMean,2.0));
  }

  return;
}

void NoiseFilter::ChirpFilterAlg(TH1F* fHist, float & isChirpFrac)
{
  const Int_t windowSize = 20;
  const Double_t chirpMinRMS = 0.9;
  const Double_t maxNormalNeighborFrac = 0.20;
  bool recoverChirpingWaveforms = true;
  const Int_t maxTicks = fMaxTicks;

  Int_t counter = 0;
  Double_t ADCval;
  Double_t runningAmpMean = 0.0;
  Double_t runningAmpRMS = 0.0;
  Int_t numLowRMS = 0;
  Int_t firstLowRMSBin = -1;
  Int_t lastLowRMSBin = -1;
  Bool_t lowRMSFlag = false;
  Double_t RMSfirst = 0.0;
  Double_t RMSsecond = 0.0;
  Double_t RMSthird = 0.0;
  Int_t numNormalNeighbors = 0;
  Int_t numBins = fHist->GetNbinsX();
  for(Int_t i = 0; i < numBins; i++)
  {
    ADCval = fHist->GetBinContent(i+1);
    runningAmpMean += ADCval;
    runningAmpRMS += TMath::Power(ADCval,2.0);

    counter++;
    if(counter == windowSize)
    {
      runningAmpMean /= (Double_t)windowSize;
      runningAmpRMS /= (Double_t)windowSize;
      runningAmpRMS = TMath::Sqrt(runningAmpRMS-TMath::Power(runningAmpMean,2.0));

      RMSfirst = RMSsecond;
      RMSsecond = RMSthird;
      RMSthird = runningAmpRMS;

      if(runningAmpRMS < chirpMinRMS)
      {
        numLowRMS++;
        if(lowRMSFlag == false)
 	{
          lowRMSFlag = true;
          firstLowRMSBin = i-windowSize+1;
          lastLowRMSBin = i-windowSize+1;
	}
	else
	{
          lastLowRMSBin = i-windowSize+1;
	}
      }

      if(i >= 3*windowSize-1)
      {
        if((RMSsecond < chirpMinRMS) && ((RMSfirst > chirpMinRMS) || (RMSthird > chirpMinRMS)))
          numNormalNeighbors++;
      }

      counter = 0;
      runningAmpMean = 0.0;
      runningAmpRMS = 0.0;
    }
  }

  Double_t chirpFrac = ((Double_t) numLowRMS)/(((Double_t) maxTicks)/((Double_t) windowSize));
  Double_t normalNeighborFrac = ((Double_t) numNormalNeighbors)/((Double_t) numLowRMS);

  if(((normalNeighborFrac < maxNormalNeighborFrac) || ((numLowRMS < 2.0/maxNormalNeighborFrac) && (lastLowRMSBin-firstLowRMSBin == numLowRMS*windowSize))) && (numLowRMS > 4))
  {
    firstLowRMSBin = TMath::Max(1,firstLowRMSBin-windowSize);
    lastLowRMSBin = TMath::Min(numBins,lastLowRMSBin+2*windowSize);
    isChirpFrac = chirpFrac;

    if((numBins-lastLowRMSBin) < windowSize)
    {
      lastLowRMSBin = numBins;
    }

    if(chirpFrac > 0.99)
    {
      firstLowRMSBin = 1;
      lastLowRMSBin = numBins;
    }

    if(recoverChirpingWaveforms == true)
    {
      for(Int_t i = 0; i < numBins; i++)
      {
        if((i+1 >= firstLowRMSBin) && (i+1 <= lastLowRMSBin))
        {
          fHist->SetBinContent(i+1,10000.0);
        }
      }
    }
    else
    {
      for(Int_t i = 0; i < numBins; i++)
      {
        fHist->SetBinContent(i+1,10000.0);
      }
    }
  }

  return;
}

void NoiseFilter::simpleZigzagFilterAlg(TH1F* fHist){

  // ** TU: I think the logic needs to be rethought here since this method is doing two things
  //get mean for waveform
  double meanVal, rmsVal;
  CalcMeanRMSWithFlags(fHist, meanVal, rmsVal);

  // In case we want to turn zig zag off...
  if(fDoZigZag == 1)
  {
      Double_t ADCval, nextADCval;
      Int_t numBins = fHist->GetNbinsX();
      for(Int_t i = 0; i < numBins-1; i++)
      {
          ADCval = fHist->GetBinContent(i+1);
          nextADCval = fHist->GetBinContent(i+2);
          if( ADCval < 4096 && nextADCval < 4096 )
              fHist->SetBinContent(i+1, (ADCval + nextADCval)/2. - meanVal);
      }
      fHist->SetBinContent(numBins, fHist->GetBinContent(numBins) - meanVal );
  }
  // Otherwise, we still need to baseline subtract here
  else
  {
      for(Int_t idx = 0; idx < fHist->GetNbinsX(); idx++)
          fHist->SetBinContent(idx, fHist->GetBinContent(idx) - meanVal);
  }
}

void NoiseFilter::ZigzagFilterAlg(TH1F* fHist)
{
  const Int_t startFiltBin = -1;
  const Int_t endFiltBin = 3700;
  const Int_t maxTicks = fMaxTicks;

  double meanVal, rmsVal;
  CalcMeanRMSWithFlags(fHist, meanVal, rmsVal);
  if(rmsVal < 0.5) return;

  Double_t waveformMean = 0.0;
  Double_t ADCval;
  Int_t counter = 0;
  Int_t numBins = fHist->GetNbinsX();
  for(Int_t i = 0; i < numBins; i++)
  {    
    ADCval = fHist->GetBinContent(i+1);
    if(ADCval != 10000.0)
    {
      waveformMean += (Double_t) ADCval;
      counter++;
    }
  }
  if(counter > 0)
  {
    waveformMean /= ((Double_t) counter);
  }

  TH1F *currentHist = new TH1F("","",numBins,-0.5,numBins-0.5);
  for(Int_t i = 0; i < numBins; i++)
  {
    ADCval = fHist->GetBinContent(i+1);
    if(ADCval != 10000.0)
    {
      currentHist->SetBinContent(i+1,ADCval);
    }
    else
    {
      currentHist->SetBinContent(i+1,waveformMean);
    }
  }

  Double_t *reFilt = new Double_t[maxTicks];
  Double_t *imFilt = new Double_t[maxTicks];

  TH1 *hm = 0;
  hm = currentHist->FFT(0,"MAG");

  TH1 *hp = 0;
  hp = currentHist->FFT(0,"PH");

  for(Int_t i = 0; i < numBins; i++)
  {
      Double_t rho = hm->GetBinContent(i+1);
      Double_t phi = hp->GetBinContent(i+1);

    if((TMath::Min(i+1,numBins-i) > startFiltBin) && (TMath::Min(i+1,numBins-i) < endFiltBin))
    {
      reFilt[i] = (rho*cos(phi))/numBins;
      imFilt[i] = (rho*sin(phi))/numBins;
    }
    else
    {
      reFilt[i] = 0.0;
      imFilt[i] = 0.0;
    }
  }
  reFilt[0] = 0.0;
  imFilt[0] = 0.0;

  Int_t nFreqBins = numBins;
  TVirtualFFT *invCurrentFFTObject = TVirtualFFT::FFT(1,&nFreqBins,"C2R M K");
  invCurrentFFTObject->SetPointsComplex(reFilt,imFilt);
  invCurrentFFTObject->Transform();
  TH1F *newHist = new TH1F("","",numBins,-0.5,numBins-0.5);
  newHist = (TH1F*)TH1::TransformHisto(invCurrentFFTObject,0,"Re");

  for(Int_t i = 0; i < numBins; i++)
  {
    ADCval = fHist->GetBinContent(i+1);

    if(ADCval != 10000.0)
    {
      fHist->SetBinContent(i+1,newHist->GetBinContent(i+1));
    }
  }

  delete hm;
  delete hp;
  delete currentHist;
  delete newHist;
  delete invCurrentFFTObject;
  delete[] reFilt;
  delete[] imFilt;

  return;
}

void NoiseFilter::SignalFilterAlg(TH1F* fHist)
{
  const Double_t sigFactor = 4.0;
  const Int_t padBins = 8;

  double meanVal, rmsVal;
  CalcMeanRMSWithFlags(fHist, meanVal, rmsVal);
  Double_t sigThreshold = sigFactor*rmsVal;

  Double_t ADCval;
  std::vector<Bool_t> signalRegions;
  Int_t numBins = fHist->GetNbinsX();
  for(Int_t i = 0; i < numBins; i++)
  {
    ADCval = fHist->GetBinContent(i+1);

    if(((ADCval > sigThreshold + meanVal) || (ADCval < -1.0*sigThreshold + meanVal)) && (ADCval < 4096.0))
    {
      signalRegions.push_back(true);
    }
    else
    {
      signalRegions.push_back(false);
    }
  }

  for(Int_t i = 0; i < numBins; i++)
  {
    if(signalRegions[i] == true)
    {
      for(Int_t j = TMath::Max(0,i-padBins); j < TMath::Min(numBins,i+padBins); j++)
      {
        ADCval = fHist->GetBinContent(j+1);
        if(ADCval < 4096.0)
	{
          fHist->SetBinContent(j+1,ADCval+20000.0);
	}
      }
    }
  }  

  return;
}

void NoiseFilter::NoisyFilterAlg(TH1F* fHist, int planeNum)
{
  double meanVal, rmsVal;
  CalcMeanRMSWithFlags(fHist, meanVal, rmsVal);
  const Double_t maxRMSCut[3] = {10.0,10.0,5.0}; //hardcoded...
  //const Double_t maxRMSCut[3] = {6.0,6.0,4.0}; //hardcoded...
  if( planeNum < 0 || planeNum >= 3 )
	return;
  if(rmsVal > maxRMSCut[planeNum])
  {
    Int_t numBins = fHist->GetNbinsX();
    for(Int_t i = 0; i < numBins; i++)
    {
      fHist->SetBinContent(i+1,10000.0);
    }                          
  }

  return;
}

void NoiseFilter::FinalNoisyFilterAlg(TH1F* fHist, int planeNum)
{
  double meanVal, rmsVal;
  CalcMeanRMSWithFlags(fHist, meanVal, rmsVal);
  const Double_t maxRMSCut[3] = {2.5,2.5,1.6}; //hardcoded...
  //const Double_t maxRMSCut[3] = {6.0,6.0,4.0}; //hardcoded...
  if( planeNum < 0 || planeNum >= 3 )
	return;
  if(rmsVal > maxRMSCut[planeNum])
  {
    Int_t numBins = fHist->GetNbinsX();
    for(Int_t i = 0; i < numBins; i++)
    {
      fHist->SetBinContent(i+1,10000.0);
    }                          
  }

  return;
}

void NoiseFilter::TransientNoiseFilterAlg(TH1F* fHist, int planeNum)
{
  /*
  double meanVal, rmsVal;
  const Double_t maxRMSCut[3] = {10.0,10.0,5.0}; //hardcoded...
  if( planeNum < 0 || planeNum >= 3 )
	return;

  //loop over waveform, caluclaute RMS in 

  if(rmsVal > maxRMSCut[planeNum])
  {
    Int_t numBins = fHist->GetNbinsX();
    for(Int_t i = 0; i < numBins; i++)
    {
      fHist->SetBinContent(i+1,10000.0);
    }                          
  }
  */

  return;
}

void NoiseFilter::WaveFilterAlg(TH1F **filtHists)
{
  Int_t numBins = filtHists[0]->GetNbinsX();
  Double_t ADCval;
  std::vector<Double_t> corrVals;
  Double_t correction;
  Int_t corrValSize;
  int waveNoiseGroupNum = 48;

  //get mean value for each hist
  //double meanValuePerCh[waveNoiseGroupNum];
  //for(Int_t j = 0; j < waveNoiseGroupNum; j++)
  //{
  //	double meanVal, rmsVal;
  //	CalcMeanRMSWithFlags(filtHists[j], meanVal, rmsVal);
  //	meanValuePerCh[j] = meanVal;
  //}

  for(Int_t i = 0; i < numBins; i++)
  {
    corrVals.clear();

    for(Int_t j = 0; j < waveNoiseGroupNum; j++)
    {
      ADCval = filtHists[j]->GetBinContent(i+1);
      if(ADCval < 4096.0)
      {
        corrVals.push_back(ADCval);
      }
    }

    corrValSize = corrVals.size();
    sort(corrVals.begin(),corrVals.end());

    if(corrValSize < 2)
    {
      correction = 0.0;
    }
    else if((corrValSize % 2) == 0)
    {
      correction = (corrVals[corrValSize/2] + corrVals[(corrValSize/2)-1])/2.0;
    }
    else
    {
      correction = corrVals[(corrValSize-1)/2];
    }

    for(Int_t j = 0; j < waveNoiseGroupNum; j++)
    {
      ADCval = filtHists[j]->GetBinContent(i+1);
      if(ADCval != 10000.0)
      {
	//std::cout << corrValSize << "\t" << ADCval << "\t" << correction << "\t" << TMath::Nint(ADCval-correction+ meanValuePerCh[j]) << std::endl;
        //filtHists[j]->SetBinContent(i+1,TMath::Nint(ADCval-correction + meanValuePerCh[j]));
	filtHists[j]->SetBinContent(i+1,TMath::Nint(ADCval-correction ));
      }
    }//end loop over wave noise group
  }//end loop over bins

  return;
}

void NoiseFilter::RawAdaptiveBaselineAlg(TH1F *filtHist)
{
  const Int_t windowSize = 20;

  Int_t numBins = filtHist->GetNbinsX();
  Int_t minWindowBins = windowSize/2;

  Double_t baselineVec[numBins];
  Bool_t isFilledVec[numBins];

  Int_t numFlaggedBins = 0;
  for(Int_t j = 0; j < numBins; j++)
  {
    if(filtHist->GetBinContent(j+1) == 10000.0)
    {
      numFlaggedBins++;
    }
  }
  if(numFlaggedBins == numBins) return; // Eventually replace this with flag check

  Double_t baselineVal = 0.0;
  Int_t windowBins = 0;
  //Int_t index;
  Double_t ADCval;
  for(Int_t j = 0; j <= windowSize/2; j++)
  {
    ADCval = filtHist->GetBinContent(j+1);
    if(ADCval < 4096.0)
    {
      baselineVal += ADCval;
      windowBins++;
    }
  }

  if(windowBins == 0)
    baselineVec[0] = 0.0;
  else
    baselineVec[0] = baselineVal/((Double_t) windowBins);
  
  if(windowBins < minWindowBins)
    isFilledVec[0] = false;
  else
    isFilledVec[0] = true;
  
  Int_t oldIndex;
  Int_t newIndex;
  for(Int_t j = 1; j < numBins; j++)
  {
    oldIndex = j-windowSize/2-1;
    newIndex = j+windowSize/2;

    if(oldIndex >= 0)
    {
      ADCval = filtHist->GetBinContent(oldIndex+1);
      if(ADCval < 4096.0)
      {
        baselineVal -= filtHist->GetBinContent(oldIndex+1);
        windowBins--;
      }
    }

    if(newIndex < numBins)
    {  
      ADCval = filtHist->GetBinContent(newIndex+1);
      if(ADCval < 4096)
      {
        baselineVal += filtHist->GetBinContent(newIndex+1);
        windowBins++;
      }
    }

    if(windowBins == 0)
      baselineVec[j] = 0.0;
    else
      baselineVec[j] = baselineVal/windowBins;
  
    if(windowBins < minWindowBins)
      isFilledVec[j] = false;
    else
      isFilledVec[j] = true;
  }

  //get mean for waveform
  double meanVal, rmsVal;
  CalcMeanRMSWithFlags(filtHist, meanVal, rmsVal);

  Int_t downIndex;
  Int_t upIndex;
  Bool_t downFlag;
  Bool_t upFlag;
  for(Int_t j = 0; j < numBins; j++)
  {
    downFlag = false;
    upFlag = false;

    ADCval = filtHist->GetBinContent(j+1);
    if(ADCval != 10000.0)
    //if(ADCval < 10000.0)
    {
      if(isFilledVec[j] == false)
      {
        downIndex = j;
        while((isFilledVec[downIndex] == false) && (downIndex > 0) && (filtHist->GetBinContent(downIndex+1) != 10000.0))
        {
          downIndex--;
        }
      
        if(isFilledVec[downIndex] == false)
          downFlag = true;
      
        upIndex = j;
        while((isFilledVec[upIndex] == false) && (upIndex < numBins-1) && (filtHist->GetBinContent(upIndex+1) != 10000.0))
        {
          upIndex++;
        }
      
        if(isFilledVec[upIndex] == false)
          upFlag = true;
      
        if((downFlag == false) && (upFlag == false))
          baselineVec[j] = ((j-downIndex)*baselineVec[downIndex]+(upIndex-j)*baselineVec[upIndex])/((Double_t) upIndex-downIndex);
        else if((downFlag == true) && (upFlag == false))
          baselineVec[j] = baselineVec[upIndex];
        else if((downFlag == false) && (upFlag == true))
          baselineVec[j] = baselineVec[downIndex];
        else
          baselineVec[j] = 0.0;
      }

      //std::cout << ADCval << "\t" << baselineVec[j] << "\t" << meanVal << std::endl;
      //filtHist->SetBinContent(j+1,ADCval-baselineVec[j]+meanVal);
      filtHist->SetBinContent(j+1,ADCval-baselineVec[j]);
    }
  }
}

void NoiseFilter::RemoveFilterFlags(TH1F *filtHist)
{
  Double_t ADCval;
  Int_t numBins = filtHist->GetNbinsX();
  for(Int_t i = 0; i < numBins; i++)
  {
    ADCval = filtHist->GetBinContent(i+1);

    if(ADCval > 4096.0)
    {
      if(ADCval > 10000.0)
        filtHist->SetBinContent(i+1,ADCval-20000.0);
      else
        filtHist->SetBinContent(i+1,0.0);
    }
  }

  return;
}


  DEFINE_ART_MODULE(NoiseFilter)

} //end namespace noise

#endif //NOISEFILTER_H

