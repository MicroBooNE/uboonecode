
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
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TTree.h>
#include <THStack.h>
#include <TGraph.h>
#include <TMath.h>
#include <TComplex.h>

#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Handle.h" 
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

#include "Utilities/LArFFT.h"
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
    void MakeHisto( TH1D* & histo1D,
		    TH2D* & histo2D,
		    std::vector< std::vector<float> > Data);
    void GetWaveforms( std::vector<raw::RawDigit> const& rawDigit,
		       std::vector<TH1I*> & Waveforms,
		       std::map< unsigned int, uint32_t > chanmap );

    void GetNoiseSpec( std::vector<raw::RawDigit> const& rawDigit,
		       std::vector<TH1D*> & NoiseSpec,
		       std::map< unsigned int, uint32_t > chanmap );
    
  private:

    //Keep Track of event number
    int fEvtNum; //Number of current event

    //******************************
    //Variables Taken from FHICL File
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
    std::vector<int>   fVsteps;
    std::vector<float> fVmin;
    std::vector<float> fVmax;
    std::vector<float> fVstep;
    std::vector<float> fLinFitVmaxColl;
    std::vector<float> fLinFitVmaxInd;
    int               fLinBegin;
    std::vector<int>  fSubRunGain;
    std::vector<int>  fSubRunShape;
    std::vector<float>fSubRunGainVal;
    std::vector<float>fSubRunShapeVal;
    int               fAverages;
    int               fGetNoiseSpec;
    int               fFFTSize;
    std::vector<int>  fEvtList;
    std::vector<int>  fChanList;
    //******************************

    // these are containers for the calibration results
    // Intended design: each of these is reinitialized at subrun begin
    // Thus, fPedestalData[ie][ich] = mean pedestal for event ie, channel ich
    std::vector< std::vector<float> > fPedestalData;
    std::vector< std::vector<float> > fNoiseData;
    std::vector< std::vector<float> > fGainData;
    std::vector< std::vector<float> > fGainAreaData;
    std::vector< std::vector<float> > fTimeData;
    //Create a Map to link ordering of channel [ich] in above vector
    //to channel number.
    std::vector< std::map< unsigned int, uint32_t> > fChanMap;
    std::vector< std::map< unsigned int, uint32_t> > fChanMapEndRun;
    
    //Vector to hold ADCs, freqSpec to then make average Waveform Plot
    std::vector< std::vector<double> > fWFADCHolder;
    std::vector< std::vector<double> > fWFFreqHolder;

    //Vector for Noise Spectrum (a vector of floats) per event per channel
    std::vector< std::vector< std::vector<float> > > fNoiseSpectra;

    //*************************************************
    //Output Histograms and MRT-Numbering 2D Maps etc..

    //TTrees:
    //------
    //TTree for baseline and noise values
    TTree* fDistributions;
    TTree* fBaselineTree;
    //Tree for gain distributions: one per Setting combination
    std::vector< TTree* > fGainDistributions;
    //Variables for TTrees
    Double_t fBaselineU, fBaselineV, fBaselineY;
    Double_t fNoiseU, fNoiseV, fNoiseY;
    Double_t fGainAmplitudeU, fGainAmplitudeV, fGainAmplitudeY;
    Double_t fGainAreaU, fGainAreaV, fGainAreaY;
    Double_t fBaselineCh_1;
    //Histograms one per subRun
    //-------------------------
    //Histo to count Events in subRun per channel
    std::vector<TH1I*> fAllChan;
    //Baseline Histogram and Map
    std::vector<TH1D*> fChannelBaseline;
    std::vector<TH2D*> fBaselineMap;
    //Noise Histogram and Map
    std::vector<TH1D*> fChannelNoise;
    std::vector<TH2D*> fNoiseMap;
    //Max ADC and Pulse Area Histograms and Maps
    std::vector<TH1D*> fChannelMaxADC;
    std::vector<TH1D*> fChannelArea;
    std::vector<TH2D*> fMaxADCMap;
    std::vector<TH2D*> fAreaMap;
    //Time-Tick of Max ADC Histograms and Maps
    std::vector<TH1D*> fChannelTimeMax;
    std::vector<TH2D*> fTimeMaxMap;
    //Gain Histograms and Maps
    std::vector<TH1D*> fChannelMaxADCAmplitude;
    std::vector<TH1D*> fChannelMaxADCArea;
    std::vector<TH2D*> fGainAmplitudeMap;
    std::vector<TH2D*> fGainAreaMap;
    //Vector for Waveforms and Noise Spectra Plots
    std::vector<TH1I*> fWaveforms; 
    std::vector<TH1D*> fNoiseSpec; 
    std::vector<THStack*> fWFOverlay;
    std::vector<THStack*> fNoiseSpecOverlay;
    //Linearity Histograms and Fits for Gain Measurements
    TH1D* fLinearityAmplitudeTmp;
    TH1D* fLinearityAreaTmp;
    TF1* fFitAmplitudeTmp;
    TF1* fFitAreaTmp;
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
    fVsteps              = pset.get<std::vector <int> >("Vsteps"); 
    fVmin                = pset.get<std::vector <float> >("Vmin"); 
    fVmax                = pset.get<std::vector <float> >("Vmax");
    fVstep               = pset.get<std::vector <float> >("Vstep");
    fLinFitVmaxColl      = pset.get<std::vector <float> >("LinFitVmaxColl");
    fLinFitVmaxInd       = pset.get<std::vector <float> >("LinFitVmaxInd");
    fLinBegin            = pset.get<int>("LinBegin");
    fSubRunGain          = pset.get< std::vector <int> >("SubRunGain");
    fSubRunShape         = pset.get< std::vector <int> >("SubRunShape");
    fSubRunGainVal       = pset.get< std::vector <float> >("SubRunGainVal");
    fSubRunShapeVal      = pset.get< std::vector <float> >("SubRunShapeVal");
    fAverages            = pset.get< int >("Averages");
    fGetNoiseSpec        = pset.get< int >("GetNoiseSpec");
    fFFTSize             = pset.get< int >("FFTSize");
    fEvtList             = pset.get< std::vector <int> >("EvtList");
    fChanList             = pset.get< std::vector <int> >("ChanList");
  }


  //-------------------------------------------------------------------
  void CalibrationTPC::beginJob(){
    
    art::ServiceHandle<geo::Geometry> geo;
    fChanBegin = 0;
    fNChannels = fNChanMax;
    fsubRunNum = 0;

  }

  //-------------------------------------------------------------------
  void CalibrationTPC::endJob(){
    //Output Message for who is running this module
    std::cout << std::endl << std::endl << std::endl
	      << "********************************************************************\n"
	      << "Module has endend. Now several scripts can be run to make additional\n"
	      << "plots such as overlays, residual plots, etc.\n"
	      << "if fShortRun was on then for event 1 a histogram of each channel was\n"
	      << "made. You can overlay all these waveforms by running the Overlay.cxx\n"
	      << "ROOT macro. It takes Calibrations.root as an input. If this module's\n"
	      << "output has a different name please change the script's input accord-\n"
	      << "-ingly.\n"
	      << "More to come...\n" 
	      << "********************************************************************\n"
	      << std::endl << std::endl << std::endl;
  }

  //-------------------------------------------------------------------
  void CalibrationTPC::beginRun(art::Run const& run){

    //In this function create some histograms and such that are 
    //filled once per run, such as gain measurements, etc...

    art::ServiceHandle<art::TFileService> tfs;
    gStyle->SetOptStat(0);
      
    if ( fEmptyRun ){

      fBaselineTree = tfs->make<TTree>("BaselineTree","BaselineTree");
      fBaselineTree->Branch("BaselineCh_1", &fBaselineCh_1, "BaselineCh_1/D");
    }

    if ( fGainRun ){
      
      for ( unsigned int count=0; count < fSubRunShape.size(); count++ ) {
	
	fGainDistributions.push_back( tfs->make<TTree>( Form("GainDistributions_ASICGain_%f_ShapingT_%f",
							     fSubRunGainVal.at(count),
							     fSubRunShapeVal.at(count) ) , 
							Form("GainDistributions_ASICGain_%f_ShapingT_%f",
							     fSubRunGainVal.at(count),
							     fSubRunShapeVal.at(count) ) ) );

	fGainDistributions.back()->Branch("GainAmplitudeU", &fGainAmplitudeU, "GainAmplitudeU/D");
	fGainDistributions.back()->Branch("GainAreaU", &fGainAreaU, "GainAreaU/D");
	fGainDistributions.back()->Branch("GainAmplitudeV", &fGainAmplitudeV, "GainAmplitudeV/D");
	fGainDistributions.back()->Branch("GainAreaV", &fGainAreaV, "GainAreaV/D");
	fGainDistributions.back()->Branch("GainAmplitudeY", &fGainAmplitudeY, "GainAmplitudeY/D");
	fGainDistributions.back()->Branch("GainAreaY", &fGainAreaY, "GainAreaY/D");
      }
    }
  }

  //-------------------------------------------------------------------
  void CalibrationTPC::endRun(art::Run const& run){

    art::ServiceHandle<art::TFileService> tfs;

    //Make Baseline-stability plots...go through all baseline values for all subruns
    if ( fGainRun || fEmptyRun ){
      Double_t avgBaseline = 0;
      Double_t RMSBaseline = 0;
      //Make a histogram per-channel of the 
      for ( int ch = 0; ch < fNChannels; ch++){//loop over all channels
	avgBaseline = 0;
	RMSBaseline = 0;
	if ( fAllChan.at(0)->GetBinContent(ch) > 0 ){//if this channel has data for this subrun
	  //temporary vector to hold baseline values for this channel
	  std::vector<float> AllBaselines;
	  // loop through subruns, find avg and RMS baseline
	  for ( uint32_t sr=0; sr < fsubRunNum; sr++ )
	    avgBaseline += fChannelBaseline.at(sr)->GetBinContent(ch);
	  avgBaseline /= fsubRunNum;
	  for ( uint32_t sr=0; sr < fsubRunNum; sr++ )
	    RMSBaseline += ( fChannelBaseline.at(sr)->GetBinContent(ch) - avgBaseline )*
	      (fChannelBaseline.at(sr)->GetBinContent(ch) - avgBaseline );
	  RMSBaseline = sqrt( RMSBaseline / (fsubRunNum-1) );
	}
      }
    }

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
	    fChannelMaxADCAmplitude.push_back( tfs->make<TH1D>( Form("GainFomAmplitude_ASICGain_%f_ShapingT_%f",
								     fSubRunGainVal.at(lincounter),
								     fSubRunShapeVal.at(lincounter)),
								"Gain from Pulse Height [ADCs/Volt]; Ch Num; Gain [ADCs/Volt]",
								fNChanMax+1, fChanBegin, fChanBegin+fNChanMax+1) );
	    //Gain 2D Histo -> see docDB 2007 slide 15
	    fGainAmplitudeMap.push_back( tfs->make<TH2D>( Form("GainMap_ASICGain_%f_ShapingT_%f",
							       fSubRunGainVal.at(lincounter),
							       fSubRunShapeVal.at(lincounter)),
							  "Gain from Pulse Height [ADCs/Volt]; Ch Num; FEM Num",
							  64, 0, 64,
							  12, 0, 12) );
	    //Prepare Gain Histograms: fill with linear fit value (slope)
	    fChannelMaxADCArea.push_back( tfs->make<TH1D>( Form("AreaVal_ASICGain_%f_ShapingT_%f",
								fSubRunGainVal.at(lincounter),
								fSubRunShapeVal.at(lincounter)),
							   "Gain from Pulse Area [ADCs x Ticks/Volt]; Ch Num; Gain [ADCs x Ticks / Volt]",
							   fNChanMax+1, fChanBegin, fChanBegin+fNChanMax+1) );
	    //Gain 2D Histo -> see docDB 2007 slide 15
	    fGainAreaMap.push_back( tfs->make<TH2D>( Form("AreaMap_ASICGain_%f_ShapingT_%f",
							  fSubRunGainVal.at(lincounter),
							  fSubRunShapeVal.at(lincounter)),
						     "Gain from Pulse Area [ADCs x Ticks/Volt]; Ch Num; FEM Num",
						     64, 0, 64,
						     12, 0, 12) );
	    //Number of steps in voltage:
	    int nsteps = fVsteps.at(lincounter);

	    for ( int ch = 0; ch < fNChannels; ch++){
	      //check and see if there are events for this chan number
	      if ( fAllChan.at(0)->GetBinContent(ch) > 0 ){ //means there are events

		fLinearityAmplitudeTmp = new TH1D("LinearityAmpTmp",
					    "Linearity Amp Tmp", 
					    nsteps,
					    fVmin.at(lincounter)-(fVstep.at(lincounter)/2.),
					    fVmax.at(lincounter)+(fVstep.at(lincounter)/2.) );
		fLinearityAreaTmp = new TH1D("LinearityAmpTmp",
					    "Linearity Amp Tmp", 
					    nsteps,
					    fVmin.at(lincounter)-(fVstep.at(lincounter)/2.),
					    fVmax.at(lincounter)+(fVstep.at(lincounter)/2.) );

		//now loop through next 11 subruns and fill histogram
		for (int i = 0; i < nsteps; i++){
		  //Need to avoid points with Vin < 50 mV from being plotted
		  //That is because the pulser used on MRT sends huge pulses if input set to < 50 mV.
		  //This will distort linearity plot and fit therefore affecting the gain measurement
		  //Vin = i*fVstep.at(lincounter)
		  if ( (i == 0) or (fVmin.at(lincounter)+i*fVstep.at(lincounter) >= 0.05) ) {
		    //if gain=3 then do not look at bin with V=0 to make life easier

		    if ( lincounter >= 6 ){
		      fLinearityAmplitudeTmp->SetBinContent( i+1, fChannelMaxADC.at(sr+i+1)->GetBinContent(ch) );
		      fLinearityAmplitudeTmp->SetBinError( i+1,  fChannelMaxADC.at(sr+i+1)->GetBinError(ch) );
		      fLinearityAreaTmp->SetBinContent( i+1, fChannelArea.at(sr+i+1)->GetBinContent(ch) );
		      fLinearityAreaTmp->SetBinError( i+1,  fChannelArea.at(sr+i+1)->GetBinError(ch) );
		    }
		    else {
		      fLinearityAmplitudeTmp->SetBinContent( i+1, fChannelMaxADC.at(sr+i)->GetBinContent(ch) );
		      fLinearityAmplitudeTmp->SetBinError( i+1,  fChannelMaxADC.at(sr+i)->GetBinError(ch) );
		      fLinearityAreaTmp->SetBinContent( i+1, fChannelArea.at(sr+i)->GetBinContent(ch) );
		      fLinearityAreaTmp->SetBinError( i+1,  fChannelArea.at(sr+i)->GetBinError(ch) );
		    }

		  }
		}//for all voltages
		
		//fit histogram
		//determine where max is (from LinFitVmax)
		//also depends on induction vs collection wire
		double fitmax;
		if ( (ch%100) < 32 ) { fitmax = fLinFitVmaxInd.at(lincounter);  }
		else                 { fitmax = fLinFitVmaxColl.at(lincounter); }

		//Create Fit objects for amplitude and area fits
		fFitAmplitudeTmp = new TF1(Form("Linearity_Ch_%d_Gain_%d_ShT_%d", ch, fSubRunGain.at(lincounter), fSubRunShape.at(lincounter))
					   ,"[0]+x*[1]", 0, fitmax );
		fFitAreaTmp = new TF1(Form("LinearityA_Ch_%d_Gain_%d_ShT_%d", ch, fSubRunGain.at(lincounter), fSubRunShape.at(lincounter))
				      ,"[0]+x*[1]", 0, fitmax );
		gStyle->SetOptFit(1);
		fFitAreaTmp->SetParName(0,"BaselineFit");
		fFitAreaTmp->SetParName(1,"GainFit");
		fFitAmplitudeTmp->SetParName(0,"BaselineFit");
		fFitAmplitudeTmp->SetParName(1,"GainFit");
		//Perform Fits
		fLinearityAmplitudeTmp->Fit( fFitAmplitudeTmp, "R" );
		fLinearityAreaTmp->Fit( fFitAreaTmp, "R" );
		//Add Gain Values to Histograms and Maps
		fChannelMaxADCAmplitude.back()->SetBinContent( ch+1, fFitAmplitudeTmp->GetParameter(1) );
		fGainAmplitudeMap.back()->SetBinContent( ch%100, int(ch/100)+1 , fFitAmplitudeTmp->GetParameter(1) ); 
		fChannelMaxADCArea.back()->SetBinContent( ch+1, fFitAreaTmp->GetParameter(1) );
		fGainAreaMap.back()->SetBinContent( ch%100, int(ch/100)+1 , fFitAreaTmp->GetParameter(1) ); 

		if ( ((ch%100) < 32) and ((ch%100)%2==0) ){
		  fGainAmplitudeU = fFitAmplitudeTmp->GetParameter(1);
		  fGainAreaU = fFitAreaTmp->GetParameter(1);
		}
		if ( ((ch%100) < 32) and ((ch%100)%2==1) ){
		  fGainAmplitudeV = fFitAmplitudeTmp->GetParameter(1);
		  fGainAreaV = fFitAreaTmp->GetParameter(1);
		}
		if ( (ch%100) >= 32 ){
		  fGainAmplitudeY = fFitAmplitudeTmp->GetParameter(1);
		  fGainAreaY = fFitAreaTmp->GetParameter(1);
		}

		fGainDistributions.at(lincounter)->Fill();
		
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
    //Mainly create histograms that are changed once per subrun (noise, baseline, maxADC, etc...)

    art::ServiceHandle<art::TFileService> tfs;
    gStyle->SetOptStat(0);

    fWFADCHolder.clear();
    fWFFreqHolder.clear();

    //fWaveforms.clear();
    //fNoiseSpec.clear();

    //First Make a channel map: number of events per channel
    fAllChan.push_back( tfs->make<TH1I>(Form("AllChan_SubRun_%d",fsubRunNum),
					"Events Recorded for Each Channel ;Chan Num; Num Evts", 
					fNChanMax+1, fChanBegin, fChanBegin+fNChanMax+1) );
    

    //if subrunnumber indicates pulse-type subrun
    if ( fGainRun or fShortRun ){

      fChannelMaxADC.push_back( tfs->make<TH1D>(Form("ChannelMaxADC_SubRun_%d",fsubRunNum),
					      "Max ADC count on Waveform by Channel ;Chan Num; Max [ADCs]",
					      fNChanMax+1, fChanBegin, fChanBegin+fNChanMax+1) );

      fChannelArea.push_back( tfs->make<TH1D>(Form("ChannelArea_SubRun_%d",fsubRunNum),
					      "ADC-Area by channel - Positive Pulse Only;Chan Num; Max [ADCs X Time-Ticks]",
					      fNChanMax+1, fChanBegin, fChanBegin+fNChanMax+1) );

      fMaxADCMap.push_back( tfs->make<TH2D>(Form("MaxADCMap_SubRun_%d",fsubRunNum),
					    "Max ADC Count on Waveform [ADCs] ;Chan Num; FEM Num",
					    64, 0, 64,
					    12, 0, 12) );

      fAreaMap.push_back( tfs->make<TH2D>(Form("AreaMap_SubRun_%d",fsubRunNum),
					  "ADC-Area by channel - [ADCs X Time-Ticks] ;Chan Num; FEM Num",
					  64, 0, 64,
					  12, 0, 12) );
      
      fChannelTimeMax.push_back( tfs->make<TH1D>(Form("ChannelTimeMax_SubRun_%d",fsubRunNum),
					      "Time-Tick of Max ADC count ;Chan Num; Max Time [Tick Number]",
					      fNChanMax+1, fChanBegin, fChanBegin+fNChanMax+1) );

      fTimeMaxMap.push_back( tfs->make<TH2D>(Form("TimeMaxMap_SubRun_%d",fsubRunNum),
					     "Time-Tick of Max ADC count [Time-Ticks];Chan Num; FEM Num",
					     64, 0, 64,
					     12, 0, 12) );

      fChannelBaseline.push_back( tfs->make<TH1D>(Form("Channelbaseline_SubRun_%d",fsubRunNum),
						  "Baseline By Channel ;Chan Num; Baseline [ADCs]",
						  fNChanMax+1, fChanBegin, fChanBegin+fNChanMax+1) );

      fBaselineMap.push_back( tfs->make<TH2D>(Form("BaselineMap_SubRun_%d",fsubRunNum),
					      "Baseline By Channel [ADCs];Chan Num; FEM Num",
					      64, 0, 64,
					      12, 0, 12) );
      
      fChannelNoise.push_back( tfs->make<TH1D>(Form("ChannelNoise_SubRun_%d",fsubRunNum),
					       "RMS Noise by Channel ;Chan Num; RMS Noise [ADCs]",
					       fNChanMax+1, fChanBegin, fChanBegin+fNChanMax+1) );
      
      fNoiseMap.push_back( tfs->make<TH2D>(Form("NoiseMap_SubRun_%d",fsubRunNum),
					   "RMS Noise by Channel [ADCs];Chan Num; FEM Num",
					   64, 0, 64,
					   12, 0, 12) );

    }

    if ( fEmptyRun ){

      fDistributions = tfs->make<TTree>(Form("Distributions_SubRun_%d",fsubRunNum),
					Form("Distributions_SubRun_%d",fsubRunNum) );

      fChannelBaseline.push_back( tfs->make<TH1D>(Form("Channelbaseline_SubRun_%d",fsubRunNum),
						  "Baseline By Channel ;Chan Num; Baseline [ADCs]",
						  fNChanMax+1, fChanBegin, fChanBegin+fNChanMax+1) );

      fBaselineMap.push_back( tfs->make<TH2D>(Form("BaselineMap_SubRun_%d",fsubRunNum),
					      "Baseline By Channel [ADCs];Chan Num; FEM Num",
					      64, 0, 64,
					      12, 0, 12) );

      fDistributions->Branch("BaselineU", &fBaselineU, "BaselineU/D");
      fDistributions->Branch("BaselineV", &fBaselineV, "BaselineV/D");
      fDistributions->Branch("BaselineY", &fBaselineY, "BaselineY/D");

      fChannelNoise.push_back( tfs->make<TH1D>(Form("ChannelNoise_SubRun_%d",fsubRunNum),
					       "RMS Noise by Channel ;Chan Num; RMS Noise [ADCs]",
					       fNChanMax+1, fChanBegin, fChanBegin+fNChanMax+1) );

      fNoiseMap.push_back( tfs->make<TH2D>(Form("NoiseMap_SubRun_%d",fsubRunNum),
					   "RMS Noise by Channel [ADCs];Chan Num; FEM Num",
					   64, 0, 64,
					   12, 0, 12) );

      fDistributions->Branch("NoiseU", &fNoiseU, "NoiseU/D");
      fDistributions->Branch("NoiseV", &fNoiseV, "NoiseV/D");
      fDistributions->Branch("NoiseY", &fNoiseY, "NoiseY/D");
            
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
      for ( size_t chanIndex=0; chanIndex < fChanMap[evtN].size(); chanIndex++ ){
	currentChan = fChanMap[evtN].at(chanIndex);
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
      MakeHisto( fChannelBaseline.at(fsubRunNum-1), fBaselineMap.at(fsubRunNum-1), fPedestalData);
      MakeHisto( fChannelNoise.at(fsubRunNum-1), fNoiseMap.at(fsubRunNum-1), fNoiseData);

      //Loop through vectors and fill TTree:
      for ( size_t ev = 0; ev < fPedestalData.size(); ev++ ){
	for ( size_t ch = 0; ch < fPedestalData.at(ev).size(); ch++) {
	  currentChan = fChanMap[ev].at(ch);
	  //if Channel == 1 then fill Branch with this channels baseline spread
	  if ( currentChan == 1 ){
	    fBaselineCh_1 = fPedestalData.at(ev).at(ch);
	    fBaselineTree->Fill();
	  }
	  //Make three different distributons for three different wire planes
	  if ( (currentChan%100) >= 32 ) {
	    fBaselineY = fPedestalData.at(ev).at(ch);
	    fNoiseY = fNoiseData.at(ev).at(ch);
	  }
	  //even ch = U plane. Odd = V plane
	  if ( ((currentChan%100) < 32) and ((currentChan%100)%2==0) ) {
	    fBaselineU = fPedestalData.at(ev).at(ch);
	    fNoiseU = fNoiseData.at(ev).at(ch);
	  }
	  if ( ((currentChan%100) < 32) and ((currentChan%100)%2==1) ) {
	    fBaselineV = fPedestalData.at(ev).at(ch);
	    fNoiseV = fNoiseData.at(ev).at(ch);
	  }
	  fDistributions->Fill();
	}
      }

    }//end emptyrun
    
    if ( fGainRun or fShortRun ){
      MakeHisto( fChannelMaxADC.at(fsubRunNum-1), fMaxADCMap.at(fsubRunNum-1), fGainData);
      MakeHisto( fChannelArea.at(fsubRunNum-1), fAreaMap.at(fsubRunNum-1), fGainAreaData);
      MakeHisto( fChannelBaseline.at(fsubRunNum-1), fBaselineMap.at(fsubRunNum-1), fPedestalData);
      MakeHisto( fChannelNoise.at(fsubRunNum-1), fNoiseMap.at(fsubRunNum-1), fNoiseData);
      MakeHisto( fChannelTimeMax.at(fsubRunNum-1), fTimeMaxMap.at(fsubRunNum-1), fTimeData);
    }//end gainrun and shortrun

    //Clear per-subrun vectors
    fChanMapEndRun = fChanMap;

    fPedestalData.clear();
    fNoiseData.clear();
    fNoiseSpectra.clear();
    fChanMap.clear();
    fGainData.clear();
    fGainAreaData.clear();

  }


  //-------------------------------------------------------------------
  void CalibrationTPC::MakeHisto( TH1D* & histo1D,
				  TH2D* & histo2D,
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
      DataErrPerChan.at(chanN) = sqrt( DataErrPerChan.at(chanN) / EvtsPerChan[chanN] );

    //fill histogram chan-by-chan
    for ( uint32_t chanN=0; chanN < fNChanMax; chanN++){
      if ( EvtsPerChan[chanN] > 0 ){
	histo1D->SetBinContent( chanN+1, DataPerChan.at(chanN));
	histo2D->SetBinContent( (chanN+1)%100, int(chanN/100)+1, DataPerChan.at(chanN) ); 
	histo1D->SetBinError( chanN+1, DataErrPerChan.at(chanN));
      }
    } 

  }

 
  //-------------------------------------------------------------------------
  //Plot Waveforms in Histograms
  void CalibrationTPC::GetWaveforms( std::vector<raw::RawDigit> const& rawDigit,
				     std::vector<TH1I*> & Waveforms,
				     std::map< unsigned int, uint32_t > chanmap ){


    if (fAverages and (fEvtNum == 1) ){
      std::vector<double> tmp(rawDigit.at(0).fADC.size(), 0.);
      for (unsigned int ch=0; ch < chanmap.size(); ch++)
	fWFADCHolder.push_back(tmp);
    }
    
    if ( fAverages ){
      for (unsigned int ich=0; ich < chanmap.size(); ich++){
	for ( unsigned int tick=0; tick < rawDigit.at(ich).fADC.size(); tick++)
	  fWFADCHolder.at(ich).at(tick) += float(rawDigit.at(ich).fADC.at(tick))/100.;
      }
    }
    

    if ( (fAverages) and (fEvtNum == 100) ){

      art::ServiceHandle<art::TFileService> tfs;
      //Make Histogram Stack
      fWFOverlay.push_back( tfs->make<THStack>( Form("WF_subrun_%d_ev_%d", fsubRunNum, fEvtNum ),
						       "Stack of Histograms" ) );

      for (unsigned int ich=0; ich < chanmap.size(); ich++){
	Waveforms.push_back( tfs->make<TH1I>(Form("WF_subRun_%d_ev_%d_chan_%d", fsubRunNum, fEvtNum, chanmap.at(ich) ), "Waveform Overlay (100 events); Time-Ticks [2 MHz - 500 ns]; ADCs", 
					     rawDigit.at(ich).fADC.size(),
					     0, rawDigit.at(ich).fADC.size() ) );

	for ( unsigned int tick=0; tick < fWFADCHolder.at(ich).size() ; tick++)
	  Waveforms.back()->SetBinContent( tick+1, fWFADCHolder.at(ich).at(tick) );
	fWFOverlay.back()->Add(Waveforms.back(), "HIST P");
      }
      fWFOverlay.back()->Draw();
    }



    if ( !fAverages ){
      
      art::ServiceHandle<art::TFileService> tfs;
      //Make Histogram Stack
      fWFOverlay.push_back( tfs->make<THStack>( Form("WF_subrun_%d_ev_%d", fsubRunNum, fEvtNum ),
						"Stack of Histograms" ) );
      
      for (unsigned int ich=0; ich < chanmap.size(); ich++){
	//only certain channels to make things faster
	if ( (std::find(fChanList.begin(), fChanList.end(), chanmap.at(ich)) != fChanList.end()) 
	     or (fChanList.at(0) == -1) ){
	  
	  Waveforms.push_back( tfs->make<TH1I>(Form("WF_subRun_%d_ev_%d_chan_%d", fsubRunNum, fEvtNum, chanmap.at(ich) ), "Waveform", 
					       rawDigit.at(ich).fADC.size(),
					       0, rawDigit.at(ich).fADC.size() ) );
	  
	  for ( unsigned int tick=0; tick < rawDigit.at(ich).fADC.size(); tick++){
	    Waveforms.back()->SetBinContent( tick+1, rawDigit.at(ich).fADC.at(tick) );
	  }
	  fWFOverlay.back()->Add(Waveforms.back(), "HIST P");
	}
      }//for all channels
      
      fWFOverlay.back()->Draw();
    }
    
    return;
    
  }
  


  //-----------------------------------------------------------------------------
  //Plot Noise Frequency Spectrum
  void CalibrationTPC::GetNoiseSpec( std::vector<raw::RawDigit> const& rawDigit,
				     std::vector<TH1D*> & NoiseSpec,
				     std::map< unsigned int, uint32_t > chanmap ){

    art::ServiceHandle<util::LArFFT> fFFT;
    std::string options = fFFT->FFTOptions();
    int fitbins = fFFT->FFTFitBins();
    fFFT->ReinitializeFFT(9600, options, fitbins);

    if (fAverages and (fEvtNum == 1) ){
      std::vector<double> tmp(rawDigit.at(0).fADC.size(), 0.);
      for (unsigned int ch=0; ch < chanmap.size(); ch++)
	fWFFreqHolder.push_back(tmp);
    }
    
    if ( fAverages ){
      
      for (unsigned int ich=0; ich < chanmap.size(); ich++){

	//make vectors for FFT analysis
	std::vector<float> noiseTime(rawDigit.at(ich).fADC.size(), 0.);
	std::vector<TComplex> noiseFrequency(rawDigit.at(ich).fADC.size(), 0.);
	std::vector<double> noiseFrequencyMag(rawDigit.at(ich).fADC.size(), 0.);

	for ( unsigned int tick=0; tick < rawDigit.at(ich).fADC.size(); tick++)
	  noiseTime.at(tick) += float(rawDigit.at(ich).fADC.at(tick));
	
	//do FFT
	fFFT->DoFFT( noiseTime, noiseFrequency);
	
	//get magnitude
	for (unsigned int f=0; f < noiseFrequency.size(); f++){
	    fWFFreqHolder.at(ich).at(f) += float(noiseFrequency.at(f).Rho());
	}
	
      }
    }
    

    if ( (fAverages) and (fEvtNum == 100) ){

      art::ServiceHandle<art::TFileService> tfs;
      //Make Histogram Stack
      fNoiseSpecOverlay.push_back( tfs->make<THStack>( Form("NoiseFreqSpec_subrun_%d_ev_%d", fsubRunNum, fEvtNum ),
						       "Stack of Histograms" ) );

      for (unsigned int ich=0; ich < chanmap.size(); ich++){
	std::cout << "Filling vector for chan: " << ich << std::endl;
	NoiseSpec.push_back( tfs->make<TH1D>(Form("NoiseFreqSpec_subRun_%d_ev_%d_chan_%d", fsubRunNum, fEvtNum, chanmap.at(ich) ),
					     "Avg Noise Frequency Spectrum over 100 Events; MHz; Amplitude", 
					     fFFT->FFTSize(), 0, 2) );
	//rawDigit.at(ich).fADC.size(), 0, 1) );

	for (int tick=0; tick < fFFT->FFTSize()/2+1 ; tick++){
	  if (tick == 100) { std::cout << "fWADCHolder at tick 100: " << fWFFreqHolder.at(ich).at(tick) << std::endl; }
	  NoiseSpec.back()->SetBinContent( tick+1, fWFFreqHolder.at(ich).at(tick)/100. );
	}

	fNoiseSpecOverlay.back()->Add(NoiseSpec.back(), "HIST P");

      }

      fNoiseSpecOverlay.back()->Draw();
    }


    if ( !fAverages ){
      
      art::ServiceHandle<art::TFileService> tfs;
      
      for (unsigned int ich=0; ich < chanmap.size(); ich++){
	//only selected channels
	if ( (std::find(fChanList.begin(), fChanList.end(), chanmap.at(ich)) != fChanList.end() )
	     or (fChanList.at(0) == -1) ){

	  NoiseSpec.push_back( tfs->make<TH1D>(Form("NoiseFreqSpec_subRun_%d_ev_%d_chan_%d", fsubRunNum, fEvtNum, chanmap.at(ich) ),
					       "Noise Frequency Spectrum; MHz; Amplitude", 
					       fFFT->FFTSize(), 0, 2 ) );

	  //make vectors for FFT analysis
	  std::vector<float> noiseTime(rawDigit.at(ich).fADC.size(), 0.);
	  std::vector<TComplex> noiseFrequency(rawDigit.at(ich).fADC.size(), 0.);
	  std::vector<float> noiseFrequencyMag(rawDigit.at(ich).fADC.size(), 0.);
	  
	  for ( unsigned int tick=0; tick < rawDigit.at(ich).fADC.size(); tick++)
	    noiseTime.at(tick) = rawDigit.at(ich).fADC.at(tick);//sin(2*3.14*0.1*tick);//

	  //do FFT
	  fFFT->DoFFT( noiseTime, noiseFrequency);

	  //get magnitude
	  for (unsigned int f=0; f < noiseFrequency.size(); f++)
	    noiseFrequencyMag.at(f) =  noiseFrequency.at(f).Rho();
	  
	  for (int tick=0; tick < fFFT->FFTSize()/2+1; tick++)
	    NoiseSpec.back()->SetBinContent( tick+1, noiseFrequencyMag.at(tick) );
	}
      }//for all channels
      
    }
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
    std::vector<float> Areas(n_channels);
    std::vector<float> minADCs(n_channels);
    std::vector<float> timeMax(n_channels);
    std::map< unsigned int, uint32_t > chanmap;

    //make map for channel index to channel number
    genChanMap( rawDigitVector, chanmap );

    //decide if empty or pulse subrun
    if (fEmptyRun){
    
      //now run the code
      analyzeEmptyEvent(rawDigitVector, pedestals, noise, noise_spectrum);
      
      fPedestalData.push_back( pedestals );
      fNoiseData.push_back( noise );
      fNoiseSpectra.push_back( noise_spectrum );
      fChanMap.push_back( chanmap );

      if ( fGetNoiseSpec ){
	//if not saving averages, make sure event number is in list of events one wants to save (fEvtNum)
	if ( ( (!fAverages) and 
	       ( (std::find(fEvtList.begin(), fEvtList.end(), fEvtNum) != fEvtList.end()) or (fEvtList.at(0) == -1) ) ) 
	     or fAverages ) {
	  GetNoiseSpec( rawDigitVector, fNoiseSpec, chanmap );
	  GetWaveforms( rawDigitVector, fWaveforms, chanmap );
	}
      }

    }//empty subrun

    if (fGainRun){ //pulse data

      //now run the code
      analyzeGainEvent(rawDigitVector, pedestals, noise,
		       maxADCs, Areas, minADCs,
		       timeMax, fprePulseTicks);

      fPedestalData.push_back( pedestals );
      fNoiseData.push_back( noise );
      fGainData.push_back( maxADCs );
      fGainAreaData.push_back ( Areas );
      fChanMap.push_back( chanmap );
      fTimeData.push_back( timeMax );
    }//gain subrun


    //if short run: generally for single-subrun run
    //plot waveform
    if (fShortRun){

      //Save waveforms to histogram
      //if ( evt.event() == 1 )
      GetWaveforms( rawDigitVector, fWaveforms, chanmap );

    }//short subrun

  }

  DEFINE_ART_MODULE(CalibrationTPC)

} //end namespace calibration

#endif //CALIBRATIONTPC_H
