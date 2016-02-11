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
#include <limits>
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

#include "lardata/Utilities/LArFFT.h"
#include "lardata/RawData/RawDigit.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/RecoBase/Wire.h"
#include "uboone/CalData/ROIAlg.h"

#include "CalibrationTPC_Algs.h"

namespace calibration {

  //  template <class Digit>

  class CalibrationTPC : public art::EDAnalyzer {

  public:

    // Configuration type: ASIC Gain and Shaping Time info
    typedef float asicGain;
    typedef float shapingT;
    typedef float voltage;
    typedef std::pair<asicGain,shapingT> Config;

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
    void getMeanAndVariance( const std::vector< std::vector<float> >& data,
			     int ch, float& mean, float& var) const;
    void prepareChannelInfoMaps();
    void prepareRunGainConfig();
    
  private:

    // ROI Alg instance
    std::unique_ptr< ::util::ROIAlg<short> > fROIAlgPtr;

    // debug mode
    int fDebug;

    //Keep Track of event number
    int fEvtNum; //Number of current event

    // An instance of CalibrationAlgs
    CalibrationAlgs _algs;

    //******************************
    //Variables Taken from FHICL File
    std::string       fRawDigitModuleLabel;   //label for rawdigit module
    uint32_t          frunNum;                //Subrun Number taken from event
    uint32_t          fsubRunNum;             //Subrun Number taken from event
    int               fsubRunCounter;         // Keep track of number of subruns passed. To figure out Configuration
    std::vector<float> fLinFitVmaxColl;
    std::vector<float> fLinFitVmaxInd;

    // A holder to keep track of the configuration for the current subrun
    Config fThisConfig;

    // Hold ASIC Gain, Shaping Time, and voltages for each subrun.
    // read in from fhicl file
    std::vector<asicGain> fSubRunASICGains;
    std::vector<shapingT> fSubRunShapingTs;
    std::vector<voltage>  fSubRunVoltageIn;
    // All possible values for the ASIC Gain. Taken from fhicl
    std::vector<asicGain> fASICGains;
    // All possible values for the Shaping Time. From fhicl
    std::vector<shapingT> fShapingTs;
    // A vector holding all possible configurations (fASICGains x gShapingTs)
    std::vector<Config> fConfigList;
    // A map holding all gain runs. For each Config store
    // a vector of all voltages used for this config set
    std::map< Config, std::vector<voltage> > fGainRuns;
    // Maps linking Configurations with vectors holding
    // amplitude and area measurements for each voltage per channel
    // each entry in the vector is a different channel
    std::map< Config, std::vector<std::vector<int  > > > fSignalRegion;
    std::map< Config, std::vector<std::vector<float> > > fAreaGainData;
    std::map< Config, std::vector<std::vector<float> > > fAmpGainData;
    std::map< Config, std::vector<std::vector<float> > > fAreaGainDataErr;
    std::map< Config, std::vector<std::vector<float> > > fAmpGainDataErr;
    //******************************

    // these are containers for the calibration results
    // Intended design: each of these is reinitialized at subrun begin
    // Thus, fPedestalData[ie][ich] = mean pedestal for event ie, channel ich
    std::vector< std::vector<float> > fPedestalData;
    std::vector< std::vector<float> > fNoiseData;
    std::vector< std::vector<float> > fPulseAmpData;
    std::vector< std::vector<float> > fPulseAreaData;
    std::vector< std::vector<float> > fTimeData;
    //Create a Map to link ordering of channel [ich] in above vector
    //to channel number.
    std::vector< std::map< unsigned int, uint32_t> > fChanMap;
    // Map that links channel number to plane number
    std::map< uint32_t, unsigned int > fChanToPlaneNumMap;
    // Map that links channel number to plane type
    std::map< uint32_t, unsigned int > fChanToPlaneTypeMap;


    //*************************************************
    //Output Histograms and MRT-Numbering 2D Maps etc..

    // Chan Data TTree:
    TTree* fDataTree;
    // Tree Variables
    int _chNum, _chIdx, _plane, _signalType;
    int _run, _subrun, _numEvents;
    float _Vin, _ASICgain, _shapingTime;
    float _pedestal, _pedestalErr;
    float _RMSnoise, _RMSnoiseErr;
    float _maxADC, _maxADCErr, _minADC, _minADCErr;
    float _areaADC, _areaADCErr;
    // Gain Tree
    TTree* fGainTree;
    // Tree Variables
    //* channel index the same
    float _asicGain;
    float _shapingT;
    float _areaConst;
    float _areaConstErr;
    float _areaGain;
    float _areaGainErr;
    float _areaChisquared;
    float _areaResidualSumSquared;
    float _areaResidualSum;
    float _ampConst;
    float _ampConstErr;
    float _ampGain;
    float _ampGainErr;
    float _ampChisquared;
    float _ampResidualSumSquared;
    float _ampResidualSum;
    // Event Tree info (once per run)
    TTree *fRunTree;
    uint64_t _timeMin;
    uint64_t _timeMax;
    int _subRunMin;
    int _subRunMax;


  }; //end class CalibrationTPC


  //-------------------------------------------------------------------
  CalibrationTPC::CalibrationTPC(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset){ 
    this->reconfigure(pset); 
    _algs.PrepareGainModel();
  }


  //-------------------------------------------------------------------
  CalibrationTPC::~CalibrationTPC(){}


  //-------------------------------------------------------------------
  void CalibrationTPC::reconfigure(fhicl::ParameterSet const& pset){

    fRawDigitModuleLabel = pset.get<std::string>("RawDigitModuleLabel");
    
    fSubRunASICGains     = pset.get< std::vector <asicGain> >("SubRunASICGains");
    fSubRunShapingTs     = pset.get< std::vector <shapingT> >("SubRunShapingTs");
    fSubRunVoltageIn     = pset.get< std::vector <voltage> > ("SubRunVoltageIn");

    fASICGains           = pset.get< std::vector<asicGain> >("ASICGains");
    fShapingTs           = pset.get< std::vector<shapingT> >("ShapingTs");

    fDebug               = pset.get<int>("Debug");

    // Now that we have read from fhicl file, prepare the gain configuration info for this run
    prepareRunGainConfig();

    // Pre-define time min & time max
    _timeMin = std::numeric_limits<uint64_t>::max();
    _timeMax = std::numeric_limits<uint64_t>::min();
    _subRunMin = std::numeric_limits<int>::max();
    _subRunMax = std::numeric_limits<int>::min();

    // Give the ROIAlg the fhicl parameters, if there are any
    std::unique_ptr< ::util::ROIAlg<short> > new_ptr(::util::ROIAlg<short>::MakeROIAlg(pset));
    fROIAlgPtr.swap(new_ptr);

    return;
  }


  //-------------------------------------------
  void CalibrationTPC::prepareChannelInfoMaps()
  {

    art::ServiceHandle<geo::Geometry> geom;    
    uint32_t nchannels = geom->Nchannels();
    unsigned int sigType = 0;
    // For each channel, make an entry in the maps
    for (uint32_t ch=0; ch < nchannels; ch++){
      
      std::vector<geo::WireID> wids = geom->ChannelToWire(ch);
      unsigned int thePlane = wids[0].Plane;
      geo::SigType_t signal = geom->SignalType(ch);

      if (signal == geo::SigType_t::kInduction)
	sigType = 2;  // 2 for by-polar pulse
      else if (signal == geo::SigType_t::kCollection) 
	sigType = 1;  // 1 for uni-polar pulse
       else
	 sigType = 0;  // 0 for unknown

       fChanToPlaneNumMap[ch]  = thePlane;
       fChanToPlaneTypeMap[ch] = sigType;
     }


     // Prepare the data-sets for Gain measurements
     // use the number of channels and the total
     // configurations to reserve the correct number
     // of entries for vectors
     fSignalRegion.clear();
     fAreaGainData.clear();
     fAmpGainData.clear();
     fAreaGainDataErr.clear();
     fAmpGainDataErr.clear();
     // generic vector of vector of floats
     // the size of the total number of channels
     std::vector<std::vector<int  > > tmpVecInt(nchannels+100);
     std::vector<std::vector<float> > tmpVec(nchannels+100);
     for (auto& config : fConfigList){
       fSignalRegion[config] = tmpVecInt;
       fAreaGainData[config] = tmpVec;
       fAmpGainData[config]  = tmpVec;
       fAreaGainDataErr[config] = tmpVec;
       fAmpGainDataErr[config]  = tmpVec;
     }

   return;
   }


   //-----------------------------------------
   void CalibrationTPC::prepareRunGainConfig()
   {

     // clear the list of configurations
     fConfigList.clear();

     // Given the ASIC gains and the shaping times
     // Make all possible combinations. Each combination
     // is a potential configuration used during the run
     std::vector<voltage> empty;
     for (auto &a : fASICGains){
       for (auto &s : fShapingTs){
	 Config tmp = std::make_pair(a,s);
	 if (fDebug) { std::cout << "Preparing Configuration: [" << a << ", " << s << "]" << std::endl; }
	 // Add configuration to list
	 fConfigList.push_back(tmp);
	 // And add to map linking to voltage values for this run
	 fGainRuns[tmp] = empty;
       }
     }

     // Now run through all the subruns, match the
     // ASIC gain and shaping times to the possible
     // configurations and add the voltage to the list
     // later we will re-order the voltages.
     if ( (fSubRunASICGains.size() != fSubRunShapingTs.size()) || (fSubRunASICGains.size() != fSubRunVoltageIn.size()) ){
       std::cout << "ERROR. These three lists should be the same length!" << std::endl;
       return;
     }


     for (size_t n=0; n < fSubRunASICGains.size(); n++){
       Config tmp( fSubRunASICGains[n] , fSubRunShapingTs[n] );
       // this configuration better exist in the list we just created!
       if ( fGainRuns.find(tmp) != fGainRuns.end() )
	 fGainRuns[tmp].push_back( fSubRunVoltageIn[n] );
     }

     if (fDebug){
       std::cout << "Configurations found for this run:" << std::endl;
       for (size_t i=0; i < fConfigList.size(); i++){
	 if ( fGainRuns.find(fConfigList.at(i)) != fGainRuns.end() ){
	   if ( fGainRuns[fConfigList.at(i)].size() != 0 ){
	     std::cout << "[" << fConfigList.at(i).first << ", " << fConfigList.at(i).second << "]  ->  [";
	     for (size_t j=0; j < fGainRuns[fConfigList.at(i)].size(); j ++)
	       std::cout << fGainRuns[fConfigList.at(i)][j] << ", ";
	     std::cout << "]" << std::endl;
	   }// if this configuration has entries
	 }// if this configuration was found
       }//for all configurations
     }//if debug

     // For each Configuration found, fGainRuns now links to
     // a list of voltages used in this run for that config.

     return;
   }


   //-------------------------------------------------------------------
     void CalibrationTPC::beginJob(){

     art::ServiceHandle<art::TFileService> tfs;
     fsubRunNum = 0;
     fsubRunCounter = 0;

     // Use the geometry to fill a map from LArSoft Channel Number
     // to Plane Type, and one from Channel number to signal type
     prepareChannelInfoMaps();

     // Setup Tree for subrun-by-subrun measurements
     fDataTree = tfs->make<TTree>("fDataTree","Per Channel Data Holder");
     // Tree Branches
     fDataTree->Branch("_chNum",&_chNum,"chNum/I");
     fDataTree->Branch("_chIdx",&_chIdx,"chIdx/I");
     fDataTree->Branch("_plane",&_plane,"plane/I");
     fDataTree->Branch("_signalType",&_signalType,"signalType/I");
     fDataTree->Branch("_run",&_run,"run/I");
     fDataTree->Branch("_subrun",&_subrun,"subrun/I");
     fDataTree->Branch("_numEvents",&_numEvents,"numEvents/I");
     fDataTree->Branch("_Vin",&_Vin,"Vin/F");
     fDataTree->Branch("_ASICgain",&_ASICgain,"ASICgain/F");
     fDataTree->Branch("_shapingTime",&_shapingTime,"shapingTime/F");
     fDataTree->Branch("_pedestal",&_pedestal,"pedestal/F");
     fDataTree->Branch("_pedestalErr",&_pedestalErr,"pedestalErr/F");
     fDataTree->Branch("_RMSnoise",&_RMSnoise,"RMSnoise/F");
     fDataTree->Branch("_RMSnoiseErr",&_RMSnoiseErr,"RMSnoiseErr/F");
     fDataTree->Branch("_maxADC",&_maxADC,"maxADC/F");
     fDataTree->Branch("_maxADCErr",&_maxADCErr,"maxADCErr/F");
     fDataTree->Branch("_areaADC",&_areaADC,"areaADC/F");
     fDataTree->Branch("_areaADCErr",&_areaADCErr,"areaADCErr/F");
     fDataTree->Branch("_minADC",&_minADC,"minADC/F");
     fDataTree->Branch("_minADCErr",&_minADCErr,"minADCErr/F");

     // Setup Tree for Gain measurements
     fGainTree = tfs->make<TTree>("fGainTree","Per Channel Gain Holder");
     // Tree Branches
     fGainTree->Branch("_chNum",&_chNum,"chNum/I");
     fGainTree->Branch("_chIdx",&_chIdx,"chIdx/I");
     fGainTree->Branch("_plane",&_plane,"plane/I");
     fGainTree->Branch("_signalType",&_signalType,"signalType/I");
     fGainTree->Branch("_asicGain",&_asicGain,"asicGain/F");
     fGainTree->Branch("_shapingT",&_shapingT,"shapingT/F");
     fGainTree->Branch("_run",&_run,"run/I");
     // Area-measurements
     fGainTree->Branch("_areaGain",&_areaGain,"areaGain/F");
     fGainTree->Branch("_areaGainErr",&_areaGainErr,"areaGainErr/F");
     fGainTree->Branch("_areaConst",&_areaConst,"areaConst/F");
     fGainTree->Branch("_areaConstErr",&_areaConstErr,"areaConstErr/F");
     fGainTree->Branch("_areaChisquared",&_areaChisquared,"areaChisquared/F");
     fGainTree->Branch("_areaResidualSum",&_areaResidualSum,"areaResidualSum/F");
     fGainTree->Branch("_areaResidualSumSquared",&_areaResidualSumSquared,"areaResidualSumSquared/F");
     // Amplitude-measurements
     fGainTree->Branch("_ampGain",&_ampGain,"ampGain/F");
     fGainTree->Branch("_ampGainErr",&_ampGainErr,"ampGainErr/F");
     fGainTree->Branch("_ampConst",&_ampConst,"ampConst/F");
     fGainTree->Branch("_ampConstErr",&_ampConstErr,"ampConstErr/F");
     fGainTree->Branch("_ampChisquared",&_ampChisquared,"ampChisquared/F");
     fGainTree->Branch("_ampResidualSum",&_ampResidualSum,"ampResidualSum/F");
     fGainTree->Branch("_ampResidualSumSquared",&_ampResidualSumSquared,"ampResidualSumSquared/F");

     // Setup Run-info TTree
     fRunTree = tfs->make<TTree>("fRunTree","Per Run Info Holder");
     fRunTree->Branch("_run",&_run,"run/I");
     fRunTree->Branch("_timeMin",&_timeMin,"timeMin/l");
     fRunTree->Branch("_timeMax",&_timeMax,"timeMax/l");
     fRunTree->Branch("_subRunMin",&_subRunMin,"subRunMin/I");
     fRunTree->Branch("_subRunMax",&_subRunMax,"subRunMax/I");
   }

   //-------------------------------------------------------------------
   void CalibrationTPC::endJob(){}

   //-------------------------------------------------------------------
   void CalibrationTPC::beginRun(art::Run const& run){

     return;
   }

   //-------------------------------------------------------------------
   void CalibrationTPC::endRun(art::Run const& run){

     // Save Run Tree info
     if (fDebug) { std::cout << "End of Run started." << std::endl; }
     fRunTree->Fill();

     // All data has been acquired. Time to calculate gains
     // and write them to their tree

     // Loop over configurations that this run took data for
     for (auto& config : fConfigList){ //size_t i=0; i < fConfigList.size(); i++){
       if (fDebug) { std::cout << "new config." << std::endl; }
       if ( fGainRuns.find(config) != fGainRuns.end() ){
	 if (fDebug) { std::cout << "config found." << std::endl; }
	 if ( fGainRuns[config].size() != 0 ){
	   if (fDebug) { std::cout << "Configuration: [" << config.first << "," << config.second << "]" << std::endl; }
	   // We found this config. Now go to data-holer map
	   // and go through channel-by-channel data for this config
	   // and calculate gain
	   for ( size_t chanIndex=0; chanIndex < fChanMap[0].size(); chanIndex++ ){
	     _chNum   = fChanMap[0].at(chanIndex);
	     _chIdx   = chanIndex;
	     // Assign plane num
	     try { _plane = fChanToPlaneNumMap.at(_chNum); } catch (const std::out_of_range& oor) { _plane = 4; }
	     // Assign plane type
	     try { _signalType = fChanToPlaneTypeMap.at(_chNum); } catch (const std::out_of_range& oor) { _signalType = 0; }
	     //for (size_t ch=0; ch < fAreaGainData[config].size(); ch++){
	     // Do not proceed if voltages, areas, errors not all same size
	     // and not all at least 2 elements long
	     if ( (fAreaGainDataErr[config][_chIdx].size() < 2) ||
		  (fGainRuns[config].size() != fAreaGainData[config][_chIdx].size()) ||
		  (fGainRuns[config].size() != fAreaGainDataErr[config][_chIdx].size()) ){
	       /*
	       if (fDebug)
		 std::cout << "Not enough voltages / sizes do not match for Config: [" 
			   << config.first << ", " << config.second << "] and Channel: "
			   << ch << std::endl; 
	       */
	       continue;
	     }

	     // Voltages: fGainRuns[config]
	     // Amplitude/Area and errors in: fAreaGainData[config][ch] and similar
	     // so...calculate gain!
	     if (fDebug){
	       std::cout << "Channel: " << _chIdx << std::endl << "Voltages: [";
	       for (auto& v : fGainRuns[config])
		 std::cout << v << ",";
	       std::cout << "]\tArea data: [";
	       for (auto& d : fAreaGainData[config][_chIdx])
		 std::cout << d << ",";
	       std::cout << "]\tArea errors: [";
	       for (auto& e : fAreaGainDataErr[config][_chIdx])
		 std::cout << e << ",";
	       std::cout << "]" << std::endl;
	     }
	     _algs.calcGain( fGainRuns[config], fAreaGainData[config][_chIdx], fAreaGainDataErr[config][_chIdx],
			     _areaChisquared, _areaGain, _areaGainErr, _areaConst, _areaConstErr,
			     _areaResidualSum, _areaResidualSumSquared);
	     _algs.calcGain( fGainRuns[config], fAmpGainData[config][_chIdx], fAmpGainDataErr[config][_chIdx],
			     _ampChisquared, _ampGain, _ampGainErr, _ampConst, _ampConstErr,
			     _ampResidualSum, _ampResidualSumSquared);

	     _asicGain = config.first;
	     _shapingT = config.second;
	     _run      = frunNum;
	     
	     fGainTree->Fill();

	     if (fDebug)
	       std::cout << "Area Gain: " << _areaGain << "\tErr: " << _areaGainErr << "\tConst: " << _areaConst
			 << "\tConst Err: " << _areaConstErr << "\tChi: " << _areaChisquared << std::endl << std::endl;
	   }// for all channels
	 }// if this configuration was used in this subrun
       }// if this configuration exists
     }// loop over configurations

     return;
   }


   //-------------------------------------------------------------------
   void CalibrationTPC::beginSubRun(art::SubRun const& subrun){

     // Get the current Configuration from the subrun number
     fThisConfig = std::make_pair( fSubRunASICGains.at(fsubRunCounter), fSubRunShapingTs.at(fsubRunCounter) );

     // Clear all subrun-by-subrun data holders
     fPedestalData.clear();
     fNoiseData.clear();
     fPulseAmpData.clear();
     fPulseAreaData.clear();
     fTimeData.clear();
     fChanMap.clear();
     fTimeData.clear();

     return;
   }

   //-------------------------------------------------------------------
   void CalibrationTPC::endSubRun(art::SubRun const& subrun){

    if (fDebug) { std::cout << "subrun end started." << std::endl; }

     // Now that the subrun has ended, take all parameters that were
     // calculated and fill the tree with one entry per channel

     // Make sure the channel map is identical for all events
     // if not maybe something went wrong when reading out the data
     // or swizzling?
     // ************************  TO DO  **************************

     // Get sub-run specific quantities (same for all channels)
     _numEvents = fChanMap.size();
     _subrun  = fsubRunNum;
     if (_subrun < _subRunMin) { _subRunMin = _subrun; }
     if (_subrun > _subRunMax) { _subRunMax = _subrun; }
     fsubRunCounter += 1;
     _Vin         = fSubRunVoltageIn[_subrun-1];
     _ASICgain    = fSubRunASICGains[_subrun-1];
     _shapingTime = fSubRunShapingTs[_subrun-1];
     // Loop over all channels
     if (fDebug) { std::cout << "Loop over channels." << std::endl; }
     for ( size_t chanIndex=0; chanIndex < fChanMap[0].size(); chanIndex++ ){
       _chNum   = fChanMap[0].at(chanIndex);
       _chIdx   = chanIndex;
       // Assign plane num
       try { _plane = fChanToPlaneNumMap.at(_chNum); } catch (const std::out_of_range& oor) { _plane = 4; }
       // Assign plane type
       try { _signalType = fChanToPlaneTypeMap.at(_chNum); } catch (const std::out_of_range& oor) { _signalType = 0; }
       // Now for each quantity measured save the mean and standard deviation
       // across all events and save it to the tree
       getMeanAndVariance( fPedestalData, _chIdx, _pedestal, _pedestalErr  );
       getMeanAndVariance( fNoiseData,    _chIdx, _RMSnoise, _RMSnoiseErr  );
       getMeanAndVariance( fPulseAmpData,  _chIdx, _maxADC,   _maxADCErr   );
       getMeanAndVariance( fPulseAreaData,  _chIdx, _areaADC,   _areaADCErr);
       // Fill Gain Data maps
       if ( fAreaGainData.find(fThisConfig) == fAreaGainData.end() ||
	    fAmpGainData.find(fThisConfig) == fAreaGainData.end() )
	 std::cerr << "Error in CalibrationTPC_module::endSubRun" 
		   << "\nCurrent configuration not found in configuration list";
       fAreaGainData[fThisConfig][_chNum].push_back(_areaADC);
       fAreaGainDataErr[fThisConfig][_chNum].push_back(_areaADCErr);
       fAmpGainData[fThisConfig][_chNum].push_back(_maxADC);
       fAmpGainDataErr[fThisConfig][_chNum].push_back(_maxADCErr);
       // Finally, fill tree per-event
       fDataTree->Fill();
     }

     if (fDebug) { std::cout << "subrun end ended." << std::endl; }

     return;
   }


   //------------------------------------------------------------------------------------
   void CalibrationTPC::getMeanAndVariance( const std::vector< std::vector<float> >& data,
					    int ch, float& mean, float& var) const
   {
     // in data vector first index is event number, second is channel number
     int events = data.size();
     mean = 0.;
     var  = 0.;
     for (size_t i=0; i < data.size(); i++)
       mean += data.at(i).at(ch);
     mean /= events;
     
    for (size_t i=0; i < data.size(); i++)
      var += (data.at(i).at(ch)-mean)*(data.at(i).at(ch)-mean);
    var = sqrt ( var / events );

    return;
  }

  //--------------------------------------------------
  void CalibrationTPC::analyze(art::Event const& evt){

    fEvtNum = evt.event();

    frunNum    = evt.run();
    fsubRunNum = evt.subRun();

    // Get event TimeStamp
    art::Timestamp evtTime = evt.time();
    uint64_t evttime = evtTime.value();
    if (evttime < _timeMin) { _timeMin = evttime; }
    if (evttime > _timeMax) { _timeMax = evttime; }

    
    LOG_INFO ("CalibrationTPC")
      << "Processing Run " << frunNum
      << ", Subrun " << fsubRunNum 
      << ", Event " << evt.event();

    art::Handle< std::vector<raw::RawDigit> > rawDigitHandle;
    evt.getByLabel(fRawDigitModuleLabel,rawDigitHandle);
    std::vector<raw::RawDigit> const& rawDigitVector(*rawDigitHandle);

    
    //if(hasCompressedRawDigit(rawDigitVector))
    //  throw "ERORR! You can't run the CalibrationTPC analyzer with compressed rawDigits!";

   //initialize per-event vectors
    const size_t n_channels = rawDigitVector.size();
    std::vector<float> pedestals(n_channels,0);
    std::vector<float> noise(n_channels,0);
    std::vector<float> maxADCs(n_channels,0);
    std::vector<float> Areas  (n_channels,0);
    std::vector<float> minADCs(n_channels,0);
    std::vector<float>   timeMax(n_channels,0);
    std::map< unsigned int, uint32_t > chanmap;

    //make map for channel index to channel number
    //if (fDebug) { std::cout << "Generating channel map" << std::endl; }
    _algs.genChanMap( rawDigitVector, chanmap );

    // Loop over all channels
    //if (fDebug) { std::cout << "About to loop over all channels" << std::endl; }
    for (size_t i=0; i < rawDigitVector.size(); i++){

      auto wf = rawDigitVector.at(i);
      // Run ROIAlg to find quiet and pulse regions
      //if (fDebug) { std::cout << "About to process waveform of size " << wf.fADC.size() << std::endl; }
      fROIAlgPtr->ProcessWaveform(wf.ADCs());
      // now we can access the baseline regions
      // to get noise and baseline measurements for this channel
      _algs.calcPedestal_BaselineRegion(fROIAlgPtr->GetBaselineRegions(),pedestals.at(i));
      _algs.calcNoise_BaselineRegion(fROIAlgPtr->GetBaselineRegions(),pedestals.at(i),noise.at(i));
      // And similarly for the Signal Region
      // If more than one Signal region found, stop!
      // We are trying to measure gain from a single pulse
      // if multiple are found something is not right!
      if (fROIAlgPtr->GetNSignalRegions() == 0){
	//if (fDebug) { std::cout << "No signal regions found." << std::endl; }
	maxADCs.at(i) = pedestals.at(i);
      }
      // only look at first ROI region -> leading edge of pulse
      else{
	//if (fDebug) { std::cout << "No signal regions found." << std::endl; }
	_algs.calcSignal_SignalRegion(fROIAlgPtr->GetSignalRegions(), pedestals.at(i),
				    maxADCs.at(i), Areas.at(i), minADCs.at(i), timeMax.at(i));
      }

    }

    // Add values to event vectors
    fPedestalData.push_back( pedestals );
    fNoiseData.push_back( noise );
    fPulseAmpData.push_back( maxADCs );
    fPulseAreaData.push_back ( Areas );
    fChanMap.push_back( chanmap );
    fTimeData.push_back( timeMax );

    if (fDebug) { std::cout << "Event has ended." << std::endl; }
    
    return;
  }

  DEFINE_ART_MODULE(CalibrationTPC)

} //end namespace calibration

#endif //CALIBRATIONTPC_H
