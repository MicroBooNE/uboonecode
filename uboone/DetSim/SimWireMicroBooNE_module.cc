////////////////////////////////////////////////////////////////////////
// $Id: SimWireMicroBooNE.cxx,v 1.22 2010/04/23 20:30:53 seligman Exp $
//
// SimWireMicroBooNE class designed to simulate signal on a wire in the TPC
//
// katori@fnal.gov
//
// - Revised to use sim::RawDigit instead of rawdata::RawDigit, and to
// - save the electron clusters associated with each digit.
//
////////////////////////////////////////////////////////////////////////

// C/C++ standard library
#include <stdexcept> // std::range_error
#include <vector>
#include <string>
#include <algorithm> // std::fill()
#include <functional>
#include <fstream>
#include <random>
#include <chrono>

// CLHEP libraries
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

// ROOT libraries
#include "TMath.h"
#include "TComplex.h"
#include "TString.h"
#include "TH2.h"
#include "TH1D.h"
#include "TFile.h"

// art library and utilities
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art extensions
#include "larsim/RandomUtils/LArSeedService.h"

// LArSoft libraries
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h" // FIXME: this is not portable
#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"
#include "lardataobj/Simulation/sim.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"


///Detector simulation of raw signals on wires
namespace detsim {

  // Base class for creation of raw signals on wires.
  class SimWireMicroBooNE : public art::EDProducer {

  public:

    explicit SimWireMicroBooNE(fhicl::ParameterSet const& pset);
    virtual ~SimWireMicroBooNE();

    // read/write access to event
    void produce (art::Event& evt);
    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& p);

  private:

    void GenNoiseInTime(std::vector<float> &noise, double noise_factor) const;
    void GenNoiseInFreq(std::vector<float> &noise, double noise_factor) const;
    void GenNoisePostFilter(std::vector<float> &noise, double noise_factor, size_t view, int chan);
    void MakeADCVec(std::vector<short>& adc, std::vector<float> const& noise, 
                    std::vector<double> const& charge, float ped_mean) const;

    std::string             fDriftEModuleLabel; ///< module making the ionization electrons
    raw::Compress_t         fCompression;       ///< compression type to use

    double                  fNoiseWidth;        ///< exponential noise width (kHz)
    double                  fNoiseRand;         ///< fraction of random "wiggle" in noise in freq. spectrum
    double                  fLowCutoff;         ///< low frequency filter cutoff (kHz)
    double                  fSampleRate;        ///< sampling rate in ns

    size_t                  fNTicks;	        ///< number of ticks of the clock    
    unsigned int            fNTimeSamples;      ///< number of ADC readout samples in all readout frames (per event)	 	   

    std::vector<TH1D*>      fNoiseDist;     ///< distribution of noise counts, one per plane
    bool                    fGetNoiseFromHisto; ///< if True -> Noise from Histogram of Freq. spectrum
    unsigned short          fGenNoise;          ///< 0 -> no noise, 1: time domain, 2: freq domain, 3: postfilter
    std::string             fNoiseFileFname;
    std::string             fNoiseHistoName;
    TH1D*                   fNoiseHist;         ///< distribution of noise counts

    std::map< double, int > fShapingTimeOrder;
    std::string             fTrigModName;       ///< Trigger data product producer name
    
    bool                    fSimDeadChannels;   ///< if True, simulate dead channels using the ChannelStatus service.  If false, do not simulate dead channels

    bool fMakeNoiseDists;

    bool        fTest; // for forcing a test case
    std::vector<sim::SimChannel> fTestSimChannel_v;
    size_t      fTestWire;
    std::vector<size_t> fTestIndex;
    std::vector<double> fTestCharge;

    int         fSample; // for histograms, -1 means no histos

    //define max ADC value - if one wishes this can
    //be made a fcl parameter but not likely to ever change
    const float adcsaturation = 4095;

    ::detinfo::ElecClock fClock; ///< TPC electronics clock

    TH1D* hTest[5] = {0, 0, 0, 0, 0};

    // little helper class to hold the params of each charge dep
    class ResponseParams {
    public:
      ResponseParams(double charge, size_t time) : m_charge(charge), m_time(time) {}
      double getCharge() { return m_charge; }
      size_t getTime()   { return m_time; }
    private:
      double m_charge;
      size_t m_time;
    };

    //
    // Needed for post-filter noise (pfn) generator by Jyoti 
    //
    std::vector<double> _pfn_shaping_time_v;
    TF1  *_pfn_f1;
    TF1  *_pfn_MyPoisson;
    TVirtualFFT *_pfn_ifft;
    std::vector<double> _pfn_rho_v;
    std::vector<double> _pfn_value_re;
    std::vector<double> _pfn_value_im;

    //
    // Needed for calculating Noise baseline
    //
    
    //std::vector<double> rms;

  }; // class SimWireMicroBooNE

  namespace {
    size_t _ch = 0;
    size_t _wr = 0;
  }

  DEFINE_ART_MODULE(SimWireMicroBooNE)

  //-------------------------------------------------
  SimWireMicroBooNE::SimWireMicroBooNE(fhicl::ParameterSet const& pset)
  : fNoiseHist(0)
    , _pfn_shaping_time_v()
    , _pfn_f1(nullptr)
    , _pfn_MyPoisson(nullptr)
    , _pfn_ifft(nullptr)
    , _pfn_rho_v()
    , _pfn_value_re()
    , _pfn_value_im()
  {
    this->reconfigure(pset);

    produces< std::vector<raw::RawDigit>   >();

    fCompression = raw::kNone;
    TString compression(pset.get< std::string >("CompressionType"));
    if(compression.Contains("Huffman",TString::kIgnoreCase)) fCompression = raw::kHuffman;

    // create a default random engine; obtain the random seed from LArSeedService,
    // unless overridden in configuration with key "Seed" and "SeedPedestal"
    art::ServiceHandle<sim::LArSeedService> Seeds;
    Seeds->createEngine(*this, "HepJamesRandom", "noise", pset, "Seed");
    Seeds->createEngine(*this, "HepJamesRandom", "pedestal", pset, "SeedPedestal");

  }

  //-------------------------------------------------
  SimWireMicroBooNE::~SimWireMicroBooNE()
  {
    delete fNoiseHist;
  }

  //-------------------------------------------------
  void SimWireMicroBooNE::reconfigure(fhicl::ParameterSet const& p)
  {
    fDriftEModuleLabel= p.get< std::string         >("DriftEModuleLabel");
    fNoiseWidth       = p.get< double              >("NoiseWidth");
    fNoiseRand        = p.get< double              >("NoiseRand");
    fLowCutoff        = p.get< double              >("LowCutoff");
    fGetNoiseFromHisto= p.get< bool                >("GetNoiseFromHisto");
    fGenNoise         = p.get< unsigned short      >("GenNoise");
    fSimDeadChannels  = p.get< bool                >("SimDeadChannels");

    fMakeNoiseDists   = p.get< bool                >("MakeNoiseDists", false);

    fTrigModName      = p.get< std::string         >("TrigModName");
    fTest             = p.get<bool                 >("Test");
    fTestWire         = p.get< size_t              >("TestWire");
    fTestIndex        = p.get< std::vector<size_t> >("TestIndex");
    fTestCharge       = p.get< std::vector<double> >("TestCharge");
    if(fTestIndex.size() != fTestCharge.size())
      throw cet::exception(__FUNCTION__)<<"# test pulse mismatched: check TestIndex and TestCharge fcl parameters...";
    fSample           = p.get<int                  >("Sample");


    //Map the Shaping Times to the entry position for the noise ADC
    //level in fNoiseFactInd and fNoiseFactColl
    fShapingTimeOrder = { {0.5, 0}, {1.0, 1}, {2.0, 2}, {3.0, 3} };

    if(fGetNoiseFromHisto)
    {
      fNoiseHistoName= p.get< std::string         >("NoiseHistoName");

      cet::search_path sp("FW_SEARCH_PATH");
      sp.find_file(p.get<std::string>("NoiseFileFname"), fNoiseFileFname);

      TFile * in=new TFile(fNoiseFileFname.c_str(),"READ");
      TH1D * temp=(TH1D *)in->Get(fNoiseHistoName.c_str());

      if(temp!=NULL)
      {
        fNoiseHist=new TH1D(fNoiseHistoName.c_str(),fNoiseHistoName.c_str(),temp->GetNbinsX(),0,temp->GetNbinsX());
        temp->Copy(*fNoiseHist);
      }
      else
        throw cet::exception("SimWireMicroBooNE") << "Could not find noise histogram in Root file\n";
      in->Close();

    }
    //detector properties information
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fSampleRate    = detprop->SamplingRate();
    fNTimeSamples  = detprop->NumberTimeSamples();

    // make the histos if not already made
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;

    if(hTest[0] == 0) {
      char buffer[80];
      //char buffer1[80];

      if(fSample>=0) {
        for(size_t i=0;i<5;++i) {
          sprintf(buffer, "hTest%i",(int)i);
          hTest[i] = tfs->make<TH1D>(buffer, buffer, 500, -250., 250.);
        }
      }
    }

    return;
  }

  //-------------------------------------------------
  void SimWireMicroBooNE::beginJob()
  {

    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
    
    char buff0[80], buff1[80];

    if(fMakeNoiseDists) {
      fNoiseDist.resize(3,0);
      for(int view = 0; view<3; ++view) {
        sprintf(buff0, "Noise%i", view);
        sprintf(buff1, ";Noise on Plane %i(ADC);", view);
        fNoiseDist[view]  = tfs->make<TH1D>(buff0, buff1, 1000,   -30., 30.);
      }
    }

    if(fTest){
      art::ServiceHandle<geo::Geometry> geo;  
      if(geo->Nchannels()<=fTestWire)
        throw cet::exception(__FUNCTION__)<<"Invalid test wire channel: "<<fTestWire;

      std::vector<unsigned int> channels;
      channels.reserve(geo->PlaneIDs().size());
      for(auto const& plane_id : geo->PlaneIDs())
        channels.push_back(geo->PlaneWireToChannel(plane_id.Plane,fTestWire));

      double xyz[3] = { std::numeric_limits<double>::max() };
      for(auto const& ch : channels) {

        fTestSimChannel_v.push_back(sim::SimChannel(ch));

        for(size_t i=0; i<fTestIndex.size(); ++i){

          fTestSimChannel_v.back().AddIonizationElectrons(-1,
                                                          fTestIndex[i],
                                                          fTestCharge[i],
                                                          xyz,
                                                          std::numeric_limits<double>::max());
        }
      }
    }

    return;

  }

  //-------------------------------------------------
  void SimWireMicroBooNE::endJob()
  {
  
  
  }

  void SimWireMicroBooNE::produce(art::Event& evt)
  {

    
    //--------------------------------------------------------------------
    //
    // Get all of the services we will be using
    //
    //--------------------------------------------------------------------
    
    //get pedestal conditions
    const lariov::DetPedestalProvider& pedestalRetrievalAlg 
       = art::ServiceHandle<lariov::DetPedestalService>()->GetPedestalProvider();
    
    //get rng for pedestals
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine("pedestal");   
    
    //channel status for simulating dead channels
    const lariov::ChannelStatusProvider& ChannelStatusProvider
       = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();

    //get the FFT
    art::ServiceHandle<util::LArFFT> fFFT;
    fFFT->ReinitializeFFT(fNTimeSamples,fFFT->FFTOptions(),fFFT->FFTFitBins());
    fNTicks = fFFT->FFTSize();

    if ( fNTicks%2 != 0 )
      LOG_DEBUG("SimWireMicroBooNE") << "Warning: FFTSize " << fNTicks << " not a power of 2. "
      << "May cause issues in (de)convolution.\n";

    if ( fNTimeSamples > fNTicks )
      mf::LogError("SimWireMicroBooNE") << "Cannot have number of readout samples "
      << fNTimeSamples << " greater than FFTSize " << fNTicks << "!";

    // TFileService
    art::ServiceHandle<art::TFileService> tfs;

    //TimeService
    art::ServiceHandle<detinfo::DetectorClocksServiceStandard> tss;
    // In case trigger simulation is run in the same job...
    tss->preProcessEvent(evt);
    auto const* ts = tss->provider();

    // Check if trigger data product exists or not. If not, throw a warning
    art::Handle< std::vector<raw::Trigger> > trig_array;
    evt.getByLabel(fTrigModName, trig_array);
    if(!trig_array.isValid())

      std::cout << std::endl << "  "
      << "\033[95m" << "<<" << __PRETTY_FUNCTION__ << ">>" << "\033[00m"
      << std::endl << "  "
      << "\033[93m"
      << " No trigger data exists => will use the default trigger time set in TimeService..."
      << "\033[00m"
      << std::endl;

    // get the geometry to be able to figure out signal types and chan -> plane mappings
    art::ServiceHandle<geo::Geometry> geo;
    const size_t N_CHANNELS = geo->Nchannels();

    //Get N_RESPONSES from SignalShapingService, on the fly
    art::ServiceHandle<util::SignalShapingServiceMicroBooNE> sss;
    const std::vector<std::vector<size_t> > N_RESPONSES = sss->GetNActiveResponses();
    const size_t N_VIEWS = N_RESPONSES[0].size();




    //--------------------------------------------------------------------
    //
    // Get the SimChannels, which we will use to produce RawDigits
    //
    //--------------------------------------------------------------------

    // make a vector of const sim::SimChannel* that has same number
    // of entries as the number of channels in the detector
    // and set the entries for the channels that have signal on them
    // using the chanHandle
    std::vector<const sim::SimChannel*> channels(N_CHANNELS,nullptr);
    if(!fTest){
      std::vector<const sim::SimChannel*> chanHandle;
      evt.getView(fDriftEModuleLabel,chanHandle);
      for(size_t c = 0; c < chanHandle.size(); ++c){
        channels.at(chanHandle.at(c)->Channel()) = chanHandle.at(c);
      }
    }else{
      for(size_t c = 0; c<fTestSimChannel_v.size(); ++c)
        channels.at(fTestSimChannel_v[c].Channel()) = &(fTestSimChannel_v[c]);
    }

    // make a unique_ptr of sim::SimDigits that allows ownership of the produced
    // digits to be transferred to the art::Event after the put statement below
    std::unique_ptr< std::vector<raw::RawDigit>> digcol(new std::vector<raw::RawDigit>);
    digcol->reserve(N_CHANNELS);
   
    std::vector<std::vector<std::vector<std::unique_ptr<ResponseParams> > > > responseParamsVec(N_CHANNELS);

    _wr = 0;
    _ch = 0;
    for (auto& channel : responseParamsVec) {
      size_t view = (size_t)geo->View(_ch);
      channel.resize(2*N_RESPONSES[0][view]-1);
    }



    //--------------------------------------------------------------------
    //
    // I'm not sure about the purpose of this first for-loop: experts please update this comment!
    //
    //--------------------------------------------------------------------

    //LOOP OVER ALL CHANNELS
    // In this version we assume that adjacent channels <-> adjacent wires, in the same plane/view
    // Is this always true?
    std::vector<int> first_channel_in_view(N_VIEWS,-1);
    for(unsigned int chan = 0; chan < N_CHANNELS; chan++) {
      auto wid = geo->ChannelToWire(chan);
      size_t view = (size_t)geo->View(chan);
      
      if (first_channel_in_view[view] == -1) {
        first_channel_in_view[view] = chan;
      }

      // get the sim::SimChannel for this channel
      const sim::SimChannel* sc = channels.at(chan);
      if( !sc ) continue;

      // remove the time offset
      int time_offset = 0;//sss->FieldResponseTOffset(chan);

      // loop over the tdcs and grab the number of electrons for each
      for(int t = 0; t < (int)fNTicks; ++t) {

        int tdc = ts->TPCTick2TDC(t);
        // continue if tdc < 0
        if( tdc < 0 ) continue;
        double charge = sc->Charge(tdc);
        if(charge==0) continue;

        // Apply artificial time offset to take care of field response convolution
        // wrap the negative times to the end of the buffer
        // The offset should be taken care of in the shaping service, by shifting the response.

        int raw_digit_index =
        ( (t + time_offset) >= 0 ? t+time_offset : (fNTicks + (t+time_offset)) );

        if(raw_digit_index <= 0 || raw_digit_index >= (int)fNTicks) continue;

        // here fill ResponseParams... all the wires!
        for(int wire = -(N_RESPONSES[0][view]-1); wire<(int)N_RESPONSES[0][view]; ++wire) {
          auto wireIndex = (size_t)wire+N_RESPONSES[0][view] - 1;
          int wireChan = (int)chan + wire;
          if(wireChan<0 || wireChan>= (int)N_CHANNELS) continue;
          if ((size_t)geo->View(wireChan)!=view) continue;

          responseParamsVec[wireChan][wireIndex].emplace_back(new ResponseParams(charge, raw_digit_index));
        } // loop over wires
      } // loop over tdcs
    } // loop over channels
    


    //--------------------------------------------------------------------
    //
    // Make the noise vector that will be used for generating noise
    //
    //--------------------------------------------------------------------   
//    DoubleVec             noiseFactVec(N_VIEWS,0.);
//    auto tempNoiseVec = sss->GetNoiseFactVec();
//    for (size_t v = 0; v != N_VIEWS; ++v) {
//
//      //the current sss only allows retrieval by channel, even though these things only change by view
//      //If these ever do change by channel, then this code automatically becomes incorrect!
//      double shapingTime = sss->GetShapingTime(first_channel_in_view[v]);
//      double asicGain    = sss->GetASICGain(first_channel_in_view[v]);
//
//      if (fShapingTimeOrder.find( shapingTime ) != fShapingTimeOrder.end() ) {
//        noiseFactVec[v]  = tempNoiseVec[v].at( fShapingTimeOrder.find( shapingTime )->second );
//	      noiseFactVec[v] *= asicGain/4.7;
//      }
//      else {//Throw exception...
//        throw cet::exception("SimWireMicroBooNE")
//        << "\033[93m"
//        << "Shaping Time received from signalservices_microboone.fcl is not one of allowed values"
//        << std::endl
//        << "Allowed values: 0.5, 1.0, 2.0, 3.0 usec"
//        << "\033[00m"
//        << std::endl;
//      }
//    }

    DoubleVec             noiseFactVec(N_VIEWS,0.);
    auto tempNoiseVec = sss->GetNoiseFactVec();	
    for (size_t v = 0; v != N_VIEWS; ++v) {
      	
      //the current sss only allows retrieval by channel, even though these things only change by view
      //If these ever do change by channel, then this code automatically becomes incorrect!
      double shapingTime = sss->GetShapingTime(first_channel_in_view[v]);	
      double asicGain    = sss->GetASICGain(first_channel_in_view[v]);
      
      if (fShapingTimeOrder.find( shapingTime ) != fShapingTimeOrder.end() ) {
        if(_pfn_shaping_time_v.size()<=v) _pfn_shaping_time_v.resize(v+1,-1);
        if(_pfn_shaping_time_v[v]<0) _pfn_shaping_time_v[v]=shapingTime;
        noiseFactVec[v]  = tempNoiseVec[v].at( fShapingTimeOrder.find( shapingTime )->second );
        noiseFactVec[v] *= asicGain/4.7;
      }
      else {//Throw exception...
        throw cet::exception("SimWireMicroBooNE")
        << "\033[93m"
        << "Shaping Time received from signalservices_microboone.fcl is not one of allowed values"
        << std::endl
        << "Allowed values: 0.5, 1.0, 2.0, 3.0 usec"
        << "\033[00m"
        << std::endl;
      }
    }

    
    //--------------------------------------------------------------------
    //
    // Loop over channels a second time and produce the RawDigits by adding together 
    // pedestal, noise, and direct&induced charges
    //
    //-------------------------------------------------------------------- 
       
    // vectors for working in the following for loop
    std::vector<short>    adcvec(fNTimeSamples, 0);
    std::vector<double>   chargeWork(fNTicks,0.);
    std::vector<double>   tempWork(fNTicks,0.);
    std::vector<float>    noisetmp(fNTicks,0.);

    int step = 0;

    // various constants: not fcl-configurable
    double slope0[5] = { 0., 2.1575, 6.4725 , 13.946, 40.857};
    double t0[5] =     { 4450., 6107., 6170., 6305., 6695. };
    double wire0[3] =  { 337., 332., -0.7 };
    double factor[3] = { 2.0, 2.0, 1.0 };
    int tickCut = 250;

    // loop over the collected responses
    //   this is needed because hits generate responses on adjacent wires!
    for(unsigned int chan = 0; chan < N_CHANNELS; chan++) {


      //clean up working vectors from previous iteration of loop
      adcvec.resize(fNTimeSamples); //compression may have changed the size of this vector
      noisetmp.resize(fNTicks); //just in case
      std::fill(adcvec.begin(),     adcvec.end(),     0);
      std::fill(chargeWork.begin(), chargeWork.end(), 0.);
      std::fill(tempWork.begin(),   tempWork.end(),   0.);
      std::fill(noisetmp.begin(),   noisetmp.end(),   0.);

      // make sure chargeWork is correct size
      if (chargeWork.size() < fNTimeSamples)
        throw std::range_error("SimWireMicroBooNE: chargeWork vector too small");
	
      //use channel number to set some useful numbers
      auto wid = geo->ChannelToWire(chan);
      size_t wireNum = wid[0].Wire;
      auto& thisChan = responseParamsVec[chan];
      size_t view = (size_t)geo->View(chan);
      
      
      //Get pedestal with random gaussian variation
      CLHEP::RandGaussQ rGaussPed(engine, 0.0, pedestalRetrievalAlg.PedRms(chan));
      float ped_mean = pedestalRetrievalAlg.PedMean(chan) + rGaussPed.fire(); 

      //Generate Noise
      double             noise_factor;
      auto tempNoiseVec = sss->GetNoiseFactVec();
      double shapingTime = sss->GetShapingTime(chan);
      double asicGain    = sss->GetASICGain(chan);
      //nanostd::cout << "Sim params: " << chan << " " << shapingTime << " " << asicGain << std::endl;

      if (fShapingTimeOrder.find( shapingTime ) != fShapingTimeOrder.end() ) {
        noise_factor  = tempNoiseVec[view].at( fShapingTimeOrder.find( shapingTime )->second );
        noise_factor *= asicGain/4.7;
      }
      else {//Throw exception...
        throw cet::exception("SimWireMicroBooNE")
        << "\033[93m"
        << "Shaping Time received from signalservices_microboone.fcl is not one of allowed values"
        << std::endl
        << "Allowed values: 0.5, 1.0, 2.0, 3.0 usec"
        << "\033[00m"
        << std::endl;
      }

      if (fGenNoise){
        if (fGenNoise==1)
          GenNoiseInTime(noisetmp, noise_factor);
        else if(fGenNoise==2)
          GenNoiseInFreq(noisetmp, noise_factor);
       	else if(fGenNoise==3)
	         GenNoisePostFilter(noisetmp, noise_factor, view, chan);
      }

      //Add Noise to NoiseDist Histogram
      //geo::SigType_t sigtype = geo->SignalType(chan);
      geo::View_t vw = geo->View(chan);
      if(fMakeNoiseDists) {
        for (size_t i=step; i < fNTimeSamples; i+=1000) {
          fNoiseDist[vw]->Fill(noisetmp[i]);
        }
      }
      ++step;

      //If the channel is bad, we can stop here
      //if you are using the UbooneChannelStatusService, then this removes disconnected, "dead", and "low noise" channels
      if (fSimDeadChannels && (ChannelStatusProvider.IsBad(chan) || !ChannelStatusProvider.IsPresent(chan)) ) {
        MakeADCVec(adcvec, noisetmp, chargeWork, ped_mean);
        raw::RawDigit rd(chan, fNTimeSamples, adcvec, fCompression);
        rd.SetPedestal(ped_mean);
        digcol->push_back(std::move(rd));



        continue;
      }
      
      
      //Channel is good, so fill the chargeWork vector with charges
      int tick0 = 0;
      if(fSample>=0) tick0 = t0[fSample] - factor[view]*slope0[fSample]*(wireNum-wire0[view]) + 0.5;

      for(int wire=-(N_RESPONSES[0][view]-1); wire<(int)N_RESPONSES[0][view];++wire) {
        int wireChan = chan + wire;
        if(wireChan<0 || wireChan>= (int)N_CHANNELS) continue;
        if ((size_t)geo->View(wireChan)!=view) continue;
        size_t wireIndex = (size_t)(wire + (int)N_RESPONSES[0][view] - 1);
        auto & thisWire = thisChan[wireIndex];
        if(thisWire.empty()) continue;
        std::fill(tempWork.begin(), tempWork.end(), 0.);
        for(auto& item : thisWire) {
          auto charge = item->getCharge();
          if(charge==0) continue;
          auto raw_digit_index = item->getTime();
          if(raw_digit_index > 0 && raw_digit_index < fNTicks) {
            tempWork.at(raw_digit_index) = charge;
          }
        }

        // now we have the tempWork for the adjacent wire of interest
        // convolve it with the appropriate response function
        sss->Convolute(chan, fabs(wire), tempWork);

        // this is to generate some plots
        if(view==1 && wireNum==360 && fSample>=0) {
          if(abs(wire)>2) continue;
          size_t index = wire + 2;

          bool printWF = false;
          if(printWF)std::cout << "printout of waveform, index = " << index << std::endl;
          for(int i=tick0-tickCut; i<tick0+tickCut;++i) {
            double val = tempWork[i];
            if(printWF) {
              if((i+1)%10==0) std::cout << std::endl << i << " " << i-tick0 << " ";
              std::cout << val << " " ;
            }
            hTest[index]->Fill(i*1.-tick0, val);
          }
          if(printWF) std::cout << std::endl;

        }

        // now add the result into the "charge" vector
        for(size_t bin = 0; bin < fNTicks; ++bin) {
          chargeWork[bin] += tempWork[bin];
        }
        // or:
        //std::transform(chargeWork.begin(), chargeWork.end(), tempWork.begin(),
        //               chargeWork.begin(), std::plus<double>());

      }//end loop over response wires


      // add this digit to the collection;
      // adcvec is copied, not moved: in case of compression, adcvec will show
      // less data: e.g. if the uncompressed adcvec has 9600 items, after
      // compression it will have maybe 5000, but the memory of the other 4600
      // is still there, although unused; a copy of adcvec will instead have
      // only 5000 items. All 9600 items of adcvec will be recovered for free
      // and used on the next loop.
      MakeADCVec(adcvec, noisetmp, chargeWork, ped_mean);
      raw::RawDigit rd(chan, fNTimeSamples, adcvec, fCompression);
      rd.SetPedestal(ped_mean);
      digcol->push_back(std::move(rd)); // we do move the raw digit copy, though

    }// end of 2nd loop over channels


    // Used in calculation of baseline noise value, but hardcoded so now uneccesary.
    /* 
    double total = 0;
    for (size_t i=3499; i < 4500; i++)
     {
       std::cout << "RMS of waveform " << i << " is " << rms[i] << std::endl;
       total = total+rms[i];
     }

      double baselinerms = total/1000;

      std::cout << "\nmn\n\n\n\n\n The average RMS is:" << baselinerms << "\n\n\n\n\n\n\n" << std::endl;
    
    */
    evt.put(std::move(digcol));
    return;
  }


  //-------------------------------------------------
  void SimWireMicroBooNE::MakeADCVec(std::vector<short>& adcvec, std::vector<float> const& noisevec, 
                                     std::vector<double> const& chargevec, float ped_mean) const {


    for(unsigned int i = 0; i < fNTimeSamples; ++i) {

       float adcval = noisevec[i] + chargevec[i] + ped_mean;

      //allow for ADC saturation
      if ( adcval > adcsaturation )
	adcval = adcsaturation;
      //don't allow for "negative" saturation
      if ( adcval < 0 )
	adcval = 0;

      adcvec[i] = (unsigned short)TMath::Nint(adcval);

    }// end loop over signal size

    // compress the adc vector using the desired compression scheme,
    // if raw::kNone is selected nothing happens to adcvec
    // This shrinks adcvec, if fCompression is not kNone.
    raw::Compress(adcvec, fCompression);
  }


  //-------------------------------------------------
  void SimWireMicroBooNE::GenNoiseInTime(std::vector<float> &noise, double noise_factor) const
  {
    //ART random number service
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine("noise");
    CLHEP::RandGaussQ rGauss(engine, 0.0, noise_factor);

    //In this case noise_factor is a value in ADC counts
    //It is going to be the Noise RMS
    //loop over all bins in "noise" vector
    //and insert random noise value
    for (unsigned int i=0; i<noise.size(); i++)
      noise.at(i) = rGauss.fire();
  }


  //-------------------------------------------------
  void SimWireMicroBooNE::GenNoiseInFreq(std::vector<float> &noise, double noise_factor) const
  {
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine("noise");
    CLHEP::RandFlat flat(engine,-1,1);

    if(noise.size() != fNTicks)
      throw cet::exception("SimWireMicroBooNE")
      << "\033[93m"
      << "Frequency noise vector length must match fNTicks (FFT size)"
      << " ... " << noise.size() << " != " << fNTicks
      << "\033[00m"
      << std::endl;

    // noise in frequency space
    std::vector<TComplex> noiseFrequency(fNTicks/2+1, 0.);

    double pval = 0.;
    double lofilter = 0.;
    double phase = 0.;
    double rnd[2] = {0.};

    // width of frequencyBin in kHz
    double binWidth = 1.0/(fNTicks*fSampleRate*1.0e-6);
    for(size_t i=0; i< fNTicks/2+1; ++i){
      // exponential noise spectrum
      flat.fireArray(2,rnd,0,1);
      //if not from histo or in time --> then hardcoded freq. spectrum
      if( !fGetNoiseFromHisto )
      {
        pval = noise_factor*exp(-(double)i*binWidth/fNoiseWidth);
        // low frequency cutoff
        lofilter = 1.0/(1.0+exp(-(i-fLowCutoff/binWidth)/0.5));
        // randomize 10%
        
        pval *= lofilter*((1-fNoiseRand)+2*fNoiseRand*rnd[0]);
      }
      
      
      else
      {
        // histogram starts in bin 1!
        pval = fNoiseHist->GetBinContent(i+1)*((1-fNoiseRand)+2*fNoiseRand*rnd[0])*noise_factor;
        //mf::LogInfo("SimWireMicroBooNE")  << " pval: " << pval;
      }
      phase = rnd[1]*2.*TMath::Pi();
      TComplex tc(pval*cos(phase),pval*sin(phase));
      noiseFrequency.at(i) += tc;
    }
    
    
    // mf::LogInfo("SimWireMicroBooNE") << "filled noise freq";
    
    // inverse FFT MCSignal
    art::ServiceHandle<util::LArFFT> fFFT;
    fFFT->DoInvFFT(noiseFrequency, noise);
    
    noiseFrequency.clear();
    
    // multiply each noise value by fNTicks as the InvFFT
    // divides each bin by fNTicks assuming that a forward FFT
    // has already been done.
    // Also need to scale so that noise RMS matches that asked
    // in fhicl parameter (somewhat arbitrary scaling otherwise)
    // harcode this scaling factor (~20) for now
    for(unsigned int i = 0; i < noise.size(); ++i) noise.at(i) *= 1.*(fNTicks/20.);
    
  }

  //---------------------------------------------------------
  Double_t rpfn_f1(double waveform_size, Double_t fitpar[], Double_t x) {
    return fitpar[0] + ( fitpar[1] + fitpar[2] * x * waveform_size/2.) * TMath::Exp(-fitpar[3] * TMath::Power( x * waveform_size/2., fitpar[4]));
  }
  
  void SimWireMicroBooNE::GenNoisePostFilter(std::vector<float> &noise, double noise_factor, size_t view, int chan)
  {
    // noise is a vector of size fNTicks, which is the number of ticks
    const size_t waveform_size = noise.size();
    

    if(_pfn_shaping_time_v.size()<=view || _pfn_shaping_time_v[view]<0)
      throw cet::exception("SimWireMicroBooNE")
	<< "GenNoisePostFilter encounters unknown view!" << std::endl;

    Double_t ShapingTime = _pfn_shaping_time_v[view];

    if(_pfn_rho_v.empty()) _pfn_rho_v.resize(waveform_size);
    if(_pfn_value_re.empty()) _pfn_value_re.resize(waveform_size);
    if(_pfn_value_im.empty()) _pfn_value_im.resize(waveform_size);

    //**Setting lambda**//
    Double_t params[1] = {0.};
    Double_t fitpar[5] = {0.};
    
    if(ShapingTime==2.0) {
      params[0] = 3.3708; //2us
      fitpar[0] = 4.27132e+01;
      fitpar[1] = 6.22750e+02;
      fitpar[2] = -2.53535e-01;
      fitpar[3] = 8.07578e-05;
      fitpar[4] = 1.35510e+00;
    }
    else if(ShapingTime==1.0) {
      params[0] = 3.5125; //1us 
      fitpar[0] = 14.4;
      fitpar[1] = 35.1;
      fitpar[2] = 0.049;
      fitpar[3] = 6.0e-9;
      fitpar[4] = 2.4;
    }else
      throw cet::exception("SimWireMicroBooNE") << "<<" << __FUNCTION__ << ">> not supported shaping time " << ShapingTime << std::endl;

    Int_t n = waveform_size;

    TVirtualFFT::SetTransform(0);
   
    // seed gamma-distibuted random number with mean params[0]
    // replacing continuous Poisson distribution from ROOT
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::gamma_distribution<double> distribution(params[0]);

    // For every tick...
    for(size_t i=0; i<waveform_size; i++){

      Double_t freq;
      if (i < n/2.){
        freq = (i)*2./n; //2 MHz digitization
      }else{
        freq = (n-i)*2./n;
      }

      // Draw gamma-dstributed random number
      gammaRand = distribution(generator); 

      // Define FFT parameters
      _pfn_rho_v[i] = rpfn_f1(waveform_size, fitpar, freq) * gammaRand /params[0];
      Double_t rho = _pfn_rho_v[i];
      Double_t phi = gRandom->Uniform(0,1) * 2. * TMath::Pi();

      _pfn_value_re[i] = rho*cos(phi)/((double)waveform_size);
      _pfn_value_im[i] = rho*sin(phi)/((double)waveform_size);
    }

    // Inverse FFT
    if(!_pfn_ifft) _pfn_ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
    _pfn_ifft->SetPointsComplex(&_pfn_value_re[0],&_pfn_value_im[0]);
    _pfn_ifft->Transform();

    // Produce fit histogram from the FFT for the real part
    TH1 *fb = 0;
    fb = TH1::TransformHisto(_pfn_ifft,fb,"Re");

    // Get wire length 
    geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();
    std::vector<geo::WireID> wireIDs = geom->ChannelToWire(chan);
    geo::WireGeo const& wire = geom->Wire(wireIDs.front());
    double wirelength = wire.HalfL() * 2;
    
    // Calculate RMS -----------------------------------------------------
    // Calculating using the 16th, 50th, and 84th percentiles.
    // Because the signal is expected to be above the 84th percentile, this 
    // effectively vetos the signal.
  
    Double_t min = fb->GetMinimum();
    Double_t max = fb->GetMaximum();	 
    TH1F* h_rms = new TH1F("h_rms", "h_rms", Int_t(10*(max-min+1)), min, max+1);

    for(size_t i=0; i < waveform_size; ++i){
      h_rms->Fill(fb->GetBinContent(i+1));
    }
    
    double par[3];
    double rms_quantilemethod = 0.0;
    if (h_rms->GetSum()>0){
    double xq = 0.5-0.34;
    h_rms->GetQuantiles(1, &par[0], &xq);
    
    xq = 0.5;
    h_rms->GetQuantiles(1, &par[1], &xq);

    xq = 0.5+0.34;
    h_rms->GetQuantiles(1, &par[2], &xq);
    
    rms_quantilemethod = sqrt((pow(par[1]-par[0],2)+pow(par[2]-par[1],2))/2.);
    }
    
    // Used for calculation of baseline
    // rms.push_back(std::move(rms_quantilemethod));

    // Scaling noise RMS with wire length dependance
    double baseline = 1.79349;
    double scalefactor = (rms_quantilemethod/baseline)*(1.019+0.00173*(wirelength));
    for(size_t i=0; i<waveform_size; ++i) {
      noise[i] = fb->GetBinContent(i+1)*scalefactor;
    }
    
    delete fb;
    delete h_rms;
  }
  
}
