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
#include "artextensions/SeedService/SeedService.hh"

// LArSoft libraries
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "RawData/TriggerData.h"
#include "Simulation/SimChannel.h"
#include "Geometry/Geometry.h"
#include "Utilities/LArFFT.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/TimeService.h"
#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"
#include "Simulation/sim.h"
#include "CalibrationDBI/WebDBI/DetPedestalRetrievalAlg.h"
#include "uboone/Database/UBooneIOVTimeStamp.h"



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

    void GenNoiseInTime(std::vector<float> &noise);
    void GenNoiseInFreq(std::vector<float> &noise);

    size_t fNChannels;
    size_t fNViews;
    std::vector<std::vector<size_t> > fNResponses;       ///< this is the number of *active* responses, <= number of responses
    std::string            fDriftEModuleLabel;///< module making the ionization electrons
    raw::Compress_t        fCompression;      ///< compression type to use

    double                 fNoiseFact;        ///< Flexible noise factor, taken from induction and collection as wires change
    DoubleVec              fNoiseFactVec;
    //double                 fNoiseFactInd;     ///< noise scale factor for Induction channels (ADC RMS for gaussian noise)
    //double                 fNoiseFactColl;    ///< noise scale factor for Collection channels (ADC RMS for gaussian noise)
    double                 fShapingTime;      ///< Shaping time in usec
    double                 fNoiseWidth;       ///< exponential noise width (kHz)
    double                 fNoiseRand;        ///< fraction of random "wiggle" in noise in freq. spectrum
    double                 fLowCutoff;        ///< low frequency filter cutoff (kHz)
    size_t                 fNTicks;           ///< number of ticks of the clock
    double                 fSampleRate;       ///< sampling rate in ns
    
    unsigned int           fNTimeSamples;     ///< number of ADC readout samples in all readout frames (per event)                 

    TH1D*                  fNoiseDistColl;    ///< distribution of noise counts
    TH1D*                  fNoiseDistInd;     ///< distribution of noise counts
    bool fGetNoiseFromHisto;                  ///< if True -> Noise from Histogram of Freq. spectrum
    bool fGenNoiseInTime;                     ///< if True -> Noise with Gaussian dsitribution in Time-domain
    bool fGenNoise;                           ///< if True -> Gen Noise. if False -> Skip noise generation entierly

    std::string fNoiseFileFname;
    std::string fNoiseHistoName;
    TH1D*             fNoiseHist;             ///< distribution of noise counts
    float                  fASICGain;

    std::map< double, int > fShapingTimeOrder;
    std::string fTrigModName;                 ///< Trigger data product producer name
    
    lariov::DetPedestalRetrievalAlg fPedestalRetrievalAlg;      ///< Pedestal Retrieval algorithm
    

    bool        fTest; // for forcing a test case
    std::vector<sim::SimChannel> fTestSimChannel_v;
    size_t      fTestWire;
    std::vector<size_t> fTestIndex;
    std::vector<double> fTestCharge;

    int         fSample; // for histograms, -1 means no histos

    //define max ADC value - if one wishes this can
    //be made a fcl parameter but not likely to ever change
    const float adcsaturation = 4095;

    ::util::ElecClock fClock; ///< TPC electronics clock

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

  }; // class SimWireMicroBooNE

  namespace {
    size_t _ch = 0;
    size_t _wr = 0;
  }

  DEFINE_ART_MODULE(SimWireMicroBooNE)

  //-------------------------------------------------
  SimWireMicroBooNE::SimWireMicroBooNE(fhicl::ParameterSet const& pset)
  : fNoiseHist(0)
    , fPedestalRetrievalAlg(pset.get<fhicl::ParameterSet>("DetPedestalRetrievalAlg"))
  {
    this->reconfigure(pset);

    produces< std::vector<raw::RawDigit>   >();

    fCompression = raw::kNone;
    TString compression(pset.get< std::string >("CompressionType"));
    if(compression.Contains("Huffman",TString::kIgnoreCase)) fCompression = raw::kHuffman;

    // create a default random engine; obtain the random seed from SeedService,
    // unless overridden in configuration with key "Seed" and "SeedPedestal"
    art::ServiceHandle<artext::SeedService> Seeds;
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
    fGenNoiseInTime   = p.get< bool                >("GenNoiseInTime");
    fGenNoise         = p.get< bool                >("GenNoise");

    fTrigModName      = p.get< std::string         >("TrigModName");
    fTest             = p.get<bool                 >("Test");
    fTestWire         = p.get< size_t              >("TestWire");
    fTestIndex        = p.get< std::vector<size_t> >("TestIndex");
    fTestCharge       = p.get< std::vector<double> >("TestCharge");
    if(fTestIndex.size() != fTestCharge.size())
      throw cet::exception(__FUNCTION__)<<"# test pulse mismatched: check TestIndex and TestCharge fcl parameters...";
    fSample           = p.get<int                  >("Sample");
    fPedestalRetrievalAlg.Reconfigure(p.get<fhicl::ParameterSet>("DetPedestalRetrievalAlg"));


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
    art::ServiceHandle<util::DetectorProperties> detprop;
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

    fNoiseDistColl  = tfs->make<TH1D>("NoiseCollection", ";Noise on Collection Wires (ADC);", 1000,   -30., 30.);
    fNoiseDistInd  = tfs->make<TH1D>("NoiseInduction", ";Noise on Induction Wires (ADC);", 1000,   -30., 30.);

    art::ServiceHandle<util::LArFFT> fFFT;
    fNTicks = fFFT->FFTSize();

    if ( fNTicks%2 != 0 )
      LOG_DEBUG("SimWireMicroBooNE") << "Warning: FFTSize " << fNTicks << " not a power of 2. "
      << "May cause issues in (de)convolution.\n";

    if ( fNTimeSamples > fNTicks )
      mf::LogError("SimWireMircoBooNE") << "Cannot have number of readout samples "
      << fNTimeSamples << " greater than FFTSize " << fNTicks << "!";

       
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
  {}

  void SimWireMicroBooNE::produce(art::Event& evt)
  {

    //update database cache
    fPedestalRetrievalAlg.Update( evt.time().value() );

    art::ServiceHandle<util::LArFFT> fFFT;
    fFFT->ReinitializeFFT(fNTimeSamples,fFFT->FFTOptions(),fFFT->FFTFitBins());
    fNTicks = fFFT->FFTSize();

    if ( fNTicks%2 != 0 )
      LOG_DEBUG("SimWireMicroBooNE") << "Warning: FFTSize " << fNTicks << " not a power of 2. "
      << "May cause issues in (de)convolution.\n";

    if ( fNTimeSamples > fNTicks )
      mf::LogError("SimWireMicroBooNE") << "Cannot have number of readout samples "
      << fNTimeSamples << " greater than FFTSize " << fNTicks << "!";


    art::ServiceHandle<art::TFileService> tfs;

    art::ServiceHandle<util::TimeService> ts;
    // In case trigger simulation is run in the same job...
    ts->preProcessEvent(evt);

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
    //unsigned int signalSize = fNTicks;

    //Get fIndShape and fColShape from SignalShapingService, on the fly
    art::ServiceHandle<util::SignalShapingServiceMicroBooNE> sss;
    //get ASIC Gain and Noise in ADCatLowestGain:

    //fASICGain      = sss->GetASICGain();
    //fShapingTime   = sss->GetShapingTime();

    fNResponses = sss->GetNActiveResponses();
    fNViews     = fNResponses[0].size();



    fNChannels = geo->Nchannels();


    // make a vector of const sim::SimChannel* that has same number
    // of entries as the number of channels in the detector
    // and set the entries for the channels that have signal on them
    // using the chanHandle
    std::vector<const sim::SimChannel*> channels(fNChannels,nullptr);

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

    //const auto NChannels = geo->Nchannels();

    // vectors for working
    std::vector<short>    adcvec(fNTimeSamples, 0);
    std::vector<double>   chargeWork(fNTicks,0.);
    std::vector<double>   tempWork(fNTicks,0.);

    // make a unique_ptr of sim::SimDigits that allows ownership of the produced
    // digits to be transferred to the art::Event after the put statement below
    std::unique_ptr< std::vector<raw::RawDigit>> digcol(new std::vector<raw::RawDigit>);
    digcol->reserve(fNChannels);

    _wr = 0;
    unsigned int chan = 0;
    //art::ServiceHandle<util::LArFFT> fFFT;

    //std::vector<std::vector<std::vector<ResponseParams* > > > responseParamsVec(fNChannels);
    std::vector<std::vector<std::vector<std::unique_ptr<ResponseParams> > > > responseParamsVec(fNChannels);


    _ch = 0;
    for (auto& channel : responseParamsVec) {
      size_t view = (size_t)geo->View(_ch);
      channel.resize(2*fNResponses[0][view]-1);
    }

    size_t view = 0;


    //LOOP OVER ALL CHANNELS
    // In this version we assume that adjacent channels <-> adjacent wires, in the same plane/view
    // Is this always true?
    //std::map<int,double>::iterator mapIter;


    for(chan = 0; chan < fNChannels; chan++) {
      auto wid = geo->ChannelToWire(chan);
      view = (size_t)geo->View(chan);

      // get the sim::SimChannel for this channel
      const sim::SimChannel* sc = channels.at(chan);
      if( !sc ) continue;

      // remove the time offset
      int time_offset = 0;//sss->FieldResponseTOffset(chan);

      // loop over the tdcs and grab the number of electrons for each

      for(int t = 0; t < (int)(chargeWork.size()); ++t) {

        int tdc = ts->TPCTick2TDC(t);
        // continue if tdc < 0
        if( tdc < 0 ) continue;
        double charge = sc->Charge(tdc);
        if(charge==0) continue;

        // Apply artificial time offset to take care of field response convolution
        // wrap the negative times to the end of the buffer
        // The offset should be taken care of in the shaping service, by shifting the response.

        int raw_digit_index =
        ( (t + time_offset) >= 0 ? t+time_offset : (chargeWork.size() + (t+time_offset)) );

        if(raw_digit_index <= 0 || raw_digit_index >= (int)(chargeWork.size())) continue;

        // here fill ResponseParams... all the wires!
        for(int wire = -(fNResponses[0][view]-1); wire<(int)fNResponses[0][view]; ++wire) {
          auto wireIndex = (size_t)wire+fNResponses[0][view] - 1;
          int wireChan = (int)chan + wire;
          if(wireChan<0 || wireChan>= (int)fNChannels) continue;
          if ((size_t)geo->View(wireChan)!=view) continue;

          responseParamsVec[wireChan][wireIndex].emplace_back(new ResponseParams(charge, raw_digit_index));
        } // loop over wires
      } // loop over tdcs
    } // loop over channels

    double slope0[5] = { 0., 2.1575, 6.4725 , 13.946, 40.857};
    double t0[5] =     { 4450., 6107., 6170., 6305., 6695. };
    double wire0[3] =  { 337., 332., -0.7 };
    double factor[3] = { 2.0, 2.0, 1.0 };

    int tickCut = 250;

    // loop over the collected responses
    //   this is needed because hits generate responses on adjacent wires!

    for(chan = 0; chan < fNChannels; chan++) {
      // the effect of the commented line below was to not make waveforms or deposit induced charge on
      //   any wire that didn't originally have deposited charge.
      // So induced charge didn't work in test mode, and the ghost hits before the start of the track
      //   were missing.
      // Thanks, Brandon!
      //if(!channels.at(chan)) continue;

      auto wid = geo->ChannelToWire(chan);
      size_t wireNum = wid[0].Wire;

      std::fill(chargeWork.begin(), chargeWork.end(), 0.);



      //const sim::SimChannel* sc = channels.at(chan);

      fASICGain      = sss->GetASICGain(chan);     //Jyoti - to read different gain for U,V & Y planes
      fShapingTime   = sss->GetShapingTime(chan);  //Jyoti - to read different shaping time for U,V & Y planes

      fNoiseFactVec.resize(fNViews);
      auto tempNoiseVec = sss->GetNoiseFactVec();
      if ( fShapingTimeOrder.find( fShapingTime ) != fShapingTimeOrder.end() ){

        size_t i = 0;
        for (auto& item : tempNoiseVec) {
          fNoiseFactVec[i]   = item.at( fShapingTimeOrder.find( fShapingTime )->second );
          fNoiseFactVec[i] *= fASICGain/4.7;
          ++i;
        }
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
      //to be moved

      //Take into account ASIC Gain
      //    fNoiseFactInd *= fASICGain/4.7;
      //        fNoiseFactColl *= fASICGain/4.7;
      // get the sim::SimChannel for this channel

      auto& thisChan = responseParamsVec[chan];
      view = (size_t)geo->View(chan);

      int tick0 = 0;
      if(fSample>=0) tick0 = t0[fSample] - factor[view]*slope0[fSample]*(wireNum-wire0[view]) + 0.5;

      for(int wire=-(fNResponses[0][view]-1); wire<(int)fNResponses[0][view];++wire) {
        int wireChan = chan + wire;
        if(wireChan<0 || wireChan>= (int)fNChannels) continue;
        if ((size_t)geo->View(wireChan)!=view) continue;
        size_t wireIndex = (size_t)(wire + (int)fNResponses[0][view] - 1);
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
        // i	size_t		convolve it with the appropriate response function

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


      }



      //Generate Noise:
      //Pedestal determination and random gaussian variation
      lariov::DetPedestal pedestal = fPedestalRetrievalAlg.Pedestal(chan);
      art::ServiceHandle<art::RandomNumberGenerator> rng;
      CLHEP::HepRandomEngine &engine = rng->getEngine("pedestal");
      CLHEP::RandGaussQ rGaussPed(engine, 0.0, pedestal.PedRms());
      float ped_mean = pedestal.PedMean() + rGaussPed.fire();

      //Generate Noise
      //geo::SigType_t sigtype = geo->SignalType(chan);
      fNoiseFact = fNoiseFactVec[view];

      std::vector<float> noisetmp(fNTicks,0.);
      if (fGenNoise){
        if (fGenNoiseInTime)
          GenNoiseInTime(noisetmp);
        else
          GenNoiseInFreq(noisetmp);
      }

      // resize the adcvec to be the correct number of time samples
      // (the compression from the previous loop might have changed size)
      // almost no-op after the first loop, since the memory is already there
      adcvec.resize(fNTimeSamples);

      // rather than checking each single one of the sequential accesses we are
      // going to perform, we do it once for all here:
      if (noisetmp.size() < fNTimeSamples)
        throw std::range_error("SimWireMicroBooNE: noisetmp vector too small");
      if (chargeWork.size() < fNTimeSamples)
        throw std::range_error("SimWireMicroBooNE: chargeWork vector too small");
      if (adcvec.size() < fNTimeSamples)
        throw std::range_error("SimWireMicroBooNE: adcvec vector too small");

      for(unsigned int i = 0; i < fNTimeSamples; ++i){
        float adcval = noisetmp[i] + chargeWork[i] + ped_mean;

        //Add Noise to NoiseDist Histogram
        geo::SigType_t sigtype = geo->SignalType(chan);
        if (i%1000==0){
          if (sigtype == geo::kCollection)
            fNoiseDistColl->Fill(noisetmp[i]);
          if (sigtype == geo::kInduction)
            fNoiseDistInd->Fill(noisetmp[i]);
        }

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

      // add this digit to the collection;
      // adcvec is copied, not moved: in case of compression, adcvec will show
      // less data: e.g. if the uncompressed adcvec has 9600 items, after
      // compression it will have maybe 5000, but the memory of the other 4600
      // is still there, although unused; a copy of adcvec will instead have
      // only 5000 items. All 9600 items of adcvec will be recovered for free
      // and used on the next loop.
      raw::RawDigit rd(chan, fNTimeSamples, adcvec, fCompression);
      rd.SetPedestal(ped_mean);
      digcol->push_back(std::move(rd)); // we do move the raw digit copy, though
      // std::cout<< "Xin2 " << chan << " " << fNChannels << std::endl;
    }// end loop over channels



    evt.put(std::move(digcol));
    return;
  }

  //-------------------------------------------------
  void SimWireMicroBooNE::GenNoiseInTime(std::vector<float> &noise)
  {
    //ART random number service
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine("noise");
    CLHEP::RandGaussQ rGauss(engine, 0.0, fNoiseFact);

    //In this case fNoiseFact is a value in ADC counts
    //It is going to be the Noise RMS
    //loop over all bins in "noise" vector
    //and insert random noise value
    for (unsigned int i=0; i<noise.size(); i++)
      noise.at(i) = rGauss.fire();
  }


  //-------------------------------------------------
  void SimWireMicroBooNE::GenNoiseInFreq(std::vector<float> &noise)
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
        pval = fNoiseFact*exp(-(double)i*binWidth/fNoiseWidth);
        // low frequency cutoff
        lofilter = 1.0/(1.0+exp(-(i-fLowCutoff/binWidth)/0.5));
        // randomize 10%
        
        pval *= lofilter*((1-fNoiseRand)+2*fNoiseRand*rnd[0]);
      }
      
      
      else
      {
        // histogram starts in bin 1!
        pval = fNoiseHist->GetBinContent(i+1)*((1-fNoiseRand)+2*fNoiseRand*rnd[0])*fNoiseFact;
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
    //Also need to scale so that noise RMS matches that asked
    //in fhicl parameter (somewhat arbitrary scaling otherwise)
    //harcode this scaling factor (~20) for now
    for(unsigned int i = 0; i < noise.size(); ++i) noise.at(i) *= 1.*(fNTicks/20.);
    
  }
  
  
}

