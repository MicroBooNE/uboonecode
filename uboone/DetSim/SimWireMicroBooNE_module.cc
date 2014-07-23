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

    std::string            fDriftEModuleLabel;///< module making the ionization electrons
    raw::Compress_t        fCompression;      ///< compression type to use

    double                 fNoiseFact;        ///< Flexible noise factor, taken from induction and collection as wires change
    double                 fNoiseFactInd;     ///< noise scale factor for Induction channels (ADC RMS for gaussian noise)
    double                 fNoiseFactColl;    ///< noise scale factor for Collection channels (ADC RMS for gaussian noise)
    double                 fShapingTime;      ///< Shaping time in usec
    double                 fNoiseWidth;       ///< exponential noise width (kHz) 
    double                 fNoiseRand;        ///< fraction of random "wiggle" in noise in freq. spectrum
    double                 fLowCutoff;        ///< low frequency filter cutoff (kHz)
    size_t                 fNTicks;           ///< number of ticks of the clock
    double                 fSampleRate;       ///< sampling rate in ns
    unsigned int           fNTimeSamples;     ///< number of ADC readout samples in all readout frames (per event)
    float                  fCollectionPed;    ///< ADC value of baseline for collection plane
    float                  fInductionPed;     ///< ADC value of baseline for induction plane
    float                  fBaselineRMS;      ///< ADC value of baseline RMS within each channel                    
    TH1D*                  fNoiseDist;        ///< distribution of noise counts
    bool fGetNoiseFromHisto;                  ///< if True -> Noise from Histogram of Freq. spectrum
    bool fGenNoiseInTime;                     ///< if True -> Noise with Gaussian dsitribution in Time-domain 
    bool fGenNoise;                           ///< if True -> Gen Noise. if False -> Skip noise generation entierly
    std::string fNoiseFileFname; 
    std::string fNoiseHistoName; 
    TH1D*             fNoiseHist;             ///< distribution of noise counts
    float                  fASICGain;

    
    std::map< double, int > fShapingTimeOrder;
    std::string fTrigModName;                 ///< Trigger data product producer name
    //define max ADC value - if one wishes this can
    //be made a fcl parameter but not likely to ever change
    const float adcsaturation = 4095;
    
    ::util::ElecClock fClock; ///< TPC electronics clock
    
  }; // class SimWireMicroBooNE
  
  DEFINE_ART_MODULE(SimWireMicroBooNE)

  //-------------------------------------------------
  SimWireMicroBooNE::SimWireMicroBooNE(fhicl::ParameterSet const& pset)
  : fNoiseHist(0)
  {
    this->reconfigure(pset);

    produces< std::vector<raw::RawDigit>   >();

    fCompression = raw::kNone;
    TString compression(pset.get< std::string >("CompressionType"));
    if(compression.Contains("Huffman",TString::kIgnoreCase)) fCompression = raw::kHuffman;    

    // get the random number seed, use a random default if not specified    
    // in the configuration file.  
    unsigned int seed = pset.get< unsigned int >("Seed", sim::GetRandomNumberSeed());

    createEngine(seed);
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
    fCollectionPed    = p.get< float               >("CollectionPed");
    fInductionPed     = p.get< float               >("InductionPed");
    fBaselineRMS      = p.get< float               >("BaselineRMS");
    fTrigModName      = p.get< std::string         >("TrigModName");

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

    return;
  }

  //-------------------------------------------------
  void SimWireMicroBooNE::beginJob() 
  { 

    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;

    fNoiseDist  = tfs->make<TH1D>("Noise", ";Noise  (ADC);", 1000,   -30., 30.);

    art::ServiceHandle<util::LArFFT> fFFT;
    fNTicks = fFFT->FFTSize();

   if ( fNTicks%2 != 0 ) 
      LOG_DEBUG("SimWireMicroBooNE") << "Warning: FFTSize not a power of 2. "
                                     << "May cause issues in (de)convolution.\n";

    if ( fNTimeSamples > fNTicks ) 
      mf::LogError("SimWireMircoBooNE") << "Cannot have number of readout samples "
                                        << "greater than FFTSize!";
    
    return;

  }

  //-------------------------------------------------
  void SimWireMicroBooNE::endJob() 
  {}

  void SimWireMicroBooNE::produce(art::Event& evt)
  {

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

    std::vector<const sim::SimChannel*> chanHandle;
    evt.getView(fDriftEModuleLabel,chanHandle);

    //Get fIndShape and fColShape from SignalShapingService, on the fly
    art::ServiceHandle<util::SignalShapingServiceMicroBooNE> sss;
    //get ASIC Gain and Noise in ADCatLowestGain:
    fASICGain      = sss->GetASICGain();
    fShapingTime   = sss->GetShapingTime();
    //Check that shaping time is an allowed value
    //If so, Pick out noise factor 
    //If not, through exception
    if ( fShapingTimeOrder.find( fShapingTime ) != fShapingTimeOrder.end() ){
      fNoiseFactInd   = sss->GetNoiseFactInd().at( fShapingTimeOrder.find( fShapingTime )->second );
      fNoiseFactColl  = sss->GetNoiseFactColl().at( fShapingTimeOrder.find( fShapingTime )->second );
    }
    else{//Throw exception...
      throw cet::exception("SimWireMicroBooNE")
        << "\033[93m"
        << "Shaping Time received from signalservices_microboone.fcl is not one of allowed values"
	<< std::endl
        << "Allowed values: 0.5, 1.0, 2.0, 3.0 usec"
        << "\033[00m"
        << std::endl;
    }
    //Take into account ASIC Gain
    fNoiseFactInd *= fASICGain/4.7;
    fNoiseFactColl *= fASICGain/4.7;

    // make a vector of const sim::SimChannel* that has same number
    // of entries as the number of channels in the detector
    // and set the entries for the channels that have signal on them
    // using the chanHandle
    std::vector<const sim::SimChannel*> channels(geo->Nchannels(),nullptr);
    for(size_t c = 0; c < chanHandle.size(); ++c){
      channels.at(chanHandle.at(c)->Channel()) = chanHandle.at(c);
    }
    
    const auto NChannels = geo->Nchannels();

    // vectors for working
    std::vector<short>    adcvec(fNTimeSamples, 0);
    std::vector<double>   chargeWork(fNTicks,0.);

    // make a unique_ptr of sim::SimDigits that allows ownership of the produced
    // digits to be transferred to the art::Event after the put statement below
    std::unique_ptr< std::vector<raw::RawDigit>> digcol(new std::vector<raw::RawDigit>);
    digcol->reserve(NChannels);

    unsigned int chan = 0; 
    art::ServiceHandle<util::LArFFT> fFFT;

    //LOOP OVER ALL CHANNELS
    std::map<int,double>::iterator mapIter;      
    for(chan = 0; chan < geo->Nchannels(); chan++) {

      std::fill(chargeWork.begin(), chargeWork.end(), 0.);

      // get the sim::SimChannel for this channel
      const sim::SimChannel* sc = channels.at(chan);

      if( sc ){

	int time_offset = sss->FieldResponseTOffset(chan);
        // loop over the tdcs and grab the number of electrons for each
        for(int t = 0; t < (int)(chargeWork.size()); ++t) {

          int tdc = ts->TPCTick2TDC(t) + time_offset;

          // continue if tdc < 0
          if( tdc < 0 ) continue;

          chargeWork.at(t) = sc->Charge(tdc);

        }

        // Convolve charge with appropriate response function 
        sss->Convolute(chan,chargeWork);

      }

      //Pedestal determination
      float ped_mean = fCollectionPed;
      geo::SigType_t sigtype = geo->SignalType(chan);
      if (sigtype == geo::kInduction){
        ped_mean = fInductionPed;
	fNoiseFact = fNoiseFactInd;
      }
      else if (sigtype == geo::kCollection){
        ped_mean = fCollectionPed;
	fNoiseFact = fNoiseFactColl;
      }
      
      //Generate Noise:
      std::vector<float> noisetmp(fNTicks,0.);
      if (fGenNoise){
        if (fGenNoiseInTime)
          GenNoiseInTime(noisetmp);
        else
          GenNoiseInFreq(noisetmp);
      }


      //slight variation on ped on order of RMS of baselien variation
      art::ServiceHandle<art::RandomNumberGenerator> rng;
      CLHEP::HepRandomEngine &engine = rng->getEngine();
      CLHEP::RandGaussQ rGaussPed(engine, 0.0, fBaselineRMS);
      ped_mean += rGaussPed.fire();
      
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
        if (i%1000==0)
          fNoiseDist->Fill(noisetmp[i]);

        //allow for ADC saturation
        if ( adcval > adcsaturation )
          adcval = adcsaturation;
        //don't allow for "negative" saturation
        if ( adcval < 0 )
          adcval = 0;

        adcvec[i] = (unsigned short)(adcval);

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
      digcol->push_back(rd);

    }// end loop over channels      

    evt.put(std::move(digcol));
    return;
  }
  
  //-------------------------------------------------                                                                                                 
  void SimWireMicroBooNE::GenNoiseInTime(std::vector<float> &noise)
  {
    //ART random number service                                                                                                                       
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
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
    CLHEP::HepRandomEngine &engine = rng->getEngine();
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
          
          pval = fNoiseHist->GetBinContent(i)*((1-fNoiseRand)+2*fNoiseRand*rnd[0])*fNoiseFact; 
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
  
