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
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <bitset>

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Utilities/LArFFT.h"
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "RawData/TriggerData.h"
#include "Utilities/LArProperties.h"
#include "Utilities/TimeService.h"
#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"
#include "Geometry/Geometry.h"
#include "Simulation/sim.h"
#include "Simulation/SimChannel.h"
#include "Utilities/DetectorProperties.h"

#include "TMath.h"
#include "TComplex.h"
#include "TString.h"
#include "TH2.h"
#include "TH1D.h"
#include "TFile.h"
#include "TRandom.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

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

    std::vector<float>     GenNoiseInTime();
    std::vector<float>     GenNoiseInFreq();

    std::string            fDriftEModuleLabel;///< module making the ionization electrons
    raw::Compress_t        fCompression;      ///< compression type to use

    double                 fNoiseFact;        ///< noise scale factor (ADC RMS for gaussian noise)
    double                 fNoiseWidth;       ///< exponential noise width (kHz) 
    double                 fNoiseRand;        ///< fraction of random "wiggle" in noise in freq. spectrum
    double                 fLowCutoff;        ///< low frequency filter cutoff (kHz)
    int                    fNTicks;           ///< number of ticks of the clock
    double                 fSampleRate;       ///< sampling rate in ns
    unsigned int           fNTimeSamples;     ///< number of ADC readout samples in all readout frames (per event)
    float                  fCollectionPed;    ///< ADC value of baseline for collection plane
    float                  fInductionPed;     ///< ADC value of baseline for induction plane
    float                  fBaselineRMS;      ///< ADC value of baseline RMS within each channel                    
    std::vector<double>    fChargeWork;
    TH1D*                  fNoiseDist;        ///< distribution of noise counts
    bool fGetNoiseFromHisto;                  ///< if True -> Noise from Histogram of Freq. spectrum
    bool fGenNoiseInTime;                     ///< if True -> Noise with Gaussian dsitribution in Time-domain 
    bool fGenNoise;                           ///< if True -> Gen Noise. if False -> Skip noise generation entierly
    std::string fNoiseFileFname; 
    std::string fNoiseHistoName; 
    TH1D*             fNoiseHist;             ///< distribution of noise counts

    std::string fTrigModName;                 ///< Trigger data product producer name
    //define max ADC value - if one wishes this can
    //be made a fcl parameter but not likely to ever change
    const float adcsaturation = 4095;
    
    ::util::ElecClock fClock; ///< TPC electronics clock
    
  }; // class SimWireMicroBooNE
  
  DEFINE_ART_MODULE(SimWireMicroBooNE)

  //-------------------------------------------------
  SimWireMicroBooNE::SimWireMicroBooNE(fhicl::ParameterSet const& pset)
  {
    std::cout<<"\033[95m"<<"Constructor Start!!!"<<"\033[00m"<<std::endl;
    this->reconfigure(pset);

    produces< std::vector<raw::RawDigit>   >();

    fCompression = raw::kNone;
    TString compression(pset.get< std::string >("CompressionType"));
    if(compression.Contains("Huffman",TString::kIgnoreCase)) fCompression = raw::kHuffman;    

    // get the random number seed, use a random default if not specified    
    // in the configuration file.  
    unsigned int seed = pset.get< unsigned int >("Seed", sim::GetRandomNumberSeed());

    createEngine(seed);
    std::cout<<"\033[95m"<<"Constructor End!!!"<<"\033[00m"<<std::endl;
  }

  //-------------------------------------------------
  SimWireMicroBooNE::~SimWireMicroBooNE()
  {
    fChargeWork.clear();
  }

  //-------------------------------------------------
  void SimWireMicroBooNE::reconfigure(fhicl::ParameterSet const& p) 
  {
    std::cout<<"\033[95m"<<"reconfigure Start!!!"<<"\033[00m"<<std::endl;
    fDriftEModuleLabel= p.get< std::string         >("DriftEModuleLabel");
    fNoiseFact        = p.get< double              >("NoiseFact");
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
	throw cet::exception("SimWireMicroBooNE") << " Could not find noise histogram in Root file\n";
      in->Close();
    
      }
    //detector properties information
    art::ServiceHandle<util::DetectorProperties> detprop;
    fNTimeSamples  = detprop->NumberTimeSamples();

    std::cout<<"\033[95m"<<"reconfigure End!!!"<<"\033[00m"<<std::endl;
    return;
  }

  //-------------------------------------------------
  void SimWireMicroBooNE::beginJob() 
  { 
    std::cout<<"\033[95m"<<"beginJob Start!!!"<<"\033[00m"<<std::endl;    

    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;

    fNoiseDist  = tfs->make<TH1D>("Noise", ";Noise  (ADC);", 1000,   -10., 10.);

    art::ServiceHandle<util::LArFFT> fFFT;
    fNTicks = fFFT->FFTSize();

   if ( fNTicks%2 != 0 ) 
      LOG_DEBUG("SimWireMicroBooNE") << "Warning: FFTSize not a power of 2. "
				     << "May cause issues in (de)convolution.\n";

    if ( (int)fNTimeSamples > fNTicks ) 
      mf::LogError("SimWireMircoBooNE") << "Cannot have number of readout samples "
					<< "greater than FFTSize!";
    
    fChargeWork.resize(fNTicks, 0.);
    art::ServiceHandle<geo::Geometry> geo;

    std::cout<<"\033[95m"<<"beginJob End!!!"<<"\033[00m"<<std::endl;    
    return;

  }

  //-------------------------------------------------
  void SimWireMicroBooNE::endJob() 
  {}

  void SimWireMicroBooNE::produce(art::Event& evt)
  {
    std::cout<<"\033[95m"<<"produce Start!!!"<<"\033[00m"<<std::endl;

    art::ServiceHandle<util::TimeService> ts;
    int start_tdc = ts->TPCTick2TDC(0);

    ts->debugReport();
    // get the geometry to be able to figure out signal types and chan -> plane mappings
    art::ServiceHandle<geo::Geometry> geo;
    //unsigned int signalSize = fNTicks;

    std::cout<<"\033[95m"<<"Start TDC @ "<<start_tdc<<"\033[00m"<<std::endl;
    
    std::vector<const sim::SimChannel*> chanHandle;
    evt.getView(fDriftEModuleLabel,chanHandle);

    //Get fIndShape and fColShape from SignalShapingService, on the fly
    art::ServiceHandle<util::SignalShapingServiceMicroBooNE> sss;

    // make a vector of const sim::SimChannel* that has same number
    // of entries as the number of channels in the detector
    // and set the entries for the channels that have signal on them
    // using the chanHandle
    std::vector<const sim::SimChannel*> channels(geo->Nchannels());
    for(size_t c = 0; c < chanHandle.size(); ++c){
      channels[chanHandle[c]->Channel()] = chanHandle[c];
    }
    
    const auto NChannels = geo->Nchannels();

    // vectors for working
    std::vector<short>    adcvec(fNTimeSamples, 0);
    // make a unique_ptr of sim::SimDigits that allows ownership of the produced
    // digits to be transferred to the art::Event after the put statement below
    std::unique_ptr< std::vector<raw::RawDigit>> digcol(new std::vector<raw::RawDigit>);
    digcol->reserve(NChannels);

    unsigned int chan = 0; 
    fChargeWork.clear();
    fChargeWork.resize(fNTicks, 0.);
    
    art::ServiceHandle<util::LArFFT> fFFT;

    //LOOP OVER ALL CHANNELS
    std::map<int,double>::iterator mapIter;      
    for(chan = 0; chan < geo->Nchannels(); chan++) {

      // No need to clear fChargeWork as all contents will be over-written
      //fChargeWork.clear();
      //fChargeWork.resize(fNTicks,0);

      // get the sim::SimChannel for this channel
      const sim::SimChannel* sc = channels[chan];

      if( sc ){

	// loop over the tdcs and grab the number of electrons for each
	for(int t = start_tdc; t < (int)(fChargeWork.size()); ++t) {

	  // continue if tdc < 0
	  if( t < 0 ) continue;

	  fChargeWork[t] = sc->Charge(t);

	  if(sc->Charge(t))

	    std::cout<<chan<<" : "<<t<<" : "<<sc->Charge(t)<<std::endl;

	}

	std::cout<<"Finished reading Ch "<<chan<<std::endl;
	std::cout<<"Size comparison: "<<fChargeWork.size()<<" ?= "<<fNTicks<<std::endl; 
	//fChargeWork.resize(fNTicks,0);

	std::cout<<"Running convolution starting!"<<std::endl;
	
        // Convolve charge with appropriate response function 
	sss->Convolute(chan,fChargeWork);

	std::cout<<"Running convolution finished!"<<std::endl;
      }
      
      //Generate Noise:
      std::vector<float> noisetmp;
      if (fGenNoise){
	if (fGenNoiseInTime)
	  noisetmp = GenNoiseInTime();
	else
	  noisetmp = GenNoiseInFreq();
      }
      else{
	noisetmp.clear();
	noisetmp.resize(fNTicks, 0.);
      }

      //Pedestal determination
      float ped_mean = fCollectionPed;
      geo::SigType_t sigtype = geo->SignalType(chan);
      if (sigtype == geo::kInduction)
        ped_mean = fInductionPed;
      else if (sigtype == geo::kCollection)
        ped_mean = fCollectionPed;
      //slight variation on ped on order of RMS of baselien variation
      art::ServiceHandle<art::RandomNumberGenerator> rng;
      CLHEP::HepRandomEngine &engine = rng->getEngine();
      CLHEP::RandGaussQ rGaussPed(engine, 0.0, fBaselineRMS);
      ped_mean += rGaussPed.fire();
      
      for(unsigned int i = 0; i < fNTimeSamples; ++i){
 	float adcval = noisetmp[i] + fChargeWork[i] + ped_mean;
	
	//Add Noise to NoiseDist Histogram
	if (i%100==0)
	  fNoiseDist->Fill(noisetmp[i]);

	//allow for ADC saturation
	if ( adcval > adcsaturation )
	  adcval = adcsaturation;
	//don't allow for "negative" saturation
	if ( adcval < 0 )
	  adcval = 0;

	adcvec[i] = (unsigned short)(adcval);

      }// end loop over signal size

      // resize the adcvec to be the correct number of time samples, 
      // just drop the extra samples
      //adcvec.resize(fNTimeSamples);
      
      // compress the adc vector using the desired compression scheme,
      // if raw::kNone is selected nothing happens to adcvec
      // This shrinks adcvec, if fCompression is not kNone.
      raw::Compress(adcvec, fCompression); 
      
      // add this digit to the collection
      raw::RawDigit rd(chan, fNTimeSamples, adcvec, fCompression);
      rd.SetPedestal(ped_mean);
      digcol->push_back(rd);

    }// end loop over channels      

    evt.put(std::move(digcol));
    std::cout<<"\033[95m"<<"produce End!!!"<<"\033[00m"<<std::endl;    
    return;
  }
  
  //-------------------------------------------------                                                                                                 
  std::vector<float> SimWireMicroBooNE::GenNoiseInTime()
  {
    std::cout<<"\033[95m"<<"GenNoiseInTime Start!!!"<<"\033[00m"<<std::endl;
    //ART random number service                                                                                                                       
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandGaussQ rGauss(engine, 0.0, fNoiseFact);

    std::vector<float> noise;

    noise.clear();
    noise.resize(fNTicks, 0.);
    //In this case fNoiseFact is a value in ADC counts
    //It is going to be the Noise RMS
    //loop over all bins in "noise" vector
    //and insert random noise value
    for (unsigned int i=0; i<noise.size(); i++)
      noise.at(i) = rGauss.fire();
    std::cout<<"\033[95m"<<"GenNoiseInTime End!!!"<<"\033[00m"<<std::endl;
    return noise;
  }


  //-------------------------------------------------
  std::vector<float> SimWireMicroBooNE::GenNoiseInFreq()
  {
    std::cout<<"\033[95m"<<"GenNoiseInFreq Start!!!"<<"\033[00m"<<std::endl;    
    std::vector<float> noise;
    
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandFlat flat(engine,-1,1);

    noise.clear();
    noise.resize(fNTicks, 0.);
    // noise in frequency space
    std::vector<TComplex> noiseFrequency(fNTicks/2+1, 0.);

    double pval = 0.; 
    double lofilter = 0.;
    double phase = 0.;
    double rnd[2] = {0.};

    // width of frequencyBin in kHz
    double binWidth = 1.0/(fNTicks*fSampleRate*1.0e-6);

    for(int i=0; i< fNTicks/2+1; ++i){
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
      noiseFrequency[i] += tc;
    }
    
    
    // mf::LogInfo("SimWireMicroBooNE") << "filled noise freq";
    
    // inverse FFT MCSignal
    art::ServiceHandle<util::LArFFT> fFFT;
    fFFT->DoInvFFT(noiseFrequency, noise);
    
    noiseFrequency.clear();
    
    // multiply each noise value by fNTicks as the InvFFT 
    // divides each bin by fNTicks assuming that a forward FFT
    // has already been done.
    for(unsigned int i = 0; i < noise.size(); ++i) noise[i] *= 1.*fNTicks;

    std::cout<<"\033[95m"<<"GenNoiseInFreq End!!!"<<"\033[00m"<<std::endl;    
    return noise;
    
  }
  
  
}
  
