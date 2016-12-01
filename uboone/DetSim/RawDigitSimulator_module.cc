////////////////////////////////////////////////////////////////////////
// $Id: RawDigitSimulator.cxx,v 1.22 2010/04/23 20:30:53 seligman Exp $
//
// RawDigitSimulator class designed to simulate signal on a wire in the TPC
//
// dcaratelli@nevis.columbia.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef RAWDIGITSIMULATOR_H
#define RAWDIGITSIMULATOR_H

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

// ROOT includes
#include <TMath.h>
#include <TH1D.h>
#include <TFile.h>
#include "TComplex.h"
#include "TString.h"
#include "TH2.h"

// C++ includes
#include <algorithm>
#include <sstream>
#include <fstream>
#include <bitset>
#include <vector>
#include <string>


// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/EDProducer.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/search_path.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/LArFFT.h"
#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"

namespace art {
  class Event;
  class ParameterSet;
}

namespace geo { class Geometry; }

///Detector simulation of raw signals on wires
namespace detsim {

  // Base class for creation of raw signals on wires. 
  class RawDigitSimulator : public art::EDProducer {
    
  public:

    enum SignalGenType {
      kGaus,
      kSquare,
      kSignalGenTypeMax
    };
        
    explicit RawDigitSimulator(fhicl::ParameterSet const& pset); 
    virtual ~RawDigitSimulator();
    
    // read/write access to event
    void produce (art::Event& evt);
    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& p);

  private:

    std::vector<float>         GenNoiseInTime();

    raw::Compress_t            fCompression;      ///< compression type to use
    double                     fSampleRate;       ///< sampling rate in ns
    double                     fNoiseFact;        ///< noise scale factor 
    int                        fNTicks;           ///< number of ticks of the clock
    unsigned int               fNSamplesReadout;  ///< number of ADC readout samples in 1 readout frame
    unsigned int               fPedestal;         ///< pedestal amplitude [ADCs]
    std::vector<double>        fSigAmp;           ///< signal amplitude [ADCs] or [# electrons]
    std::vector<double>        fSigWidth;         ///< signal wifdth [time ticks]
    std::vector<int>           fSigTime;          ///< Time at which pulse should happen
    std::vector<unsigned char> fSigType;          ///< signal type. For now gaussian or pulse
    bool                       fSigUnit;          ///< true = ADC, false = electrons
    unsigned int               fChannel;          ///< Channel number
    bool                       fGenNoise;         ///< Boolean to generate noise
    size_t fEventCount; ///< count of event processed
    std::map<double,int> fShapingTimeOrder;
  }; // class RawDigitSimulator

}

namespace detsim{

  //-------------------------------------------------
  RawDigitSimulator::RawDigitSimulator(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);

    produces< std::vector<raw::RawDigit>   >();
 
    fCompression = raw::kNone;
    std::string compression(pset.get< std::string >("CompressionType"));
    if(compression.compare("Huffman") == 0) fCompression = raw::kHuffman;    

    // create a default random engine; obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed"
    art::ServiceHandle<rndm::NuRandomService>()
      ->createEngine(*this, pset, "Seed");

    fEventCount = 0;
  }

  //-------------------------------------------------
  RawDigitSimulator::~RawDigitSimulator()
  { }

  //-------------------------------------------------
  void RawDigitSimulator::reconfigure(fhicl::ParameterSet const& p) 
  {
   
    //fNoiseFact        = p.get< double                     >("NoiseFact");
    fNoiseFact        = 0;
    fGenNoise         = p.get< bool                       >("GenNoise");
    fNTicks           = p.get< int                        >("NTicks");
    fPedestal         = p.get< unsigned int               >("Pedestal");
    fSigAmp           = p.get< std::vector<double>        >("SigAmp");
    fSigWidth         = p.get< std::vector<double>        >("SigWidth");
    fSigType          = p.get< std::vector<unsigned char> >("SigType");
    fSigTime          = p.get< std::vector<int>           >("SigTime");
    fSigUnit          = p.get< bool                       >("SigUnit");
    fChannel          = p.get< unsigned int               >("Channel");
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fSampleRate       = detprop->SamplingRate();
    fNSamplesReadout  = detprop->NumberTimeSamples();
    fShapingTimeOrder = { {0.5, 0}, {1.0, 1}, {2.0, 2}, {3.0, 3} };

    // Simple check:
    if( fSigAmp.size() != fSigWidth.size() ||
	fSigAmp.size() != fSigType.size()  ||
	fSigAmp.size() != fSigTime.size() )

      throw cet::exception(__PRETTY_FUNCTION__) << "Input signal info vector have different length!";

    for(auto const& v : fSigType)
      
      if(v >= kSignalGenTypeMax) 

	throw cet::exception(__PRETTY_FUNCTION__) << "Invalid signal type found!";

    return;
  }

  //-------------------------------------------------
  void RawDigitSimulator::beginJob() 
  { 
    fEventCount = 0;

    art::ServiceHandle<util::LArFFT> fFFT;
    fFFT->ReinitializeFFT(fNTicks,fFFT->FFTOptions(),fFFT->FFTFitBins());
    fNTicks = fFFT->FFTSize();

    return;

  }

  //-------------------------------------------------
  void RawDigitSimulator::endJob() 
  {
  }

  //-------------------------------------------------
  void RawDigitSimulator::produce(art::Event& evt)
  {

    //std::cout << "in SimWire::produce " << std::endl;

    // get the geometry to be able to figure out signal types and chan -> plane mappings
    art::ServiceHandle<geo::Geometry> geo;
    unsigned int signalSize = fNTicks;
    std::cout << "Signal size is " << signalSize << std::endl;
    
    // vectors for working
    std::vector<short>    adcvec(signalSize, 0);	
    std::vector<double>   sigvec(signalSize, 0);
    std::vector<float>    noisevec(signalSize, 0);

    art::ServiceHandle<util::SignalShapingServiceMicroBooNE> sss;
    //Xin remove the time offset, now doing it in the SignalShapingService
    int electron_time_offset = 0.;//sss->FieldResponseTOffset(fChannel);
    std::cout<<"Offset: "<<electron_time_offset << std::endl;

    // make an unique_ptr of sim::SimDigits that allows ownership of the produced
    // digits to be transferred to the art::Event after the put statement below
    std::unique_ptr< std::vector<raw::RawDigit>   >  digcol(new std::vector<raw::RawDigit>);

    //fill all pulses according to their type and properties
    for (unsigned int n=0; n<fSigType.size(); n++){
    
      if ( fSigType.at(n) == kGaus ){//gaussian
	for (unsigned short i=0; i<signalSize; i++){
	  sigvec.at(i) += fSigAmp.at(n)*TMath::Gaus(i,fSigTime.at(n),fSigWidth.at(n),0)*TMath::Sqrt((2*TMath::Pi()));
	}
      }
      if ( fSigType.at(n) == kSquare ){//square pulse
	for (int i=0; i<fSigWidth.at(n); i++){
	  int timetmp = (int)( (fSigTime.at(n) - fSigWidth.at(n)/2.) + i );
	  if(!fSigUnit) timetmp += electron_time_offset;
	  //if(timetmp > sigvec.size()) throw cet::exception(__FUNCTION__) << "Invalid timing: "<<timetmp<<std::endl;
	  if(timetmp < 0) timetmp += sigvec.size();
	  if(timetmp < ((int)signalSize))
	    sigvec.at(timetmp) += fSigAmp.at(n);
	}
      }
      
    }//fill all pulses

    // If the unit is # electrons, run convolution
    if(!fSigUnit) {
      sss->Convolute(fChannel,sigvec);
    }


    //
    // Do noise
    //
    //get ASIC Gain and Noise in ADCatLowestGain:
    double fASICGain      = sss->GetASICGain(fChannel);    //Jyoti - to read different gain for U,V & Y planes 
    double fShapingTime   = sss->GetShapingTime(fChannel); //Jyoti - to read different shaping time for U,V & Y planes 
    //Check that shaping time is an allowed value
    //If so, Pick out noise factor 
    //If not, through exception

    if ( fShapingTimeOrder.find( fShapingTime ) != fShapingTimeOrder.end() ){
      geo::View_t view = geo->View(fChannel);
      auto noiseFactVec = sss->GetNoiseFactVec();

      fNoiseFact = noiseFactVec[(int)view].at( fShapingTimeOrder.find( fShapingTime )->second );
    }
    else{//Throw exception...
      throw cet::exception("RawDigitSimulator_module")
	<< "\033[93m"
	<< "Shaping Time received from signalservices_microboone.fcl is not one of allowed values"
	<< std::endl
	<< "Allowed values: 0.5, 1.0, 2.0, 3.0 usec"
	<< "\033[00m"
	<< std::endl;
    }
    //Take into account ASIC Gain
    fNoiseFact *= fASICGain/4.7;
    //generate noise
    if(fGenNoise)
      noisevec = GenNoiseInTime();

    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
    TH1D* fNoiseDist = tfs->make<TH1D>(Form("Noise_%04zu",fEventCount), ";Noise (ADC);", fNTicks,-0.5,fNTicks-0.5);
    TH1D* fSignalDist = tfs->make<TH1D>(Form("Signal_%04zu",fEventCount), ";Noise (ADC);", fNTicks,-0.5,fNTicks-0.5);
    TH1D* fWaveform  = tfs->make<TH1D>(Form("Waveform_%04zu",fEventCount), "; Pulse [ADC];", fNTicks, -0.5, fNTicks-0.5);

    for(unsigned int i = 0; i < signalSize; ++i){
      float adcval = fPedestal + noisevec.at(i) + sigvec.at(i);
      fWaveform->SetBinContent(i+1,(unsigned int)(adcval));
      fSignalDist->Fill(sigvec.at(i));
      fNoiseDist->Fill(noisevec.at(i));
      adcvec.at(i) = (unsigned short)(adcval);
    }
      
    raw::RawDigit rd(fChannel, signalSize, adcvec, fCompression);
    rd.SetPedestal(fPedestal);

    // Then, resize adcvec back to full length!
    adcvec.clear();
    adcvec.resize(signalSize,0.0);
    
    // add this digit to the collection
    digcol->push_back(rd);
    
    
    evt.put(std::move(digcol));

    fEventCount++;

    return;
  }

  //----------------------------------------------
  std::vector<float> RawDigitSimulator::GenNoiseInTime()
  {

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

    return noise;
  }


}


namespace detsim{

  DEFINE_ART_MODULE(RawDigitSimulator)

}

#endif // RAWDIGITSIMULATOR_H
