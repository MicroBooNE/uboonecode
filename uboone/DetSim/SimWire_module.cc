////////////////////////////////////////////////////////////////////////
// $Id: SimWire.cxx,v 1.22 2010/04/23 20:30:53 seligman Exp $
//
// SimWire class designed to simulate signal on a wire in the TPC
//
// dcaratelli@nevis.columbia.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef SIMWIRE_H
#define SIMWIRE_H

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

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Simulation/sim.h"
#include "Simulation/SimChannel.h"
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/LArFFT.h"
#include "Utilities/LArProperties.h"

namespace art {
  class Event;
  class ParameterSet;
}

namespace geo { class Geometry; }

///Detector simulation of raw signals on wires
namespace detsim {

  // Base class for creation of raw signals on wires. 
  class SimWire : public art::EDProducer {
    
  public:
        
    explicit SimWire(fhicl::ParameterSet const& pset); 
    virtual ~SimWire();
    
    // read/write access to event
    void produce (art::Event& evt);
    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& p);

  private:

    std::vector<float>         GenNoiseInTime();

    raw::Compress_t        fCompression;      ///< compression type to use
    double                 fSampleRate;       ///< sampling rate in ns
    double                 fNoiseFact;        ///< noise scale factor 
    int                    fNTicks;           ///< number of ticks of the clock
    unsigned int           fNSamplesReadout;  ///< number of ADC readout samples in 1 readout frame
    unsigned int           fPedestal;         ///< pedestal amplitude [ADCs]
    double                 fSigAmp;           ///< signal amplitude [ADCs]
    double                 fSigWidth;         ///< signal wifdth [time ticks]
    int                    fSigTime;          ///< Time at which pulse should happen
    std::string            fSigType;          ///< signal type. For now gaussian or pulse
    TH1D*                  fNoiseDist;        ///< distribution of noise counts
    TH1D*                  fWaveform;         ///< waveform histogram

  }; // class SimWire

}

namespace detsim{

  //-------------------------------------------------
  SimWire::SimWire(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);

    produces< std::vector<raw::RawDigit>   >();
 
    fCompression = raw::kNone;
    std::string compression(pset.get< std::string >("CompressionType"));
    if(compression.compare("Huffman") == 0) fCompression = raw::kHuffman;    

    // get the random number seed, use a random default if not specified    
    // in the configuration file.  
    unsigned int seed = pset.get< unsigned int >("Seed", sim::GetRandomNumberSeed());

    createEngine(seed);
  }

  //-------------------------------------------------
  SimWire::~SimWire()
  { }

  //-------------------------------------------------
  void SimWire::reconfigure(fhicl::ParameterSet const& p) 
  {
   
    fNoiseFact        = p.get< double              >("NoiseFact");
    fNTicks           = p.get< int                 >("NTicks");
    fPedestal         = p.get< unsigned int        >("Pedestal");
    fSigAmp           = p.get< double              >("SigAmp");
    fSigWidth         = p.get< double              >("SigWidth");
    fSigType          = p.get< std::string         >("SigType");
    fSigTime          = p.get< int                 >("SigTime");
    
   
    art::ServiceHandle<util::DetectorProperties> detprop;
    fSampleRate       = detprop->SamplingRate();
    fNSamplesReadout  = detprop->NumberTimeSamples();

    return;
  }

  //-------------------------------------------------
  void SimWire::beginJob() 
  { 
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;

    fNoiseDist      = tfs->make<TH1D>("Noise", ";Noise (ADC);", 1000, -20., 20.);
    fWaveform       = tfs->make<TH1D>("Waveform", "; Pulse [ADC];", fNTicks, 0, fNTicks);

    art::ServiceHandle<util::LArFFT> fFFT;
    fNTicks = fFFT->FFTSize();

    return;

  }

  //-------------------------------------------------
  void SimWire::endJob() 
  {
  }

  //-------------------------------------------------
  void SimWire::produce(art::Event& evt)
  {

    //std::cout << "in SimWire::produce " << std::endl;

    // get the geometry to be able to figure out signal types and chan -> plane mappings
    art::ServiceHandle<geo::Geometry> geo;
    unsigned int signalSize = fNTicks;
    // vectors for working
    std::vector<short>    adcvec(signalSize, 0);	
    std::vector<float>    sigvec(signalSize, 0);
    std::vector<float>    noisevec(signalSize, 0);

    // make an unique_ptr of sim::SimDigits that allows ownership of the produced
    // digits to be transferred to the art::Event after the put statement below
    std::unique_ptr< std::vector<raw::RawDigit>   >  digcol(new std::vector<raw::RawDigit>);
	  
    unsigned int chan = 1; 

    if ( fSigType == "gauss" ){//gaussian
      for (unsigned short i=0; i<signalSize; i++){
	sigvec[i] += fSigAmp*TMath::Gaus(i,fSigTime,fSigWidth,0)*TMath::Sqrt((2*TMath::Pi()));
      }
    }
    if ( fSigType == "pulse" ){//square pulse
      for (int i=0; i<fSigTime; i++){
	unsigned int timetmp = (int)( fSigTime/2.+ i );
	if ( (timetmp < signalSize) && (timetmp > 0) )
	  sigvec[timetmp] += fSigAmp;
      }
    }


    //generate noise
    noisevec = GenNoiseInTime();

    for(unsigned int i = 0; i < signalSize; ++i){
      float adcval = fPedestal + noisevec[i] + sigvec[i];
      fWaveform->SetBinContent(i+1,adcval);
      fNoiseDist->Fill(noisevec[i]);
      adcvec[i] = (unsigned short)(adcval);
    }
      
    raw::RawDigit rd(chan, signalSize, adcvec, fCompression);
    // Then, resize adcvec back to full length!
    adcvec.clear();
    adcvec.resize(signalSize,0.0);
    
    // add this digit to the collection
    digcol->push_back(rd);
    
    
    evt.put(std::move(digcol));
    
    return;
  }

  //----------------------------------------------
  std::vector<float> SimWire::GenNoiseInTime()
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

  DEFINE_ART_MODULE(SimWire)

}

#endif // SIMWIRE_H


