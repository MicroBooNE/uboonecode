/**
 * \file UBTriggerSim_module.cc
 *
 * \ingroup OpticalDetectorSim
 * 
 * \brief optdata::ChannelData producer module for microboone
 *
 * The module takes in art data products, SimPhotons and BeamGateInfo
 * It then produces an art data product consisting of a vector of ChannelData instances,
 * which represents raw ADC waveforms for each readout channel.
 * 
 * the classes UBOpticalADC and UBLogicPulseADC are processors
 * that produce the adc waveforms for this class.
 *
 * @author kazuhiro
 * @author twongjirad
 */

/** \addtogroup OpticalDetectorSim

@{*/

#ifndef UBFlasherMC_module_CC
#define UBFlasherMC_module_CC

//#define _UBFlasherMC_DEBUG_

// framework
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

/// LArSoft
#include "lardataobj/OpticalDetectorData/ChannelDataGroup.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsimobj/Simulation/SimPhotons.h"
#include "UBOpticalADC.h" // uboonecode
#include "UBLogicPulseADC.h" // uboonecode
#include "uboone/Geometry/UBOpReadoutMap.h" // uboone

/// nutools
#include "larsimobj/Simulation/BeamGateInfo.h"

// Other
#include "TF1.h"
#include "TH1S.h"
#include "TRandom3.h"
#include <algorithm>
#include <map>

namespace opdet {
  /**
     \class UBFlasherMC
  */ 
  class UBFlasherMC : public art::EDProducer{
  public:

    /// Default ctor
    UBFlasherMC(const fhicl::ParameterSet&);

    /// Default dtor
    ~UBFlasherMC();

    /// art::EDProducer::produce implementation
    virtual void produce (art::Event&); 

  protected:

    void setupFlasher( fhicl::ParameterSet const& pset);

    typedef enum { kBurst, kSequence } FlasherRunMode_t;

    /// OpticalADC processor class instance
    UBOpticalADC fOpticalGen;

    /// LogicPulseADC processor class instance
    UBLogicPulseADC fLogicGen;

    /// Duration of waveform in microseconds to be generated as digitizer output
    double fDuration;

    /// G4 time to start waveform generation (default 0)
    double fG4StartTime;

    // FLASHER FICHLE PARAMETERS
    unsigned int fTTLChannel; //<  Readout Channel where TTL pulse from Flasher will go
    unsigned int fTTLMod36Channel; //<  Readout Channel where TTL pulse from Flasher will go

    std::vector< std::string > fLEDlevels; //< Hex strings that tell us the level of the flasher

    std::map< unsigned int, double > fPElevels;

    double fHexLevelToPE;

    double fPulserRateMHz; //< rate in MHz that pulser fires

    double fPulserEfficiency; //< fraction of successful pulses (not implemented yet)

    double fPulserDelay_ns;

    bool fMakeDebugPlots;

    // FLASHER VARIABLES
    TRandom3* fRand;

    int fNumberOfPulseTrains;

    FlasherRunMode_t fMode;

    // GATE VARIABLES
    double fGlobalTimeOffset;             /// The start of a simulated "beam gate".
    double fRandomTimeOffset;             /// The width of a simulated "beam gate".
    ::sim::BeamType_t fBeamType;          /// The type of beam


  };

} 

#endif

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
namespace opdet {
  DEFINE_ART_MODULE(UBFlasherMC)
}

namespace opdet {

  //###############################################################

  UBFlasherMC::UBFlasherMC(fhicl::ParameterSet const& pset)
  {

    if(pset.get<bool>("EnableSpread")) fOpticalGen.EnableSpread(true);
    else fOpticalGen.EnableSpread(false);

    fDuration = pset.get<double>("Duration");

    fG4StartTime = pset.get<double>("G4StartTime",0);

    fLogicGen.SetAmplitude(pset.get<double>("LogicPulseAmplitude"));

    fLogicGen.SetPedestal(2048,0.3); // move to ubchannelconfig

    setupFlasher( pset );

    fRand = new TRandom3(10);

    // beam config
    fBeamType = ::sim::kBNB;
    fGlobalTimeOffset = pset.get< double >("GlobalTimeOffset",1000.0);
    fRandomTimeOffset = pset.get< double >("RandomTimeOffset",1600.0);
    
    produces< optdata::ChannelDataGroup >();
    produces< std::vector<sim::BeamGateInfo> >();
  }

  //#################################

  UBFlasherMC::~UBFlasherMC()
  {}


  //############################################

  void UBFlasherMC::setupFlasher( fhicl::ParameterSet const& pset) {

    // services
    art::ServiceHandle<geo::Geometry> geom;

    fTTLChannel = pset.get<int>("TTLChannel");
    fTTLMod36Channel = pset.get<int>("TTLMod36Channel");
    fLEDlevels  = pset.get<std::vector<std::string>>("LEDlevelsHex");
    fPulserDelay_ns = pset.get<double>("PulserDelay_ns",100.0);
    fPulserRateMHz = pset.get<double>("PulserRateMHz",1.0);
    fPulserEfficiency = pset.get<double>("PulserEfficency");
    fHexLevelToPE    = pset.get<double>("HexToPE");
    fMakeDebugPlots  = pset.get<bool>("MakeDebugPlots");
    std::string mode = pset.get<std::string>("Mode");
    std::transform(mode.begin(), mode.end(), mode.begin(), ::tolower);

    if ( mode=="burst" )
      fMode = kBurst;
    else if ( mode=="sequence" )
      fMode = kSequence;
    else {
      throw std::runtime_error("Error in UBFlasherMC::setupFlasher. Unrecognized run mode."); 
    }

    // Set pe levels
    fPElevels.clear();
    for ( unsigned int iopdet=0; iopdet<fLEDlevels.size(); iopdet++ ) {
      std::string hex = fLEDlevels.at(iopdet);
      unsigned int dec = 0;
      sscanf(hex.c_str(), "%x", &dec); 
      double level = (double)dec;
      fPElevels[ iopdet ]  = level*fHexLevelToPE;
    }
    
  }

   
  //############################################

  void UBFlasherMC::produce(art::Event& evt) 
  {

    // 
    // Get services
    //
    ::art::ServiceHandle<geo::Geometry> geom;
    ::art::ServiceHandle<geo::UBOpReadoutMap> chanmap;
    ::art::ServiceHandle<opdet::UBOpticalChConfig> ch_conf;
    auto const* ts = lar::providerFrom<detinfo::DetectorClocksService>();
    ::art::ServiceHandle<art::TFileService> tfs;

    // allocate the container
    ::std::unique_ptr< optdata::ChannelDataGroup > wfs(new optdata::ChannelDataGroup);

    // get the clock definition
    ::detinfo::ElecClock clock = ts->OpticalClock();
    clock.SetTime(ts->G4ToElecTime(fG4StartTime));
    if(clock.Time()<0)
      throw UBOpticalException(Form("Cannot start readout @ %g (before Electronics Clock start time %g)",
				    fG4StartTime,
				    (-1)*(ts->G4ToElecTime(0))
				    )
			       );
    fOpticalGen.SetTimeInfo(clock,fDuration);
    fLogicGen.SetTimeInfo(clock,fDuration);

    fNumberOfPulseTrains = fDuration*(fPulserRateMHz);

    wfs->reserve(geom->NOpChannels());
    wfs->SetFrame(clock.Frame());
    wfs->SetTimeSlice(clock.Sample());

    // ================================================================================================================================
    // Make gate
    // trigger board provides gate which forces a readout of the PMTs
    // for now, we just use the BNB gate
    std::unique_ptr< std::vector<sim::BeamGateInfo> > gateCollection(new std::vector<sim::BeamGateInfo>);
    sim::BeamGateInfo trigger(  fG4StartTime + fGlobalTimeOffset, fRandomTimeOffset, fBeamType );
#ifdef _UBFlasherMC_DEBUG_
    optdata::TimeSlice_t gateTime = ts->OpticalG4Time2TDC(trigger.Start());
    std::cout << "Start of Gate: Time slice=" << gateTime << " G4 time=" << trigger.Start() << std::endl;
#endif
    gateCollection->push_back( trigger );
    
    // ================================================================================================================================
    //
    // Loop over opdets, and make waveforms for each readout channels. process photons, make, then store waveforms
    //
    double period_ns = (1.0/fPulserRateMHz)*1000.0;
    //double led_seq_delay_ns = 100.0;
    double flasher_delay_ns = 125.0;
    
    for(unsigned int ipmt=0; ipmt<geom->NOpDets(); ipmt++) {
      
      // generate times of photon hits from LED

      std::vector<double> photon_time; // microseconds
      photon_time.reserve(10*fNumberOfPulseTrains);

      for (int ipulse=0; ipulse<fNumberOfPulseTrains; ipulse++) {
	
	if ( fMode==kSequence && ipulse%36!=(int)ipmt) 
	  continue;
	
	double pulse_start = fG4StartTime + fGlobalTimeOffset + fPulserDelay_ns; // G4time of gate start
	//if ( fMode==kBurst )
	pulse_start += double(ipulse)*period_ns + flasher_delay_ns;
	//else if ( fMode==kSequence )
	//pulse_start += double(ipulse)*period_ns + double(ipmt)*led_seq_delay_ns +  flasher_delay_ns;
	
	//unsigned int nphotons_in_pulse = fRand->Poisson( fPElevels[ ipmt ] );

	//for(size_t photon_index=0; photon_index<nphotons_in_pulse; ++photon_index) {
	//t = fRand->Gaus( pulse_start, 1.0 );
	//photon_time.push_back( pulse_start );
	//}
      }//end of pulse loop
      
      // send the hits over to the waveform generator
      fOpticalGen.SetPhotons(photon_time);
      
      // generate dark noise hits
      fOpticalGen.GenDarkNoise(ipmt,fG4StartTime);
      
      /// for each pmt we generate multiple readout streams with different gains
      /// tmw: right now adc/pe is assigned per channel. I don't like this. 
      /// this is a shaper property
      for (unsigned int ireadout=0; ireadout<geom->NOpHardwareChannels(ipmt); ireadout++) {
	unsigned int channel_num = geom->OpChannel( ipmt, ireadout );
	optdata::ChannelData adc_wfm( channel_num );
	fOpticalGen.GenWaveform(ipmt, adc_wfm );

	// save if asked to
	
	if ( fMakeDebugPlots  ) {
	  std::ostringstream hname;
	  hname << "flasherwfm_E" << evt.id().event() << "_CH" << (int)channel_num << "_PMT" << (int)ipmt;

	  TH1* wfmHist = tfs->make<TH1S>(hname.str().c_str(),
					 hname.str().c_str(),
					 adc_wfm.size(),
					 0, adc_wfm.size());
	  for ( int ibin=0; ibin<(int)adc_wfm.size(); ibin++ )
	    wfmHist->SetBinContent( ibin+1, adc_wfm.at(ibin) );
	}


	wfs->emplace_back( adc_wfm );
      }	
    } // loop over pmts
    
    // ================================================================================================================================
    //
    // Handle special readout channels (>= 40)
    //
    std::vector< unsigned int > logicchannels;
    chanmap->GetLogicChannelList( logicchannels );
    //logicchannels.push_back( (unsigned int)fTTLChannel );      // the driving pulse from external pulser
    //logicchannels.push_back( (unsigned int)fTTLMod36Channel ); // the return Mod36 pulse
    
    for ( auto logicch : logicchannels ) {
      unsigned int ch = logicch; 
      opdet::UBOpticalChannelCategory_t chcat = chanmap->GetChannelCategory( ch );

      fLogicGen.SetPedestal( ch_conf->GetFloat( kPedestalMean, ch ), ch_conf->GetFloat( kPedestalSpread, ch ) );

      std::vector<double> pulse_times;
      pulse_times.clear();

      // THE FLASHER LOGIC PULSES
      if ( ch==fTTLMod36Channel || ch==fTTLChannel ) {
	// Pulser TTL pulse train into TTLChannel
	// Flasher MOD36 output into TTLMod36Channel
	for (int ipulse=0; ipulse<fNumberOfPulseTrains; ipulse++) {
	  
	  if ( ch==fTTLMod36Channel && ipulse%36!=0)
	    continue;
	  
	  double pulse_start = fG4StartTime + fGlobalTimeOffset + fPulserDelay_ns;
	  pulse_start += ipulse*period_ns;
	  
	  pulse_times.push_back( pulse_start );
	  
	}
      }
      // TRIGGER GATE
      else if ( (chcat == opdet::BNBLogicPulse ) ) {
	pulse_times.push_back( trigger.Start() );
      }//end of if TRIGGER GATE

      fLogicGen.SetPulses( pulse_times );
      
      optdata::ChannelData logic_wfm( ch );
      fLogicGen.GenWaveform( logic_wfm );


      if ( fMakeDebugPlots  ) {
	std::ostringstream hname;
	hname << "flasherwfm_E" << evt.id().event() << "_LOGICCH" << (int)ch;
	
	TH1* wfmHist = tfs->make<TH1S>(hname.str().c_str(),
				       hname.str().c_str(),
				       logic_wfm.size(),
				       0, logic_wfm.size());
	for ( int ibin=0; ibin<(int)logic_wfm.size(); ibin++ )
	  wfmHist->SetBinContent( ibin+1, logic_wfm.at(ibin) );
      }
      
      wfs->emplace_back( logic_wfm );
      
    }//loop over logic channels
    
    // ================================================================================================================================    
    //
    // Finally, pass to art event
    // 
    if(wfs->size())
      evt.put(std::move(wfs));
    evt.put(std::move(gateCollection));

    
    // Make sure to free memory
    fOpticalGen.Reset();
    fLogicGen.Reset();
    
  }
    
}
/** @} */ // end of doxygen group 
