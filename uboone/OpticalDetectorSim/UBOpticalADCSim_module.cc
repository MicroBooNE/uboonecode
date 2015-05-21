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

#ifndef UBOpticalADCSim_module_CC
#define UBOpticalADCSim_module_CC

// framework
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

/// LArSoft
#include "OpticalDetectorData/ChannelData.h" // from lardata
#include "Geometry/Geometry.h" // larcore
#include "Simulation/SimPhotons.h" // larsim
#include "UBOpticalADC.h" // uboonecode
#include "UBLogicPulseADC.h" // uboonecode
#include "uboone/Geometry/ChannelMapUBooNEAlg.h" // uboone

/// nutools
#include "Simulation/BeamGateInfo.h"

namespace opdet {
  /**
     \class UBOpticalADCSim
  */ 
  class UBOpticalADCSim : public art::EDProducer{
  public:

    /// Default ctor
    UBOpticalADCSim(const fhicl::ParameterSet&);

    /// Default dtor
    ~UBOpticalADCSim();

    /// art::EDProducer::produce implementation
    virtual void produce (art::Event&); 

  protected:
    
    /// G4 photons producer module name 
    std::string fG4ModName;

    /// BeamGateInfo producer module name
    std::vector<std::string> fBeamModName;

    /// OpticalADC processor class instance
    UBOpticalADC fOpticalGen;

    /// LogicPulseADC processor class instance
    UBLogicPulseADC fLogicGen;

    /// Duration of waveform in ns to be generated as digitizer output
    double fDuration;

    /// G4 time to start waveform generation (default 0)
    double fG4StartTime;

  };

} 

#endif

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
namespace opdet {
  DEFINE_ART_MODULE(UBOpticalADCSim)
}

namespace opdet {

  //###############################################################
  UBOpticalADCSim::UBOpticalADCSim(fhicl::ParameterSet const& pset)
  //###############################################################
  {
    fG4ModName = pset.get<std::string>("G4ModName");

    fBeamModName = pset.get<std::vector<std::string> >("BeamModName");

    if(pset.get<bool>("EnableSpread")) fOpticalGen.EnableSpread(true);
    else fOpticalGen.EnableSpread(false);

    fDuration = pset.get<double>("Duration");

    fG4StartTime = pset.get<double>("G4StartTime",0);

    fLogicGen.SetAmplitude(pset.get<double>("LogicPulseAmplitude"));

    fLogicGen.SetPedestal(2048,0.3); // move to ubchannelconfig

    produces< std::vector<optdata::ChannelData> >();
  }

  //#################################
  UBOpticalADCSim::~UBOpticalADCSim()
  //#################################
  {}
   
  //############################################
  void UBOpticalADCSim::produce(art::Event& evt) 
  //############################################
  {

    // 
    // Get services
    //
    ::art::ServiceHandle<geo::Geometry> geom;
    ::art::ServiceHandle<util::TimeService> ts;

    // allocate the container
    ::std::unique_ptr< std::vector<optdata::ChannelData> > wfs(new ::std::vector<optdata::ChannelData>);

    // get the clock definition
    ::util::ElecClock clock = ts->OpticalClock();
    clock.SetTime(ts->G4ToElecTime(fG4StartTime));
    if(clock.Time()<0)
      throw UBOpticalException(Form("Cannot start readout @ %g (before Electronics Clock start time %g)",
				    fG4StartTime,
				    (-1)*(ts->G4ToElecTime(0))
				    )
			       );
    fOpticalGen.SetTimeInfo(clock,fDuration);
    fLogicGen.SetTimeInfo(clock,fDuration);

    //
    // Read-in data
    //
    art::Handle< std::vector<sim::SimPhotons> > pmtHandle;
    evt.getByLabel(fG4ModName, pmtHandle);
    if(!pmtHandle.isValid()) {
      std::cout << Form("Did not find any G4 photons from a prodcuer: %s",fG4ModName.c_str()) << std::endl;
      return;
    }

    // From larsim/Simulation/SimListUtils.cxx, pmtHandle possibly contain same PMT channel number more than once.
    // Create a dumb map (vector of vector) that can be used to identify indexes of pmtHandle that
    // corresponds to a specific PMT to avoid loop over the entire pmt Handle per channel.
    // this assumes consecutiely numbered opdets, which is OK, I think. Would using a map future-proof this?
    std::vector<std::vector<size_t> > pmt_indexes(geom->NOpDets(),std::vector<size_t>());
    for(auto &v : pmt_indexes) v.reserve(10); // reserve at least 10 (cost nothing in memory but help speed)

    for(size_t i=0; i<pmtHandle->size(); ++i) {
      const art::Ptr<sim::SimPhotons> pmt(pmtHandle,i);
      if(pmt->OpChannel() >= (int)geom->NOpDets())
	throw UBOpticalException(Form("Found OpChannel (%d) larger than # channels from geo (%zu)!",
				      pmt->OpChannel(),
				      pmt_indexes.size())
				 );
      pmt_indexes.at(pmt->OpChannel()).push_back(i);
    }

    //
    // Loop over opdets, and channels, process photons, make, then store waveforms
    //
    for(unsigned int ipmt=0; ipmt<geom->NOpDets(); ipmt++) {
      
      // transfer the time of each hit (in opdet 'ch') into a vector<double>
      std::vector<double> photon_time;
      for(auto const &pmt_index : pmt_indexes.at(ipmt)) {
	
	const art::Ptr<sim::SimPhotons> pmt_ptr(pmtHandle,pmt_index);

	photon_time.reserve(photon_time.size() + pmt_ptr->size());
	
	for(size_t photon_index=0; photon_index<pmt_ptr->size(); ++photon_index)
	  photon_time.push_back(pmt_ptr->at(photon_index).Time);
	
      }

      // send the hits over to the waveform generator
      fOpticalGen.SetPhotons(photon_time);
      
      // generate dark noise hits
      fOpticalGen.GenDarkNoise(ipmt,fG4StartTime);
      
      /// for each pmt we generate multiple readout streams with different gains
      for (unsigned int ireadout=0; ireadout<geom->NOpHardwareChannels(ipmt); ireadout++) {
	unsigned int channel_num = geom->OpChannel( ipmt, ireadout );
	optdata::ChannelData adc_wfm( channel_num );
	fOpticalGen.GenWaveform(ipmt, adc_wfm );
	wfs->emplace_back( adc_wfm );
      }	
    } // loop over pmts
    
    //
    // Handle special readout channels (>= 40)
    //
    std::shared_ptr< const geo::ChannelMapUBooNEAlg > chanmap = std::dynamic_pointer_cast< const geo::ChannelMapUBooNEAlg >( geom->GetChannelMapAlg() );
    art::ServiceHandle<opdet::UBOpticalChConfig> ch_conf;
    std::vector< unsigned int > logicchannels;
    chanmap->GetLogicChannelList( logicchannels );
    
    //for(size_t ch=kLogicStartChannel; ch<(kLogicStartChannel+kLogicNChannel); ++ch) {
    for ( auto logicch : logicchannels ) {
      unsigned int ch = logicch; 
      opdet::UBOpticalChannelCategory_t chcat = chanmap->GetChannelType( ch );

      fLogicGen.SetPedestal( ch_conf->GetParameter( kPedestalMean, ch ), ch_conf->GetParameter( kPedestalSpread, ch ) );

      if( (chcat == opdet::FEMBeamTriggerBNB || chcat == opdet::FEMBeamTriggerNUMI) && fBeamModName.size() ) {

	// get beam gate data product
	for(auto const& name : fBeamModName) {
	  art::Handle< std::vector<sim::BeamGateInfo> > beamHandle;
	  evt.getByLabel(name, beamHandle);
	  if(!beamHandle.isValid()) continue;
	  
	  for(size_t i=0; i<beamHandle->size(); ++i) {

	    const art::Ptr<sim::BeamGateInfo> beam_ptr(beamHandle,i);
	    
	    if( (beam_ptr->BeamType() == ::sim::kBNB  && chcat == opdet::FEMBeamTriggerBNB ) ||
		(beam_ptr->BeamType() == ::sim::kNuMI && chcat == opdet::FEMBeamTriggerNUMI ) )

	      fLogicGen.AddPulse(beam_ptr->Start());

	  }

	}

      }
      
      optdata::ChannelData logic_wfm( ch );
      fLogicGen.GenWaveform( logic_wfm );
      wfs->emplace_back( logic_wfm );
      
    }//loop over logic channels
    
    //
    // Finally, pass to art event
    // 
    if(wfs->size())
      evt.put(std::move(wfs));

    // Make sure to free memory
    fOpticalGen.Reset();
    fLogicGen.Reset();
    
  }
} 
/** @} */ // end of doxygen group 

