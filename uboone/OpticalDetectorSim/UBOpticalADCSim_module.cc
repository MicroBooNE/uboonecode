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
#include "lardata/OpticalDetectorData/ChannelDataGroup.h" // from lardata
#include "lardata/DetectorInfoServices/DetectorClocksService.h" // lardata
#include "larcore/Geometry/Geometry.h" // larcore
#include "larsim/Simulation/SimPhotons.h" // larsim
#include "UBOpticalADC.h" // uboonecode
#include "UBLogicPulseADC.h" // uboonecode
#include "uboone/Geometry/UBOpReadoutMap.h" // uboone

/// nutools
#include "larsim/Simulation/BeamGateInfo.h"

#include <algorithm>

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

    /// User-defined beamgate pulse (BNB) in G4 time ns
    std::vector<double> fUserBNBTime_v;

    /// User-defined beamgate pulse (NuMI)
    std::vector<double> fUserNuMITime_v;

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

    fUserBNBTime_v = pset.get<std::vector<double> >("UserBNBTime");

    fUserNuMITime_v = pset.get<std::vector<double> >("UserNuMITime");

    produces< optdata::ChannelDataGroup >();

    produces< std::vector<sim::BeamGateInfo > >();
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
    auto const* ts = lar::providerFrom<detinfo::DetectorClocksService>();

    // allocate the container
    ::std::unique_ptr< optdata::ChannelDataGroup > wfs(new optdata::ChannelDataGroup);
    ::std::unique_ptr< std::vector<sim::BeamGateInfo > > beam_info_ptr(new std::vector<sim::BeamGateInfo>);

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

    wfs->reserve(geom->NOpChannels());
    wfs->SetFrame(clock.Frame());
    wfs->SetTimeSlice(clock.Sample());

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

    // ================================================================================================================================
    //
    // Loop over opdets, and make waveforms for each readout channels. process photons, make, then store waveforms
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
      //std::cout<<"INPUT: pmt="<<ipmt<<" PE ="<<photon_time.size()<<std::endl;
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
	//fOpticalGen.GenWaveform(ipmt, adc_wfm );
	fOpticalGen.GenWaveform(channel_num, adc_wfm );
	/*
	if(channel_num<32) {
	  double pe=0;
	  double ped_mean=0;
	  double max=0;
	  double min=1e12;
	  for(auto const& v : adc_wfm) {
	    ped_mean += v;
	    if(v > max) max = (double)v;
	    if(v < min) min = (double)v;
	  }
	  ped_mean /= ((double)(adc_wfm.size()));
	  for(auto const& v : adc_wfm)
	    if(v > ped_mean + 3) pe += (v - ped_mean);

	  std::cout<<"Digit Ch: "<<channel_num
		   <<" Pedestal: "<<ped_mean
		   <<" ... " << min << " => "<<max
		   <<" PE: "<< pe / 119.76 << " or PE: " << (max - ped_mean)/20. << std::endl<<std::endl;
	}
	*/
	wfs->emplace_back( adc_wfm );
      }	
    } // loop over pmts
    
    // ================================================================================================================================
    //
    // Handle special readout channels (>= 40)
    //
    art::ServiceHandle<geo::UBOpReadoutMap> chanmap;
    art::ServiceHandle<opdet::UBOpticalChConfig> ch_conf;
    std::vector< unsigned int > logicchannels;
    chanmap->GetLogicChannelList( logicchannels );
    
    //for(size_t ch=kLogicStartChannel; ch<(kLogicStartChannel+kLogicNChannel); ++ch) {
    for ( auto logicch : logicchannels ) {
      unsigned int ch = logicch; 
      opdet::UBOpticalChannelCategory_t chcat = chanmap->GetChannelCategory( ch );

      fLogicGen.SetPedestal( ch_conf->GetFloat( kPedestalMean, ch ), ch_conf->GetFloat( kPedestalSpread, ch ) );

      if( chcat == opdet::BNBLogicPulse || chcat == opdet::NUMILogicPulse ) {

	if(fBeamModName.size()) {
	  // get beam gate data product
	  for(auto const& name : fBeamModName) {
	    art::Handle< std::vector<sim::BeamGateInfo> > beamHandle;
	    evt.getByLabel(name, beamHandle);
	    if(!beamHandle.isValid()) continue;
	    
	    for(size_t i=0; i<beamHandle->size(); ++i) {
	      
	      const art::Ptr<sim::BeamGateInfo> beam_ptr(beamHandle,i);
	      
	      if( (beam_ptr->BeamType() == ::sim::kBNB  && chcat == opdet::BNBLogicPulse ) ||
		  (beam_ptr->BeamType() == ::sim::kNuMI && chcat == opdet::NUMILogicPulse ) )
		
		fLogicGen.AddPulse(beam_ptr->Start());
	      
	    }
	  }
	}

	// open user-defined beamgate open (BNB)
	if(chcat == opdet::BNBLogicPulse) {
	  for(auto const& t : fUserBNBTime_v)
	    fLogicGen.AddPulse(t);
	}

	if(chcat == opdet::NUMILogicPulse) {

	  for(auto const& t : fUserNuMITime_v)
	    fLogicGen.AddPulse(t);
	}
	// open user-defined beamgate open (NuMI)

      }
      
      optdata::ChannelData logic_wfm( ch );
      fLogicGen.GenWaveform( logic_wfm );
      wfs->emplace_back( logic_wfm );
      
    }//loop over logic channels

    // ================================================================================================================================    
    //
    // Finally, pass to art event
    // 
    if(wfs->size())
      evt.put(std::move(wfs));

    for(auto const& t : fUserBNBTime_v)    
      beam_info_ptr->push_back(sim::BeamGateInfo( t, 1600, ::sim::kBNB) );
    for(auto const& t : fUserNuMITime_v)
      beam_info_ptr->push_back(sim::BeamGateInfo( t, 1600*6, ::sim::kNuMI) );
    evt.put(std::move(beam_info_ptr));
    // Make sure to free memory
    fOpticalGen.Reset();
    fLogicGen.Reset();
    
  }
} 
/** @} */ // end of doxygen group 

