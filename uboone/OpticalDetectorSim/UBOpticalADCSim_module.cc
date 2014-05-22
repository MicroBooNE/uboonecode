/**
 * \file UBTriggerSim_module.cc
 *
 * \ingroup OpticalDetectorSim
 * 
 * \brief optdata::ChannelData producer module for microboone
 *
 * @author kazuhiro
 */

/** \addtogroup OpticalDetectorSim

@{*/

#ifndef UBOpticalADCSim_module_CC
#define UBOpticalADCSim_module_CC

// LArSoft includes
//#include "Geometry/Geometry.h"

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
#include "UBOpticalADC.h"
#include "UBLogicPulseADC.h"
#include "OpticalDetectorData/ChannelDataGroup.h"
#include "Simulation/SimPhotons.h"
#include "Geometry/Geometry.h"

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

    fLogicGen.SetPedestal(2048,0.3);

    produces< std::vector<optdata::ChannelDataGroup> >();
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
    // Initialize
    //
    ::art::ServiceHandle<geo::Geometry> geom;
    ::art::ServiceHandle<util::TimeService> ts;
    ::std::unique_ptr< std::vector<optdata::ChannelDataGroup> > wfs(new ::std::vector<optdata::ChannelDataGroup>);
    ::util::ElecClock clock = ts->OpticalClock();
    clock.SetTime(ts->G4ToElecTime(fG4StartTime));

    fOpticalGen.SetTimeInfo(clock,fDuration);
    fLogicGen.SetTimeInfo(clock,fDuration);

    //
    // Prepare output container
    //
    wfs->reserve(3);

    const size_t hg_index = wfs->size();
    wfs->push_back(optdata::ChannelDataGroup(::optdata::kHighGain,clock.Sample(),clock.Frame())); // high gain

    const size_t lg_index = wfs->size();
    wfs->push_back(optdata::ChannelDataGroup(::optdata::kLowGain,clock.Sample(),clock.Frame())); // low gain

    const size_t logic_index = wfs->size();
    wfs->push_back(optdata::ChannelDataGroup(::optdata::kLogicPulse,clock.Sample(),clock.Frame())); // logic pulse channels

    // Create entries ... # PMTs + 2 channels (NuMI & BNB readout channels)
    wfs->at(hg_index).reserve(geom->NOpChannels());
    wfs->at(lg_index).reserve(geom->NOpChannels());
    wfs->at(logic_index).reserve(kLogicNChannel);

    //
    // Read-in data
    //
    art::Handle< std::vector<sim::SimPhotons> > pmtHandle;
    evt.getByLabel(fG4ModName, pmtHandle);
    if(!pmtHandle.isValid())
      throw UBOpticalException(Form("Did not find any G4 photons from a prodcuer: %s",fG4ModName.c_str()));

    //std::vector<sim::SimPhotons> const& pmts(*pmtHandle);

    // From larsim/Simulation/SimListUtils.cxx, pmtHandle possibly contain same PMT channel number more than once.
    // Create a dumb map (vector of vector) that can be used to identify indexes of pmtHandle that
    // corresponds to a specific PMT to avoid loop over the entire pmt Handle per channel.
    std::vector<std::vector<size_t> > pmt_indexes(geom->NOpChannels(),std::vector<size_t>());
    for(auto &v : pmt_indexes) v.reserve(10); // reserve at least 10 (cost nothing in memory but help speed)

    for(size_t i=0; i<pmtHandle->size(); ++i) {
      const art::Ptr<sim::SimPhotons> pmt(pmtHandle,i);
      if(pmt->OpChannel() > (int)(pmt_indexes.size()))
	throw UBOpticalException(Form("Found OpChannel (%d) larger than # channels from geo (%zu)!",
				      pmt->OpChannel(),
				      pmt_indexes.size())
				 );
      pmt_indexes.at(pmt->OpChannel()).push_back(i);
    }

    //
    // Loop over channels, process photons, & store waveform 
    //
    for(size_t ch=0; ch<geom->NOpChannels(); ++ch) {
      
      wfs->at(hg_index).push_back(optdata::ChannelData(ch));
      
      wfs->at(lg_index).push_back(optdata::ChannelData(ch));

      std::vector<double> photon_time;
      for(auto const &pmt_index : pmt_indexes.at(ch)) {
	
	const art::Ptr<sim::SimPhotons> pmt_ptr(pmtHandle,pmt_index);

	photon_time.reserve(photon_time.size() + pmt_ptr->size());
	
	for(size_t photon_index=0; photon_index<pmt_ptr->size(); ++photon_index)

	  photon_time.push_back(pmt_ptr->at(photon_index).Time);

      }

      fOpticalGen.SetPhotons(photon_time);
      fOpticalGen.GenWaveform(ch,
		       wfs->at(hg_index).at(ch),
		       wfs->at(lg_index).at(ch));

    }

    //
    // Handle special readout channels (>= 40)
    //
    for(size_t ch=kLogicStartChannel; ch<(kLogicStartChannel+kLogicNChannel); ++ch) {

      wfs->at(logic_index).push_back(optdata::ChannelData(ch));

      if( (ch == kFEMChannelBNB || ch == kFEMChannelNuMI) && fBeamModName.size() ) {

	for(auto const& name : fBeamModName) {
	  art::Handle< std::vector<sim::BeamGateInfo> > beamHandle;
	  evt.getByLabel(name, beamHandle);
	  if(!beamHandle.isValid())
	    throw UBOpticalException(Form("Did not find any BeamGateInfo from a prodcuer: %s",name.c_str()));
	  
	  for(size_t i=0; i<beamHandle->size(); ++i) {

	    const art::Ptr<sim::BeamGateInfo> beam_ptr(beamHandle,i);
	    
	    if( (beam_ptr->BeamType() == ::sim::kBNB  && ch == kFEMChannelBNB) ||
		(beam_ptr->BeamType() == ::sim::kNuMI && ch == kFEMChannelNuMI) )

	      fLogicGen.AddPulse(beam_ptr->Start());

	  }

	}

      }

      fLogicGen.GenWaveform((*(wfs->at(logic_index).rbegin())));
      
    }
    //size_t bnb_index  = geom->NOpChannels();
    //size_t numi_index = bnb_index+1;

    // BNB
    

    // NuMI

    // Store
    if(wfs->size())
      evt.put(std::move(wfs));

    // Make sure to free memory
    fOpticalGen.Reset();

  }
} 
/** @} */ // end of doxygen group 

