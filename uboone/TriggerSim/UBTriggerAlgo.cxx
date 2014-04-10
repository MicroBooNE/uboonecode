////////////////////////////////////////////////////////////////////////
//
//  UBTriggerAlgo source
//
////////////////////////////////////////////////////////////////////////

#ifndef UBTRIGGERALGO_CXX
#define UBTRIGGERALGO_CXX

#include <iostream>
#include <sstream>
#include <TString.h>
#include "Utilities/TimeService.h"
#include "UBTriggerTypes.h"
#include "UBTriggerAlgo.h"

namespace trigger{

  /*
  //##################################################################
  UBTriggerAlgo::UBTriggerAlgo(fhicl::ParameterSet const& pset)
    : _tpc_clock(pset.get< fhicl::ParameterSet >("ClockTPC")),
      _pmt_clock(pset.get< fhicl::ParameterSet >("ClockPMT")),
      _trigger_clock(pset.get< fhicl::ParameterSet >("ClockTrigger"))
  //##################################################################
  {    
    _trigger_counter=0;

    // Debug flag
    SetDebugMode(pset.get<bool>("DebugMode"));

    // Deadtime 
    SetDeadtime(pset.get<uint16_t>("Deadtime") );

    // Mask
    SetMask(pset.get<std::vector<uint16_t> >("Mask"));

    // Prescale
    SetPrescale(pset.get<std::vector<bool> >("Prescale"));

    // BNB gate parameters
    SetBNBParams( pset.get<uint16_t>("GateWidthBNB"),
                  pset.get<uint16_t>("TrigDelayBNB"),
		  pset.get<uint16_t>("CosmicStartBNB"),
		  pset.get<uint16_t>("CosmicEndBNB") );

    // NuMI gate parameters
    SetNuMIParams( pset.get<uint16_t>("GateWidthNuMI"),
                   pset.get<uint16_t>("TrigDelayNuMI"),
		   pset.get<uint16_t>("CosmicStartNuMI"),
		   pset.get<uint16_t>("CosmicEndNuMI") );


    ClearPredefinedTriggers();

    // BNB trigger
    _bnb_timings = pset.get<std::vector<double> >("PermTriggerBNB",std::vector<double>())

    // NuMI trigger
    _numi_timings = pset.get<std::vector<double> >("PermTriggerNuMI",std::vector<double>())

    // Calib trigger
    _calib_timings = pset.get<std::vector<double> >("PermTriggerCalib",std::vector<double>())

    // Calib trigger
    _calib_timings = pset.get<std::vector<double> >("PermTriggerCalib",std::vector<double>())

    // Ext trigger
    _ext_timings = pset.get<std::vector<double> >("PermTriggerExt",std::vector<double>())

    // PC trigger
    _pc_timings = pset.get<std::vector<double> >("PermTriggerPC",std::vector<double>())

    ReportConfig();

    ClearInputTriggers();
  }
  */

  //##############################################################
  UBTriggerAlgo::UBTriggerAlgo() : _mask(9,0),
				   _prescale(9,false)
  //##############################################################
  {
    art::ServiceHandle<util::TimeService> ts;
    _trig_clock = ts->TriggerClock();
    _pmt_clock  = ts->OpticalClock();
    _tpc_clock  = ts->TPCClock();
						
    // LArLight version constructor
    SetDeadtime(4);
    SetBNBParams(102,256,0,0);
    SetNuMIParams(102,256,0,0);
  }

  //######################################################
  void UBTriggerAlgo::Report(const std::string &msg) const
  //######################################################
  {
    std::cout << msg.c_str() << std::endl;
  }

  //##############################################################
  void UBTriggerAlgo::RaiseTriggerException(std::string msg) const
  //##############################################################
  {
    //throw cet::exception(__FUNCTION__) << msg.c_str();
    throw std::invalid_argument(msg.c_str());
    std::exit(1);
  }
  
  //#########################################################
  void UBTriggerAlgo::SetBNBParams(unsigned short width,
				   unsigned short delay,
				   unsigned short cosmic_min,
				   unsigned short cosmic_max) 
  //#########################################################
  {
    if(_bnb_cosmic_allow_max < _bnb_cosmic_allow_min) {
      RaiseTriggerException("BNB Cosmic allow window max is smaller than min!");
      _bnb_cosmic_allow_max = _bnb_cosmic_allow_min;
    }
    _bnb_gate_width       = width;
    _bnb_delay            = delay;
    _bnb_cosmic_allow_min = cosmic_min;
    _bnb_cosmic_allow_max = cosmic_max;
  }

  //######################################################
  void UBTriggerAlgo::SetNuMIParams(unsigned short width,
			       unsigned short delay,
			       unsigned short cosmic_min,
			       unsigned short cosmic_max) 
  //######################################################
  {
    if(_numi_cosmic_allow_max < _numi_cosmic_allow_min) {

      _numi_cosmic_allow_max = _numi_cosmic_allow_min;
    }    
    _numi_gate_width       = width;
    _numi_delay            = delay;
    _numi_cosmic_allow_min = cosmic_min;
    _numi_cosmic_allow_max = cosmic_max;
  }


  //#############################################################
  void UBTriggerAlgo::SetMask(unsigned char index, uint32_t mask)
  //#############################################################
  {
    if(index>9) RaiseTriggerException("Index >9 cannot be set!");
    _mask[index]=mask;
  }

  //#################################################################
  void UBTriggerAlgo::SetPrescale(unsigned char index, bool prescale)
  //#################################################################
  {
    if(index>9) RaiseTriggerException("Index >9 cannot be set!");
    _prescale[index]=prescale;
  }

  //##################################################################
  void UBTriggerAlgo::SetMask(const std::vector<unsigned short> &mask) 
  //##################################################################
  {
    if(_mask.size()!=9) {
      RaiseTriggerException(Form("Length of masks (=%zu) is invalid! Initializing to right length...",
				 mask.size()));
      _mask.resize(9,0);
    }else
      _mask = mask;

  }

  //##################################################################
  void UBTriggerAlgo::SetPrescale(const std::vector<bool> &prescale) 
  //##################################################################
  {
    if(_prescale.size()!=9) {
      RaiseTriggerException(Form("Length of prescales (=%zu) is invalid! Initializing to right length...",
				 prescale.size()));
      _prescale.resize(9,0);
    }else
      _prescale = prescale;

  }

  //######################################
  void UBTriggerAlgo::ClearInputTriggers()
  //######################################
  {
    _candidates.clear();

    for( auto const t : _bnb_timings ) AddTriggerBNB(t);

    for( auto const t : _numi_timings ) AddTriggerNuMI(t);

    for( auto const t : _calib_timings ) AddTriggerCalib(t);

    for( auto const t : _ext_timings ) AddTriggerExt(t);

    for( auto const t : _pc_timings ) AddTriggerPC(t);

  }

  //############################################################################
  util::ElecClock UBTriggerAlgo::BNBStartTime(const util::ElecClock& time) const
  //############################################################################
  {
    art::ServiceHandle<util::TimeService> ts;
    auto clock = ts->OpticalClock(_pmt_clock.Time(time.Time()));
    //auto clock = util::TimeService::GetME().OpticalClock(_pmt_clock.Time(time.Time()));

    clock += _bnb_delay;

    return clock;
  }

  //#############################################################################
  util::ElecClock UBTriggerAlgo::NuMIStartTime(const util::ElecClock& time) const
  //#############################################################################
  {
    art::ServiceHandle<util::TimeService> ts;
    auto clock = ts->OpticalClock(_pmt_clock.Time(time.Time()));
    //auto clock = util::TimeService::GetME().OpticalClock(_pmt_clock.Time(time.Time()));

    clock += _numi_delay;

    return clock;
  }

  //###################################################################################
  util::ElecClock UBTriggerAlgo::ReadOutStartTimeTPC(const util::ElecClock& time) const
  //###################################################################################
  {
    art::ServiceHandle<util::TimeService> ts;
    auto clock = ts->OpticalClock(_pmt_clock.Time(time.Time()));
    //auto clock = util::TimeService::GetME().OpticalClock(_pmt_clock.Time(time.Time()));

    clock += _tpc_readout_offset;

    return clock;

  }

  //#######################################################################################
  util::ElecClock UBTriggerAlgo::ReadOutStartTimeOptical(const util::ElecClock& time) const
  //#######################################################################################
  {
    art::ServiceHandle<util::TimeService> ts;
    auto clock = ts->OpticalClock(_pmt_clock.Time(time.Time()));
    //auto clock = util::TimeService::GetME().OpticalClock(_pmt_clock.Time(time.Time()));

    clock += _readout_frame_offset * clock.FrameTicks();

    return clock;
  }

  //##############################################################
  void UBTriggerAlgo::AddTriggerCalib(const util::ElecClock &time)
  //##############################################################
  {
    _pmt_clock.SetTime(_pmt_clock.Time(time.Time()));
    // Calibration triggers
    AddTrigger(raw::Trigger(0,
			    _pmt_clock.Time(),
			    _pmt_clock.Time(),
			    ReadOutStartTimeTPC(_pmt_clock).Time(),
			    ReadOutStartTimeOptical(_pmt_clock).Time(),
			    ( (0x1) << trigger::kTriggerCalib )) 
	       );
    _pmt_clock.SetTime(0);
  }

  //############################################################
  void UBTriggerAlgo::AddTriggerExt(const util::ElecClock& time)
  //############################################################
  {
    _pmt_clock.SetTime(_pmt_clock.Time(time.Time()));
    // EXT triggers
    AddTrigger(raw::Trigger(0,
			    _pmt_clock.Time(),
			    _pmt_clock.Time(),
			    ReadOutStartTimeTPC(_pmt_clock).Time(),
			    ReadOutStartTimeOptical(_pmt_clock).Time(),
			    ( (0x1) << trigger::kTriggerEXT ))
	       );
    _pmt_clock.SetTime(0);
  }

  //###########################################################
  void UBTriggerAlgo::AddTriggerPC(const util::ElecClock& time)
  //###########################################################
  {
    _pmt_clock.SetTime(_pmt_clock.Time(time.Time()));
    // PC triggers
    AddTrigger(raw::Trigger(0,
			    _pmt_clock.Time(),
			    _pmt_clock.Time(),
			    ReadOutStartTimeTPC(_pmt_clock).Time(),
			    ReadOutStartTimeOptical(_pmt_clock).Time(),
			    ( (0x1) << trigger::kTriggerPC ))
	       );
    _pmt_clock.SetTime(0);
  }

  //#######################################
  void UBTriggerAlgo::ReportConfig() const
  //#######################################
  {
    std::ostringstream msg;

    msg
      << std::endl
      << " UBTriggerAlgo Configuration:              " << std::endl
      << "---------------------------------------------" << std::endl;

    msg << " Debug Mode ... " << (_debug_mode ? "enabled!" : "disabled!") << std::endl
	<< std::endl;

    for(size_t i=0; i<_mask.size(); ++i)

      msg << Form("  Mask %zu     : %d", i, _mask[i]) << std::endl;
    
    msg << std::endl;

    for(size_t i=0; i<_prescale.size(); ++i)

      msg << Form("  Prescale %zu : %d", i, _prescale[i]) << std::endl;

    msg << std::endl;

    msg
      << Form(" Trigger Deadtime    : %d frames",_deadtime) << std::endl
      << std::endl
      << Form("  NuMI Trigger Delay : %d samples",_numi_delay) << std::endl
      << Form("  NuMI Cosmic Start  : %d samples",_numi_cosmic_allow_min) << std::endl
      << Form("  NuMI Cosmic End    : %d samples",_numi_cosmic_allow_max) << std::endl
      << std::endl
      << Form("  BNB Trigger Delay  : %d samples",_bnb_delay) << std::endl
      << Form("  BNB Cosmic Start   : %d samples",_bnb_cosmic_allow_min) << std::endl
      << Form("  BNB Cosmic End     : %d samples",_bnb_cosmic_allow_max) << std::endl
      << "---------------------------------------------" << std::endl
      << std::endl;

    //mf::LogInfo(__FUNCTION__) << msg.str();
    Report(msg.str());
  }

  //#######################################################################################
  const raw::Trigger UBTriggerAlgo::CombineTriggers(const raw::Trigger &trigger1, 
						    const raw::Trigger &trigger2)
  //#######################################################################################
  {

    auto trig1_time = trigger1.TriggerTime();
    auto trig2_time = trigger2.TriggerTime();

    auto trig1_number = trigger1.TriggerNumber();
    auto trig2_number = trigger2.TriggerNumber();

    if(_debug_mode)
      
      Report(Form("  Attempting to combine two triggers: (%d,%d) and (%d,%d)",
		  _trig_clock.Sample(trig1_time),
		  _trig_clock.Frame(trig1_time),
		  _trig_clock.Sample(trig2_time),
		  _trig_clock.Frame(trig2_time))
	     );

    if( trig1_number != trig2_number ) {
      RaiseTriggerException("Cannot combine triggers with different trigger counters!");
      return raw::Trigger();
    }

    if( _trig_clock.Frame(trig1_time) != _trig_clock.Frame(trig2_time) ) {
      RaiseTriggerException("Cannot combine triggers in different frames!");
      return raw::Trigger();
    }

    if( _trig_clock.Sample(trig1_time) != _trig_clock.Sample(trig2_time) ) {
      RaiseTriggerException("Cannot combine triggers in different trigger clock sample number!");
      return raw::Trigger();
    }

    auto const res_number = trig1_number;
    auto const res_frame  = _pmt_clock.Frame(trig1_time);
    auto const trig1_sample = _pmt_clock.Sample(trig1_time);
    auto const trig2_sample = _pmt_clock.Sample(trig2_time);
    auto const res_sample = (trig1_sample < trig2_sample ? trig1_sample : trig2_sample);

    // Construct beam gate region. Note that this is relevant only if either of
    // two candidates is beam. Else do nothing.
    unsigned int beam_frame  = 0;
    unsigned int beam_sample = 0;

    // Case1: both are beam ... not supported for now
    if( (trigger1.Triggered(trigger::kTriggerBNB) || trigger1.Triggered(trigger::kTriggerNuMI)) &&
	(trigger2.Triggered(trigger::kTriggerBNB) || trigger2.Triggered(trigger::kTriggerNuMI)) ) {
      RaiseTriggerException("Combining two beam gates not supported for now!");
      return raw::Trigger();
    }
    // Case2: only trigger 1 is beam
    else if( (trigger1.Triggered(trigger::kTriggerBNB) || trigger1.Triggered(trigger::kTriggerNuMI)) ) {
      
      beam_frame  = _pmt_clock.Frame(trigger1.BeamGateTime());
      beam_sample = _pmt_clock.Sample(trigger1.BeamGateTime());

    }
    // Case3: only trigger 2 is beam
    else if( (trigger2.Triggered(trigger::kTriggerBNB) || trigger2.Triggered(trigger::kTriggerNuMI)) ) {

      beam_frame  = _pmt_clock.Frame(trigger2.BeamGateTime());
      beam_sample = _pmt_clock.Sample(trigger2.BeamGateTime());

    }
    // Case4: neither is beam ... use trigger timing
    else{

      beam_frame  = res_frame;
      beam_sample = res_sample;

    }

    if(_debug_mode) Report(Form("    Combined beam timing @ (%d,%d)", beam_sample, beam_frame));

    // Combine trigger bits
    uint32_t res_bits = 0x0;
    for(unsigned char i=0; i < 4*(sizeof(res_bits)); ++i) {
      
      if( trigger1.Triggered(i) || trigger2.Triggered(i) )
	{
	  res_bits += (uint32_t)( (0x1) << i );
	  if(_debug_mode) Report(Form("    Combined bit %d ... now %d",i,res_bits));
	}
    }
    art::ServiceHandle<util::TimeService> ts;
    auto trig_time = ts->OpticalClock( res_sample,  res_frame  );
    auto beam_time = ts->OpticalClock( beam_sample, beam_frame );
    return raw::Trigger(res_number,
			trig_time.Time(),
			beam_time.Time(),
			ReadOutStartTimeTPC(trig_time).Time(),
			ReadOutStartTimeOptical(trig_time).Time(),
			res_bits);
    
  }

  //##################################################################
  void UBTriggerAlgo::AddTrigger(const raw::Trigger &new_trigger)
  //##################################################################
  {
    auto frame  = _trig_clock.Frame(new_trigger.TriggerTime());
    auto sample = _trig_clock.Sample(new_trigger.TriggerTime());

    if(_debug_mode)

      Report(Form("Requested to add a new trigger @ (%d, %d)",sample,frame));

    auto frame_iter = _candidates.find(frame);

    if( frame_iter == _candidates.end() ) {

      _candidates[frame] = std::map<unsigned int,raw::Trigger>();
      _candidates[frame].insert(std::pair<unsigned int,raw::Trigger>(sample,new_trigger));
      
    }else{

      auto sample_iter = (*frame_iter).second.find(sample);

      if( sample_iter == (*frame_iter).second.end() )

	(*frame_iter).second.insert(std::pair<unsigned int,raw::Trigger>(sample,new_trigger));
      
      else {

	raw::Trigger combined_trigger = CombineTriggers(new_trigger, (*sample_iter).second);

	_candidates[frame][sample]=combined_trigger;

      }
    }
  }

  //##################################################################
  void UBTriggerAlgo::AddTriggerPMTCosmic(const util::ElecClock& time)
  //##################################################################
  {
    _pmt_clock.SetTime(_pmt_clock.Time(time.Time()));
    // Trigger bits
    unsigned int trig_bits=0x0;
    trig_bits += ( (0x1) << kPMTTriggerCosmic );
    trig_bits += ( (0x1) << kPMTTrigger );
    
    // Create this trigger candidate object
    raw::Trigger trig_candidate(0,
				_pmt_clock.Time(),
				_pmt_clock.Time(),
				ReadOutStartTimeTPC(_pmt_clock).Time(),
				ReadOutStartTimeOptical(_pmt_clock).Time(),
				trig_bits);
    
    // Add this trigger candidate
    AddTrigger(trig_candidate);
    _pmt_clock.SetTime(0);
  }

  //################################################################
  void UBTriggerAlgo::AddTriggerPMTBeam(const util::ElecClock& time)
  //################################################################
  {
    _pmt_clock.SetTime(_pmt_clock.Time(time.Time()));
   // Trigger bits
    unsigned int trig_bits=0x0;
    trig_bits += ( (0x1) << kPMTTriggerBeam );
    trig_bits += ( (0x1) << kPMTTrigger );
    
    // Create this trigger candidate object
    raw::Trigger trig_candidate(0,
				_pmt_clock.Time(),
				_pmt_clock.Time(),
				ReadOutStartTimeTPC(_pmt_clock).Time(),
				ReadOutStartTimeOptical(_pmt_clock).Time(),
				trig_bits);
    
    // Add this trigger candidate
    AddTrigger(trig_candidate);
    _pmt_clock.SetTime(0);
  }

  //############################################################
  void UBTriggerAlgo::AddTriggerBNB(const util::ElecClock& time)
  //############################################################
  {
    _pmt_clock.SetTime(_pmt_clock.Time(time.Time()));
    uint32_t trig_bits = ( (0x1) << trigger::kTriggerBNB );
 
    // Create this trigger candidate object
    raw::Trigger trig_candidate(0,
				_pmt_clock.Time(),
				BNBStartTime(_pmt_clock).Time(),
				ReadOutStartTimeTPC(_pmt_clock).Time(),
				ReadOutStartTimeOptical(_pmt_clock).Time(),
				trig_bits);
    
    // Add this trigger candidate
    AddTrigger(trig_candidate);
    _pmt_clock.SetTime(0);
  }

  //#############################################################
  void UBTriggerAlgo::AddTriggerNuMI(const util::ElecClock& time)
  //#############################################################
  {
    _pmt_clock.SetTime(_pmt_clock.Time(time.Time()));
    uint32_t trig_bits = ( (0x1) << trigger::kTriggerNuMI );

    // Create this trigger candidate object
    raw::Trigger trig_candidate(0,
				_pmt_clock.Time(),
				NuMIStartTime(_pmt_clock).Time(),
				ReadOutStartTimeTPC(_pmt_clock).Time(),
				ReadOutStartTimeOptical(_pmt_clock).Time(),
				trig_bits);

    // Add this trigger candidate
    AddTrigger(trig_candidate);
    _pmt_clock.SetTime(0);
  }
  
  //###############################################
  void UBTriggerAlgo::ShowCandidateTriggers() const
  //###############################################
  {
    auto frame_iter = _candidates.begin();
    std::ostringstream msg;

    msg
      << "######################################################"
      << "Trigger Inputs"
      << "######################################################" << std::endl
      << std::endl
      << Form("  %-18s : %-18s : %-18s : %-18s : %-18s : %-18s : %-32s",
	      "Trig Sample(16MHz)",
	      "Trig Frame (16MHz)",
	      "Trig Sample(64MHz)",
	      "Trig Frame (64MHz)",
	      "Beam Sample(64Mhz)",
	      "Beam Frame (64Mhz)",
	      "Bits") << std::endl
      << std::endl;

    while(frame_iter!=_candidates.end()) {
      
      auto sample_iter = (*frame_iter).second.begin();

      while(sample_iter != (*frame_iter).second.end()) {

	unsigned int decision_sample = (*sample_iter).first;
	unsigned int decision_frame  = (*frame_iter).first;

	unsigned int trigger_sample  = _pmt_clock.Sample(((*sample_iter).second.TriggerTime()));
	unsigned int trigger_frame   = _pmt_clock.Frame(((*sample_iter).second.TriggerTime()));

	unsigned int beam_sample     = _pmt_clock.Sample(((*sample_iter).second.BeamGateTime()));
	unsigned int beam_frame      = _pmt_clock.Sample(((*sample_iter).second.BeamGateTime()));
	
	msg
	  << Form("  %-18d : %-18d : %-18d : %-18d : %-18d : %-18d : ",
		  decision_sample,
		  decision_frame,
		  trigger_sample,
		  trigger_frame,
		  beam_sample,
		  beam_frame);

	for(unsigned char k=0; k<32; ++k)

	  msg << ((*sample_iter).second.Triggered(k) ? 1 : 0);

	msg << std::endl;

	sample_iter++;
      }
      frame_iter++;
    }
    msg
      << std::endl
      << "######################################################"
      << "##############"
      << "######################################################" << std::endl
      << std::endl;

    Report(msg.str());
  }

  //#########################################################################
  void UBTriggerAlgo::ProcessTrigger(std::vector<raw::Trigger> &triggers)
  //#########################################################################
  {
    triggers.clear();

    if(_debug_mode)
      
      ShowCandidateTriggers();

    // Just make sure trigger clock is only used for sample/frame getter (i.e. no time used)
    _trig_clock.SetTime(0);
    util::ElecClock trig_time = _trig_clock;

    time_window_t deadtime    (_trig_clock,_trig_clock);
    time_window_t bnb_active  (_trig_clock,_trig_clock);
    time_window_t numi_active (_trig_clock,_trig_clock);
    time_window_t bnb_gate    (_trig_clock,_trig_clock);
    time_window_t numi_gate   (_trig_clock,_trig_clock);

    auto const mask0 = _mask.at(0);
    auto const mask1 = _mask.at(1);
    auto const mask8 = _mask.at(8);
    auto const scale0 = _prescale.at(0);
    auto const scale1 = _prescale.at(1);
    auto const scale8 = _prescale.at(8);

    // Loop over candidate frames
    auto frame_iter = _candidates.begin();
    while( frame_iter != _candidates.end() ) {

      // Make sure sample exists if frame is found
      auto sample_iter = (*frame_iter).second.begin();
      
      if( sample_iter == (*frame_iter).second.end() ) {
	
	RaiseTriggerException("Logic error: found candidate frame but no associated sample & trigger!");
	triggers.clear();
	return;
      }    

      // Loop over samples
      while( sample_iter != (*frame_iter).second.end() ) {

	if(_debug_mode) Report(Form("\n  Processing candidate @ %d, %d",
				    (*sample_iter).first,
				    (*frame_iter).first));

	trig_time.SetTime(_trig_clock.Time((*sample_iter).first,(*frame_iter).first));

	// If in deadtime, continue
	if( InWindow(trig_time, deadtime) ) {

	  if(_debug_mode) Report(Form("    Skipping as deadtime is (%d,%d) => (%d,%d)",
				      deadtime.first.Sample(),
				      deadtime.first.Frame(),
				      deadtime.second.Sample(),
				      deadtime.second.Frame())
				 );
	  sample_iter++;
	  continue;
	}else if(_debug_mode){
	  
	  Report(Form("    Not in deadtime: (%d,%d) => (%d,%d)",
		      deadtime.first.Sample(),
		      deadtime.first.Frame(),
		      deadtime.second.Sample(),
		      deadtime.second.Frame())
		 );
	}

	auto const pmt0   = (*sample_iter).second.Triggered(trigger::kPMTTriggerCosmic);
	auto const pmt1   = (*sample_iter).second.Triggered(trigger::kPMTTriggerBeam);
	auto const numi   = (*sample_iter).second.Triggered(trigger::kTriggerNuMI);
	auto const bnb    = (*sample_iter).second.Triggered(trigger::kTriggerBNB);
	auto const calib  = (*sample_iter).second.Triggered(trigger::kTriggerCalib);
	auto const ext    = (*sample_iter).second.Triggered(trigger::kTriggerEXT);
	auto const pc     = (*sample_iter).second.Triggered(trigger::kTriggerPC);
	auto const active = ( InWindow(trig_time, bnb_active) || InWindow(trig_time, numi_active) );

	// If BNB or NuMI, open new beam trigger gate
	if(bnb) {
	  bnb_gate.first.SetTime((*sample_iter).first, (*frame_iter).first);
	  bnb_gate.second.SetTime(bnb_gate.first.Time() + _trig_clock.Time(_bnb_gate_width));
	}
	if(numi) {
	  numi_gate.first.SetTime((*sample_iter).first, (*frame_iter).first);
	  numi_gate.second.SetTime(numi_gate.first.Time() + _trig_clock.Time(_bnb_gate_width));
	}
	
	auto const bnb_gate_active  = InWindow(trig_time, bnb_gate);
	auto const numi_gate_active = InWindow(trig_time, numi_gate);

	if(_debug_mode)
	  {
	    std::ostringstream msg;
	    msg 
	      << Form("    PMT Cosmic? %s", (pmt0  ? "yes" : "no")) << std::endl
	      << Form("    PMT Beam?   %s", (pmt1  ? "yes" : "no")) << std::endl
	      << Form("    NuMI Gate?  %s", (numi  ? "yes" : "no")) << std::endl
	      << Form("    BNB Gate?   %s", (bnb   ? "yes" : "no")) << std::endl
	      << Form("    Calib?      %s", (calib ? "yes" : "no")) << std::endl
	      << Form("    External?   %s", (ext   ? "yes" : "no")) << std::endl
	      << Form("    PC?         %s", (pc    ? "yes" : "no")) << std::endl
	      << std::endl
	      << Form("    NuMI   Window? %s", (numi_gate_active  ? "yes" : "no")) << std::endl
	      << Form("    BNB    Window? %s", (bnb_gate_active ? "yes" : "no")) << std::endl
	      << Form("    Cosmic Window? %s", (active           ? "yes" : "no")) << std::endl;
	    Report(msg.str());
	  }



	// Evaludate trigger condition p0 ... this is for PMT Cosmic trigger
	// Equation follows Nevis FPGA code including similar notations
	bool p0 = ( (numi   && pmt0 && ((0x1) & (mask0 >> 0))           ) ||
		    (bnb    && pmt0 && ((0x1) & (mask0 >> 1))           ) ||
		    (active && pmt0 && ((0x1) & (mask0 >> 2)) && !scale0) ||
		    (          pmt0 && ((0x1) & (mask0 >> 3)) && !scale0) );

	// Evaluate trigger condition p1 ... this is for PMT Beam trigger
	// Equation follows Nevis FPGA code including similar notations
	bool p1 = ( (numi   && pmt1 && ((0x1) & (mask1 >> 0))           ) || 
		    (bnb    && pmt1 && ((0x1) & (mask1 >> 1))           ) ||
		    (active && pmt1 && ((0x1) & (mask1 >> 2)) && !scale1) ||
		    (          pmt1 && ((0x1) & (mask1 >> 3)) && !scale1) );

	// Evaluate trigger condition p8 ... this is for Beam trigger input
	// Equation follows Nevis FPGA code including similar notations
	bool p8 = ( (active && ext    && ((0x1) & (mask8 >> 0)) && !scale8) ||
		    (          ext    && ((0x1) & (mask8 >> 1)) && !scale8) ||
		    (active && pc     && ((0x1) & (mask8 >> 2))           ) ||
		    (          pc     && ((0x1) & (mask8 >> 3))           ) ||
		    (bnb_gate_active  && ((0x1) & (mask8 >> 4))           ) ||
		    (numi_gate_active && ((0x1) & (mask8 >> 5))           ) ||
		    (           calib && ((0x1) & (mask8 >> 6))           ) );

	if(_debug_mode) {
	  
	  std::ostringstream msg;

	  msg << Form("    p0 condition ... %s", (p0 ? "satisfied" : "not met")) << std::endl
	      << Form("    p1 condition ... %s", (p1 ? "satisfied" : "not met")) << std::endl
	      << Form("    p8 condition ... %s", (p8 ? "satisfied" : "not met")) << std::endl;
	  Report(msg.str());
	}


	// Ask OR condition of p0, p1, and p8
	if( p0 || p1 || p8 ) {
	  // New trigger found. 

	  deadtime.first.SetTime(trig_time.Time());
	  deadtime.second.SetTime(trig_time.Time());
	  deadtime.second += (int)(_deadtime * _trig_clock.FrameTicks() - 1);
	  if(_debug_mode) {
	    std::cout<<Form("    Dead time set: (%d,%d) => (%d,%d) with trigger @ (%d,%d)",
			    deadtime.first.Sample(),
			    deadtime.first.Frame(),
			    deadtime.second.Sample(),
			    deadtime.second.Frame(),
			    trig_time.Sample(),
			    trig_time.Frame()
			    )
		     << std::endl;
	  }

	  // Assign active window for bnb/numi
	  if(bnb) {
	    bnb_active.first.SetTime(trig_time.Time());
	    bnb_active.first  += (int)(_bnb_cosmic_allow_min);
	    bnb_active.second.SetTime(trig_time.Time());
	    bnb_active.second += (int)(_bnb_cosmic_allow_max);
	  }
	  if(numi) {
	    numi_active.first.SetTime(trig_time.Time());
	    numi_active.first  += (int)(_numi_cosmic_allow_min);
	    numi_active.second.SetTime(trig_time.Time());
	    numi_active.second += (int)(_numi_cosmic_allow_max);
	  }
	  
	  // Store trigger object
	  triggers.push_back( raw::Trigger(_trigger_counter,
					   (*sample_iter).second.TriggerTime(),
					   (*sample_iter).second.BeamGateTime(),
					   (*sample_iter).second.ReadOutStartTPC(),
					   (*sample_iter).second.ReadOutStartOptical(),
					   (*sample_iter).second.TriggerBits()) 
			      );
	  _trigger_counter++;
	} // end if trigger found
	sample_iter++;
      } // end looping over samples
      frame_iter++;
    } // end looping over frames
    
  }

}

#endif
