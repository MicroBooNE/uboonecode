#ifndef OPTICALDETECTORANA_OPDETWAVEFORMANA_CXX
#define OPTICALDETECTORANA_OPDETWAVEFORMANA_CXX

#include "OpDetWaveformAna.h"
#include <limits>
#include <climits>
#include <iostream>

namespace pmtana {

  OpDetWaveformAna::OpDetWaveformAna(const std::string name)
    : _name(name)
    , _hitana_tree (nullptr)
    , _wfana_tree  (nullptr)
    , _wf_tree     (nullptr)
  {
    ClearEvent();
    ClearWaveform();
    _preco_mgr.AddRecoAlgo(&_preco_alg);

    //_preco_mgr.SetPedAlgo(pmtana::kHEAD);
    //_preco_mgr.SePedSampleCosmic (  3 );
    //_preco_mgr.SetPedSampleBeam   ( 10 );
  }

  void OpDetWaveformAna::TickPeriod   ( const double period )
  { _period = period; }
  
  void OpDetWaveformAna::SetEventInfo ( const unsigned int run,
					const unsigned int subrun,
					const unsigned int event  )
  {
    _run    = run;
    _subrun = subrun;
    _event  = event;
  }
  
  void OpDetWaveformAna::AnaWaveform  ( const unsigned int ch,
					const double time_wrt_trigger,
					const std::vector<short>& wf)
  {

    if( !_hitana_tree && !_wfana_tree && !_wf_tree ) return;
    ClearWaveform();
    
    _ch = ch;
    _t_wrt_trigger = time_wrt_trigger;

    //_preco_mgr.RecoPulse(wf);

    /*
    for(size_t i = 0; i<_preco_alg.GetNPulse(); ++i) {
      auto const& p = _preco_alg.GetPulse(i);
      
      _tstart = time_wrt_trigger + p.t_start * _period;
      _tpeak  = time_wrt_trigger + p.t_max   * _period;
      _tend   = time_wrt_trigger + p.t_end   * _period;

      _q   = p.area;
      _amp = p.peak;

      _ped_mean = p.ped_mean;
      _ped_rms  = p.ped_sigma;

      if(_hitana_tree) _hitana_tree->Fill();
    }
    */
    if( _wfana_tree || _wf_tree ) {
      for(auto const& adc : wf) {
	if(adc > _max_adc) _max_adc = adc;
	if(adc < _min_adc) _min_adc = adc;
      }
      _wf_size = wf.size();


      if( _wf_tree ) _wf = wf;
    }

    if( _wfana_tree  ) _wfana_tree->Fill();
    if( _wf_tree     ) _wf_tree->Fill();

  }
  
  void OpDetWaveformAna::ClearEvent()
  {
    _run = _subrun = _event = std::numeric_limits<unsigned int>::max();
  }

  void OpDetWaveformAna::ClearWaveform()
  {
    _ch = std::numeric_limits<unsigned int>::max();
    _wf_size = 0;

    _ped_mean = _ped_rms = -1;
    _tstart = _tpeak = _tend = -1;
    _t_wrt_trigger = 0;
    _q = _amp = -1;
    _max_adc = 0;
    _min_adc = std::numeric_limits<unsigned short>::max();
    _wf_size = 0;
    _wf.clear();
  }

  void OpDetWaveformAna::AnaHit       ( TTree* ptr )
  {
    if(!ptr) {
      std::cerr << "<<" << __FUNCTION__ << ">>" << " Invalid ptr!" << std::endl;
      throw std::exception();
    }
    if(ptr->GetEntries()) {
      std::cerr << "<<" << __FUNCTION__ << ">>" << " Non-initialized TTree!" << std::endl;
      throw std::exception();
    }     
    _hitana_tree = ptr;
    _hitana_tree->Branch( "run",      &_run,      "run/i"      );
    _hitana_tree->Branch( "subrun",   &_subrun,   "subrun/i"   );
    _hitana_tree->Branch( "event",    &_event,    "event/i"    );
    _hitana_tree->Branch( "ch",       &_ch,       "ch/i"       );
    _hitana_tree->Branch( "ped_mean", &_ped_mean, "ped_mean/F" );
    _hitana_tree->Branch( "ped_rms",  &_ped_rms,  "ped_rms/F"  );
    _hitana_tree->Branch( "tstart",   &_tstart,   "tstart/D"   );
    _hitana_tree->Branch( "tpeak",    &_tpeak,    "tpeak/D"    );
    _hitana_tree->Branch( "tend",     &_tend,     "tend/D"     );
    _hitana_tree->Branch( "q",        &_q,        "q/D"        );
    _hitana_tree->Branch( "amp",      &_amp,      "amp/D"      );

  }
  void OpDetWaveformAna::AnaWaveform  ( TTree* ptr )
  {
    if(!ptr) {
      std::cerr << "<<" << __FUNCTION__ << ">>" << " Invalid ptr!" << std::endl;
      throw std::exception();
    }
    if(ptr->GetEntries()) {
      std::cerr << "<<" << __FUNCTION__ << ">>" << " Non-initialized TTree!" << std::endl;
      throw std::exception();
    }
    _wfana_tree = ptr;
    _wfana_tree->Branch( "run",       &_run,      "run/i"      );
    _wfana_tree->Branch( "subrun",    &_subrun,   "subrun/i"   );
    _wfana_tree->Branch( "event",     &_event,    "event/i"    );
    _wfana_tree->Branch( "ch",        &_ch,       "ch/i"       );
    _wfana_tree->Branch( "ped_mean",  &_ped_mean, "ped_mean/F" );
    _wfana_tree->Branch( "ped_rms",   &_ped_rms,  "ped_rms/F"  );
    _wfana_tree->Branch( "max_adc",   &_max_adc,  "max_adc/s"  );
    _wfana_tree->Branch( "min_adc",   &_min_adc,  "min_adc/s"  );
    _wfana_tree->Branch( "wf_size",   &_wf_size,  "wf_size/i"  );

  }
  void OpDetWaveformAna::SaveWaveform ( TTree* ptr )
  {
    if(!ptr) {
      std::cerr << "<<" << __FUNCTION__ << ">>" << " Invalid ptr!" << std::endl;
      throw std::exception();
    }
    if(ptr->GetEntries()) {
      std::cerr << "<<" << __FUNCTION__ << ">>" << " Non-initialized TTree!" << std::endl;
      throw std::exception();
    }
    _wf_tree = ptr;
    _wf_tree->Branch( "run",      &_run,      "run/i"      );
    _wf_tree->Branch( "subrun",   &_subrun,   "subrun/i"   );
    _wf_tree->Branch( "event",    &_event,    "event/i"    );
    _wf_tree->Branch( "ch",       &_ch,       "ch/i"       );
    _wf_tree->Branch( "t_wrt_trigger", &_t_wrt_trigger, "t_wrt_trigger/D"   );
    _wf_tree->Branch( "ped_mean", &_ped_mean, "ped_mean/F" );
    _wf_tree->Branch( "ped_rms",  &_ped_rms,  "ped_rms/F"  );
    _wf_tree->Branch( "max_adc",  &_max_adc,  "max_adc/s"  );
    _wf_tree->Branch( "min_adc",  &_min_adc,  "min_adc/s"  );
    _wf_tree->Branch( "wf", "std::vector<short>", &_wf  );
  }

}

#endif
