#ifndef WFALGODIGITIZEDSPE_CXX
#define WFALGODIGITIZEDSPE_CXX

#include "WFAlgoDigitizedSPE.h"

namespace opdet {
  
  //--------------------------------------------------------
  WFAlgoDigitizedSPE::WFAlgoDigitizedSPE() : WFAlgoSPEBase()
  //--------------------------------------------------------
  {
    Reset();
    fSPE.clear();
    //fSPETime = util::TimeService::GetME().OpticalClock();
    art::ServiceHandle<util::TimeService> ts;
    fSPETime = ts->OpticalClock();
  }

  //------------------------------
  void WFAlgoDigitizedSPE::Reset()
  //------------------------------
  {
    WFAlgoSPEBase::Reset();
  }

  //--------------------------------------------------------------
  void WFAlgoDigitizedSPE::Process(std::vector<float> &wf,
			  const ::util::ElecClock &start_time)
  //--------------------------------------------------------------
  {
    // Predefine variables to save time later
    ::util::ElecClock rel_spe_start = start_time;

    rel_spe_start.SetTime(0);

    double unit_time = fSPETime.TickPeriod();

    for(auto const &t : fPhotonTime) {

      //
      // Check if this photon should be added or not
      //

      // Time in electronics clock frame (with T0)
      //double time = ::util::TimeService::GetME().G4ToElecTime(t);
      art::ServiceHandle<util::TimeService> ts;
      double time = ts->G4ToElecTime(t);

      if(fEnableSpread)  time += fT0 + RandomServer::GetME().Gaus(fT0,fT0Sigma);
      else time += fT0;

      // If before waveform vector, ignore
      if(time < start_time.Time()) continue;

      // If after waveform vector, ignore
      if(time > (start_time.Time() + start_time.Time((int)(wf.size())))) continue;

      //
      // Add signal
      //

      // Figure out time stamp of the beginning of SPE
      rel_spe_start.SetTime(time - fSPETime.Time() - start_time.Time());

      for(size_t i=0; i < fSPE.size(); ++i ) {

	if(rel_spe_start.Ticks() >= (int)(wf.size())) break;

	if(rel_spe_start.Ticks() >= 0) {

	  if(fEnableSpread) 

	    wf.at(rel_spe_start.Ticks()) += (fGain * fSPE.at(i));

	  else

	    wf.at(rel_spe_start.Ticks()) += ( RandomServer::GetME().Gaus(fGain,fGainSigma*fGain) * fSPE.at(i) );

	}

	rel_spe_start += unit_time;
	
      }
      
    }

  }
  
  //----------------------------------------------------------------
  void WFAlgoDigitizedSPE::SetSPE( const std::vector<float> &wf,
				   const util::ElecClock &time_info)
  //----------------------------------------------------------------
  {
    if(time_info.Time() < 0 || time_info.Ticks() >= (int)(wf.size()))
      
      throw UBOpticalException(Form("Invalid WF index (%d) for the WF of size %zu",
				    time_info.Ticks(),
				    wf.size()
				    )
			       );

    fSPE.clear();
    fSPE.reserve(wf.size());
    double max_amp = 0;
    for(auto const &v : wf) {
      if(max_amp < v) max_amp = v;
      fSPE.push_back(v);
    }

    // Normalize amplitude
    for(auto &v : fSPE) v /= max_amp;

    fSPETime = time_info;
    
  }

}

#endif
