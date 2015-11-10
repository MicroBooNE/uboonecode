#ifndef WFALGOANALYTICALSPE_CXX
#define WFALGOANALYTICALSPE_CXX

#include "WFAlgoAnalyticalSPE.h"

#include "Utilities/DetectorClocksService.h"

namespace opdet {
  
  //----------------------------------------------------------
  WFAlgoAnalyticalSPE::WFAlgoAnalyticalSPE() : WFAlgoSPEBase()
  //----------------------------------------------------------
  {
    Reset();
  }

  //-------------------------------
  void WFAlgoAnalyticalSPE::Reset()
  //-------------------------------
  {
    WFAlgoSPEBase::Reset();
  }

  //--------------------------------------------------------------------
  void WFAlgoAnalyticalSPE::Process(std::vector<float> &wf,
				    const ::util::ElecClock &start_time)
  //--------------------------------------------------------------------
  {
    // Predefine variables to save time later
    ::util::ElecClock rel_spe_start = start_time;

    auto const* ts = lar::providerFrom<util::DetectorClocksService>();

    rel_spe_start.SetTime(0);

    for(auto const &t : fPhotonTime) {

      //
      // Check if this photon should be added or not
      //

      // Time in electronics clock frame (with T0)
      //double time = ::util::TimeService::GetME().G4ToElecTime(t);
      double time = ts->G4ToElecTime(t);

      if(fEnableSpread) time += RandomServer::GetME().Gaus(fT0,fT0Sigma) * 1.e-3;
      else time += fT0 * 1.e-3;

      // If before waveform vector, ignore
      if(time < start_time.Time()) continue;

      // If after waveform vector, ignore
      if(time > (start_time.Time() + start_time.Time((int)(wf.size())))) continue;
      
      // Figure out time stamp of the beginning of SPE
      rel_spe_start.SetTime(time - start_time.Time());

      //
      // Add signal
      //
      //bool peaked=false;
      for(size_t i=rel_spe_start.Ticks(); i<wf.size(); ++i) {

	double func_time = rel_spe_start.Time(i,0) - rel_spe_start.Time() + rel_spe_start.TickPeriod()/2.;
	
	if(func_time<0) continue;
	
	double amp = EvaluateSPE(func_time*1.e3);

	double gain = fGain;
	if(fEnableSpread) gain = RandomServer::GetME().Gaus(fGain,fGainSigma * fGain);

	amp *= gain;

	wf.at(i) += amp;
	/*
	if(!peaked && amp >0.01) peaked = true;
	else if(peaked && amp<0.01) break;
	*/
	if(func_time>0.624) break;
      }
    }
  }
  
  //------------------------------------------------------------
  double WFAlgoAnalyticalSPE::EvaluateSPE(const double x) const
  //------------------------------------------------------------
  {
    //
    // x should be in ns.
    //
    // Max @ x=62.8000 (and I believe we don't need sub pico-second accuracy) 
    //
    return (2.853e-3 * pow(x,3) * exp( -x / 20.94) - 4.988e-3 * exp( -x / 110000)) / 35.208752 / 5.9865;
    
  }

}

#endif
