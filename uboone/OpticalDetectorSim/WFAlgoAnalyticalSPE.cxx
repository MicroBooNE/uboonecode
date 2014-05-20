#ifndef WFALGOANALYTICALSPE_CXX
#define WFALGOANALYTICALSPE_CXX

#include "WFAlgoAnalyticalSPE.h"

namespace opdet {
  
  //----------------------------------------------------------
  WFAlgoAnalyticalSPE::WFAlgoAnalyticalSPE() : WFAlgoSPEBase()
  //----------------------------------------------------------
  {
    Reset();
    fSPEParam_TimeConstant = 6.; // 6 ns
    fSPEParam_PowerFactor  = 10; // "(n-1)" power factor
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

    rel_spe_start.SetTime(0);

    for(auto const &t : fPhotonTime) {

      //
      // Check if this photon should be added or not
      //

      // Time in electronics clock frame (with T0)
      double time = ::util::TimeService::GetME().G4ToElecTime(t);

      if(fEnableSpread) time += RandomServer::GetME().Gaus(fT0,fT0Sigma);
      else time += fT0;

      // If before waveform vector, ignore
      if(time < start_time.Time()) continue;

      // If after waveform vector, ignore
      if(time > (start_time.Time() + start_time.Time((int)(wf.size())))) continue;

      // Figure out time stamp of the beginning of SPE
      rel_spe_start.SetTime(time - start_time.Time());


      //
      // Add signal
      //
      bool peaked=false;
      for(size_t i=rel_spe_start.Ticks(); i<wf.size(); ++i) {

	double func_time = rel_spe_start.Time(i+1,0) - rel_spe_start.Time();

	double amp = EvaluateSPE(func_time);

	if(fEnableSpread) amp *= RandomServer::GetME().Gaus(fGain,fGainSigma * fGain);
	else amp *= fGain;

	wf.at(i) += amp;

	if(!peaked && amp >0.01) peaked = true;
	else if(amp<0.01) break;
      }
      
    }
  }
  
  //------------------------------------------------------------
  double WFAlgoAnalyticalSPE::EvaluateSPE(const double x) const
  //------------------------------------------------------------
  {
    double factorial = 1;
    for(int i=1; i<=((int)(fSPEParam_PowerFactor)); ++i)
      factorial *= (double)(i);
    return  1./(7565.e3) * pow(x,fSPEParam_PowerFactor) * exp(-x/fSPEParam_TimeConstant) / factorial;
  }

}

#endif
