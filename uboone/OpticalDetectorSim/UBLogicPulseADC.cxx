#ifndef UBLogicPulseADC_CXX
#define UBLogicPulseADC_CXX

#include "UBLogicPulseADC.h"

namespace opdet {

  //----------------------------------------
  UBLogicPulseADC::UBLogicPulseADC() : UBADCBase()
  //----------------------------------------
  {
    Reset();
  }

  //---------------------------
  void UBLogicPulseADC::Reset()
  //---------------------------
  {
    UBADCBase::Reset();
    fPulseTime.clear();
    fSPE.SetT0(0,0);
    fSPE.EnableSpread(false);
  }

  //----------------------------------------------------------------
  void UBLogicPulseADC::SetPulses(const std::vector<double>& g4time)
  //----------------------------------------------------------------
  {
    fPulseTime.clear();
    fPulseTime.reserve(g4time.size());
    for(auto const &v : g4time) fPulseTime.push_back(v);
  }

  //----------------------------------------------------------------------
  void UBLogicPulseADC::GenWaveform(std::vector<unsigned short>& logic_wf)
  //----------------------------------------------------------------------
  {
    //
    // Initialize
    //
    size_t nticks = fTimeInfo.Ticks(fDuration);
    std::vector<float> logic_tmp_wf(nticks,0);

    //
    // Generate Signal & DarkNoise
    //
    // Configure to generate high gain SPE
    fSPE.SetPhotons(fPulseTime);

    fSPE.Process(logic_tmp_wf,fTimeInfo);

    //
    // Simulate pedestal
    //
    fPED.Process(logic_tmp_wf,fTimeInfo);

    // Make sure algorithms did not alter the waveform size

    if(logic_tmp_wf.size()!=nticks)

      throw UBOpticalException("Waveform length changed (prohibited)!");

    //
    // Digitize amplitude
    //
    logic_wf.clear();
    Digitize(logic_tmp_wf,logic_wf);
    
  }
  
}

#endif
