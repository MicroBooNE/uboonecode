#ifndef WFALGOANALYTICALSPE_CXX
#define WFALGOANALYTICALSPE_CXX

#include "WFAlgoSPEBase.h"

namespace opdet {
  
  //-------------------------------------------------------
  WFAlgoSPEBase::WFAlgoSPEBase() : WFAlgoBase()
  //-------------------------------------------------------
  {
    fEnableSpread = false;
    Reset();
  }

  
  //-------------------------
  void WFAlgoSPEBase::Reset()
  //-------------------------
  {
    fPhotonTime.clear();
    fGain      = 0; // 20 ADC / p.e.
    fGainSigma = 0; // 0 spread
    fT0        = 0; // 0 ns delay
    fT0Sigma   = 0; // 0 ns spread
  }

  //---------------------------------------------------------------
  void WFAlgoSPEBase::SetPhotons(const std::vector<double> &g4time)
  //---------------------------------------------------------------
  {
    fPhotonTime.clear();
    fPhotonTime.reserve(g4time.size());
    for(auto const &v: g4time) fPhotonTime.push_back(v);
  }

}

#endif
