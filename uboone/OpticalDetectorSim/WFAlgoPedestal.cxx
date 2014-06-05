#ifndef WFALGOPEDESTAL_CXX
#define WFALGOPEDESTAL_CXX

#include "WFAlgoPedestal.h"

namespace opdet {

  //---------------------------------------------
  WFAlgoPedestal::WFAlgoPedestal() : WFAlgoBase()
  //---------------------------------------------
  {
    fPedMean  = 2048;
    fPedSigma = 0;
  }

  //---------------------------------------------------------
  void WFAlgoPedestal::SetPedestal(double mean, double sigma)
  //---------------------------------------------------------
  {
    fPedMean  = mean;
    fPedSigma = sigma;
  }


  //-------------------------------------------------------------------
  void WFAlgoPedestal::Process(std::vector<float> &wf,
			       const ::util::ElecClock &start_time)
  //-------------------------------------------------------------------
  { 
    if(fPedSigma == 0)
      
      for(auto &v : wf) v += fPedMean;

    else

      for(auto &v : wf) v += RandomServer::GetME().Gaus(fPedMean,fPedSigma);

  }   
    

}

#endif


