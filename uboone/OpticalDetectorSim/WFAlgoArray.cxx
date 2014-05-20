#ifndef WFALGOARRAY_CXX
#define WFALGOARRAY_CXX

#include "WFAlgoArray.h"

namespace opdet {

  //---------------------------------------
  WFAlgoArray::WFAlgoArray() : WFAlgoBase()
  //---------------------------------------
  {
    fAlgoArray.clear();
  }

  //------------------------------------------------------------
  void WFAlgoArray::Process(std::vector<float> &wf,
			    const ::util::ElecClock &start_time)
  //------------------------------------------------------------
  {
    fTimeArray.clear();
    fTimeArray.resize(fAlgoArray.size(),0);

    for(size_t i=0; i<fAlgoArray.size(); ++i) {

      fWatch.Start();
      fAlgoArray.at(i)->Process(wf,start_time);
      fTimeArray.at(i) = fWatch.RealTime();

    }
    
  }


}
#endif
