#ifndef WFALGOBASE_CXX
#define WFALGOBASE_CXX

#include "WFAlgoBase.h"

namespace opdet {

  void WFAlgoBase::CombineVector(std::vector<float> &out_wf,
				 const size_t start_index,
				 const std::vector<float> &in_wf) const
  {
    
    for(size_t i=start_index; i<out_wf.size(); ++i)

      out_wf.at(i) = in_wf.at(i-start_index);

  }

}

#endif
