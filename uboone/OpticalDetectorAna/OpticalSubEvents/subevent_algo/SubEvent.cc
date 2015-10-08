#include "SubEvent.hh"

ClassImp( subevent::SubEvent );

namespace subevent {

  SubEvent::SubEvent() {
    tstart_sample = -1;
    tend_sample = -1;
    totpe = 0.0;
    sumflash30 = sumfcomp_gausintegral = 0.0;
    maxamp = 0.0;
  }
  
  SubEvent::~SubEvent() {}

}
