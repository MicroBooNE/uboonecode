#ifndef __SUBEVENT__
#define __SUBEVENT__

#include "TObject.h"
#include "Flash.hh"
#include "FlashList.hh"

namespace subevent {

  class SubEvent : public TObject {
    
  public:
    
    SubEvent();
    ~SubEvent();

    int tstart_sample;
    int tend_sample;
    double tstart_ns;
    double tend_ns;

    double maxamp;
    double totpe;
    double sumflash30;
    double sumfcomp_gausintegral;
    
    FlashList flashes;
    
    ClassDef( SubEvent, 1 );
    
  };

}

#endif
