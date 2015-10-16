#ifndef __SUBEVENT__
#define __SUBEVENT__

#ifdef __BUILD_ROOT_DICT__
#include "TObject.h"
#endif

#include "Flash.hh"
#include "FlashList.hh"

namespace subevent {

#ifdef __BUILD_ROOT_DICT__
  class SubEvent : public TObject {
#else
  class SubEvent {
#endif
    
  public:
    
    SubEvent();
    ~SubEvent();

    int tstart_sample;
    int tend_sample;
    int tmax_sample;

    double tstart_ns;
    double tend_ns;
    double tmax_ns;
    
    double maxamp;
    double totpe;
    double pe30;
    double totpe_1; // first pass
    double pe30_1;  // first pass
    double sumflash30;
    double sumfcomp_gausintegral;
    
    FlashList flashes;        // first pass flashes
    FlashList flashes_pass2;  // second pass flashes
    
#ifdef __BUILD_ROOT_DICT__
    ClassDef( SubEvent, 1 );
#endif
    
  };

}

#endif
