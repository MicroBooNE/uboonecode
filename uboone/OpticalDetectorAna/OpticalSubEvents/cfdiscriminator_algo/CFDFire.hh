#ifndef __CFDFire__

#ifdef __MAKE_ROOT_DICT__
#include "TObject.h"
#endif

namespace cfd {

#ifdef __MAKE_ROOT_DICT__
  class CFDFire : public TObject {
#else
  class CFDFire {
#endif

  public:

    CFDFire();
    ~CFDFire();
    
    int tfire;
    int maxamp;
    int tmax;
    int maxdiff;

#ifdef __MAKE_ROOT_DICT__
    ClassDef( CFDFire, 1 );
#endif

  };

}

#endif
