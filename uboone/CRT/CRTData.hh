#ifndef CRTData_hh_
#define CRTData_hh_

#include<cstdint>

namespace crt {

  class CRTData {
    uint32_t fChannel;
    uint32_t fT0;
    uint32_t fT1;
   public:
    CRTData();
    CRTData(uint32_t channel, uint32_t t0, uint32_t t1);
    virtual ~CRTData();

    uint32_t Channel();
    uint32_t T0();
    uint32_t T1();
  };

}

#endif
