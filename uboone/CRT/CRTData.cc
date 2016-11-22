#include "CRT/CRTDatas.hh"

namespace crt{

  CRTData::CRTData(): fChannel(0), fT0(0), fT1(0){
  }
  CRTData::CRTData(uint32_t channel, uint32_t t0, uint32_t t1):
    fChannel(channel),
    fT0(t0),
    fT1(t1) {
    }
  CRTData::~CRTData(){
  }
  uint32_t CRTData::Channel(){
    return self->fChannel;
  }
  uint32_t CRTData::T0(){
    return self->fT0;
  }
  uint32_t CRTData::T1(){
    return self->fT1;
  }

}
