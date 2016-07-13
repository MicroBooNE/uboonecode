#ifndef OVERLAY_DATAOVERLAY_RAWDIGITADDER_HARDSATURATE_CXX
#define OVERLAY_DATAOVERLAY_RAWDIGITADDER_HARDSATURATE_CXX

#include "RawDigitAdder_HardSaturate.h"
#include <limits>
#include <stdexcept>
#include <cmath>

mix::RawDigitAdder_HardSaturate::RawDigitAdder_HardSaturate(bool t):
  RawDigitAdder(t),
  _max(std::numeric_limits<short>::max()),
  _scale1(1),
  _scale2(1)
{}

void mix::RawDigitAdder_HardSaturate::SetSaturationPoint(short x)
{
  if(x<0){
    if(_throw)
      throw std::runtime_error("Error in RawDigitAdder_HardSaturate::SetSaturationPoint : point < 0");
    return;
  }
  _max = x;
}

void mix::RawDigitAdder_HardSaturate::SetScaleInput(float f, float& _scale)
{
  if(f<0){
    if(_throw)
      throw std::runtime_error("Error in RawDigitAdder_HardSaturate::SetScaleInput : scale < 0");
    return;
  }
  _scale = f;
}

void mix::RawDigitAdder_HardSaturate::AddRawDigit(short const& d1, short const& d2, short& d_out)
{
  d_out = (short)(std::round((float)d1 * _scale1)) + 
    (short)(std::round((float)d2 * _scale2));
  FixOverflow(d_out);
  if(d_out > _max)
    d_out=_max;
}

#endif
