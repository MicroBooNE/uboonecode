/**
 * \file RawDigitAdder_HardSaturate.h
 *
 * \ingroup DataOverlay
 * 
 * \brief Defintion for a class to add two vectors together,
 *        and give an "added" waveform. Takes in a saturation point.
 *
 *
 * @author wketchum
 */

/** \addtogroup DataOverlay

    @{*/
#ifndef OVERLAY_DATAOVERLAY_RAWDIGITADDER_HARDSATURATE_H
#define OVERLAY_DATAOVERLAY_RAWDIGITADDER_HARDSATURATE_H

#include <vector>
#include <string>
#include "RawDigitAdder.h"

/**
   \class RawDigitAdder_HardSaturate
   Add two vectors together. Needs a saturation point set, where
   everything above that point is just the max.
   Allows for a scale factor to be applied to inputs: defaults to 1.
   
*/
namespace mix {
  class RawDigitAdder_HardSaturate;
}

class mix::RawDigitAdder_HardSaturate : public mix::RawDigitAdder {

public:

  RawDigitAdder_HardSaturate(bool t=true);
  void SetSaturationPoint(short x);
  
  void SetScaleFirstInput(float f1)  { SetScaleInput(f1,_scale1); }
  void SetScaleSecondInput(float f2) { SetScaleInput(f2,_scale2); }

  void SetScaleInputs(float f1, float f2)
  { SetScaleFirstInput(f1); SetScaleSecondInput(f2); }
  
  std::string Name() { return "RawDigitAdder_HardSaturate"; }
  
 private:

  short _max;
  float _scale1,_scale2;
  void SetScaleInput(float f, float& _scale);
  void AddRawDigit( short const&, short const&, short&);

};

#endif
/** @} */ // end of doxygen group 

