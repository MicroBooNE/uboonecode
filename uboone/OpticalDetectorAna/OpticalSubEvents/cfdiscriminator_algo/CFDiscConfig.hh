#ifndef __CFDiscConfig__
#define __CFDiscConfig__

namespace cpysubevent {

  class CFDiscConfig {

  public:
    CFDiscConfig();
    ~CFDiscConfig();

    int delay; // number of previous ticks from which to calculate difference
    int threshold; // ADC difference value that triggers discriminator
    int deadtime; // deadtime in ticks before next allowed discriminator firing
    int width; // time in ticks after fire to fill discriminator fires
    int gate;  // coincidence gate
  };
}
#endif
