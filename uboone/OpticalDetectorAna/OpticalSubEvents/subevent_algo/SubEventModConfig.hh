#ifndef __SubEventModConfig__
#define __SubEventModConfig__

#include "uboone/OpticalDetectorAna/OpticalSubEvents/cfdiscriminator_algo/CFDiscConfig.hh"
//#include "CFDiscConfig.hh"

namespace subevent {


  class SubEventModConfig {

  public:
    SubEventModConfig();
    ~SubEventModConfig();

    double spe_sigma;
    double fastfraction;
    double slowfraction;
    double fastconst_ns;
    double slowconst_ns;
    double noslowthreshold;

    int npresamples;    
    int pedsamples;
    double pedmaxvar;
    double nspersample;
    int maxchflashes;

    int hgslot;
    int lgslot;

    int flashgate;
    int maxsubeventloops;
    double ampthresh;
    int hitthresh;

    cpysubevent::CFDiscConfig cfdconfig;
    cpysubevent::CFDiscConfig cfdconfig_pass2;

  };

}

#endif
