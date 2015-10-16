#ifndef __FLASH_HH__
#define __FLASH_HH__

#ifdef __BUILD_ROOT_DICT__
#include "TObject.h"
#endif
#include <vector>

namespace subevent {
#ifdef __BUILD_ROOT_DICT__
  class Flash : public TObject {
#else
  class Flash {
#endif

  public:
    
    Flash();
    Flash( int ch, int tstart, int tend, int tmax, float maxamp, std::vector< double >& expectation, std::vector< double >& waveform );
    Flash( const Flash& orig ); // copy constructor
    ~Flash();
    
    void storeWaveform( const std::vector< double >& waveform );
    void storeExpectation( const std::vector< double >& expectation );

    int ch;
    int tstart;
    int tend;
    int tmax;
    double maxamp;
    double area;
    double area30;
    double fcomp_gausintegral;
    bool claimed;
    std::vector< double > expectation;
    std::vector< double > waveform;
#ifdef __BUILD_ROOT_DICT__
    ClassDef( Flash, 1 );
#endif
  };

}

#endif
