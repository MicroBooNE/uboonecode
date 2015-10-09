#include <iostream>
#include <math.h>

#include "scintresponse.hh"

namespace subevent {

  void calcScintResponseCPP( std::vector< double >& fexpectation, int tstart, int tend, int maxt, float sig, float maxamp, float fastconst, float slowconst, float nspertick ) {
    
    //slow component shape: expo convolved with gaus
    //float t_smax = 95.0;    // peak of only slow component. numerically solved for det. smearing=3.5*15.625 ns, decay time const= 1500 ns
    float t_fmax = 105.0;   // numerically solved for det. smearing=3.5*15.625 ns, decay time const= 6 ns
    float smax = exp( sig*sig/(2*slowconst*slowconst) - t_fmax/slowconst )*(1 - erf( (sig*sig - slowconst*t_fmax )/(sqrt(2.0)*sig*slowconst ) ) );
    // normalize max at fast component peak
    float As = 0.3*maxamp/smax; 
    //std::cout << "smax:" << smax << " As=" << As << std::endl;
  
    // fast component: since time const is smaller than spe response, we model as simple gaussian
    float Af = 0.8*maxamp;

    // 0.3 and 0.8 are kind of made up!
//     float fast_integral = Af*sig*sqrt(2.0*3.14159);
//     float target_slow_cdf_fract = 1.0/fast_integral;
//     float target_slow_t = -log( 1.0-target_slow_cdf_fract )*slowconst;
//     std::cout << "fast=" << fast_integral << " maxamp=" << maxamp << " t=" << target_slow_t << std::endl;

    if ( maxamp<30.0 ) { // single pe level
      Af = maxamp;
      As = 0.0;
    }

    int arrlen = tend-tstart;
    bool rising = true;
    float t = 0.0;
    float tmax_ns = maxt*nspertick;
    float tstart_ns = tstart*nspertick;

    //texpectation.clear();
    //texpectation.reserve( arrlen );
    fexpectation.clear();
    fexpectation.reserve( arrlen );

    for ( int tdc=0; tdc<arrlen; tdc++ ) {
      // convert to time
      t = (tstart_ns + float(tdc*nspertick)) - tmax_ns;
      float farg = (t)/sig;
      float f = Af*exp( -0.5*farg*farg );
      float s = As*exp( sig*sig/(2*slowconst*slowconst) - (t)/slowconst )*(1 - erf( (sig*sig - slowconst*(t) )/(sqrt(2.0)*sig*slowconst ) ) );
      float amp = f+s;
      // std::cout << tdc << " " << t << " " << f << " " << s << std::endl;
      fexpectation.push_back( amp ); // amp vs tdc
      if ( rising && amp>5.0 )
	rising = false;
      //else if ( (!rising && amp<0.1) || (t>target_slow_t && amp<5.0) )
      else if ( (!rising && amp<0.1) )
	break;
    }
  }

}
