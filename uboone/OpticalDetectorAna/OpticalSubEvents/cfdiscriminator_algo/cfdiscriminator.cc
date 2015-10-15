
#include "cfdiscriminator.hh"
#include <iostream>

namespace cpysubevent {
  void runCFdiscriminatorCPP( std::vector< int >& t_fire, std::vector< int >& amp_fire, std::vector< int >& maxt_fire, std::vector< int >& diff_fire,
			      double* waveform, int delay, int threshold, int deadtime, int width, int arrlen ) {
    
    // fill diff vector
    //std::cout << "Waveform: ";
    std::vector<float> diff( arrlen, 0.0);
    for ( int tdc=delay; tdc<arrlen-delay; tdc++ ) {
      //std::cout << waveform[tdc] << ", ";
      diff.at( tdc ) = waveform[tdc]-waveform[tdc-delay];
    }
    //std::cout << std::endl;

    // reset vectors
    t_fire.clear();
    t_fire.reserve(20);
    diff_fire.clear();
    diff_fire.reserve(20);
    amp_fire.clear();
    amp_fire.reserve(20);
    maxt_fire.clear();
    maxt_fire.reserve(20);

    // determine time
    //std::cout << "Diff: ";
    for ( int t=0; t<arrlen; t++ ) {
      //std::cout << diff.at(t) << ", ";
      if ( diff.at(t)>threshold && 
	   ( t_fire.size()==0 || ( t_fire.at( t_fire.size()-1 )+deadtime<t ) ) ) {
	t_fire.push_back( t-delay );
	//diff_fire.push_back( int(diff.at(t)) );
      }
    }
    //std::cout << std::endl;

    // determine max amp
    for ( std::vector< int >::iterator it=t_fire.begin(); it!=t_fire.end(); it++ ) {
      int trig = *it;
      int end = trig+width;
      if ( end>arrlen )
	end = arrlen;
      int maxamp = (int)waveform[trig];
      int maxt = trig;
      double maxdiff = diff.at(trig);
      for ( int t=trig; t<end; t++ ) {
	if ( maxamp < (int)waveform[t] ) {
	  maxamp = (int)waveform[t];
	  maxt = t;
	}
	if ( maxdiff < diff.at(t) )
	  maxdiff = diff.at(t);
      }
      
      amp_fire.push_back( maxamp );
      maxt_fire.push_back( maxt );
      diff_fire.push_back( maxdiff );
    }
  }

  void runCFdiscriminatorCPP( std::vector< int >& t_fire, std::vector< int >& amp_fire, std::vector< int >& maxt_fire, std::vector< int >& diff_fire,
			      double* waveform, int delay, int threshold, int deadtime, int width, int arrlen );

}

