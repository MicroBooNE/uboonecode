
#include "pedestal.hh"
#include <algorithm>
#include <cmath>
#include <iostream>

namespace subevent {


  double calcPedestal( std::vector< double >& wfm, int samplelength, double variance_threshold, double default_ped ) {
    int site = 0;
    int start = 0;
    int end = 0;
    while ( start<(int)wfm.size() ) {
      start = std::max(0,site*samplelength);
      end = std::min((int)wfm.size(),(site+1)*samplelength);
      if ( (end-start)<3) {
	site++;
	continue;
      }
      double x = 0.0;
      double xx = 0.0;
      for (int i=start; i<end; i++) {
	x += wfm.at(i);
	xx += wfm.at(i)*wfm.at(i);
      }
      
      x /= double( end-start );
      xx /= double( end-start );
      double var = sqrt( xx-x*x );
      //std::cout << "site: " << site << " pedestal: " << x << " var=" <<  var << std::endl;
      if (var<variance_threshold) {
	return x;
      }
      site++;
    }// end of while loop
    return default_ped;
  }


  double removePedestal( std::vector< double >& wfm, int samplelength, double variance_threshold, double default_ped ) {
    double ped = calcPedestal( wfm, samplelength, variance_threshold, default_ped );
    for ( int i=0; i<(int)wfm.size(); i++)
      wfm.at( i ) -= ped;
    return ped;
  }
  
}
