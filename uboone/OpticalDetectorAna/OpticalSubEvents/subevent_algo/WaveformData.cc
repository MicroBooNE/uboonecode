#include "WaveformData.hh"
#include <cmath>
#include <iostream>

namespace subevent {

  WaveformData::WaveformData() {
    waveforms.clear();
  }

  WaveformData::~WaveformData() {}

  void WaveformData::set( int ch, std::vector< double >& wfm, bool islowgain ) { 
    waveforms[ch] = std::vector< double >( wfm.begin(), wfm.end() );
    channels.insert( ch );
    is_lowgain_channel[ch] = islowgain;
  }
  
  void WaveformData::setLowGain( int ch, bool islowgain ) {
    is_lowgain_channel[ch] = islowgain;
  }

  bool WaveformData::isLowGain( int ch ) {
    if ( is_lowgain_channel.find( ch )!=is_lowgain_channel.end() )
      return is_lowgain_channel[ch];
    return false;
  }

  void WaveformData::calcBaselineInfo() {
    for ( ChannelSetIter ich=chbegin(); ich!=chend(); ich++ ) {

      std::vector<double>& wfm = waveforms[*ich];
      std::vector<double> mean( wfm.size(), 0.0 );
      std::vector<double> var( wfm.size(), 0.0 );
      // first five
      for ( int i=0; i<5; i++ ) {
	double x=0;
	double xx = 0;
	for (int j=0; j<10; j++) {
	  x += wfm.at( i+j );
	  xx += wfm.at( i+j )*wfm.at( i+j );
	}
	x /= 10.0;
	xx /= 10.0;
	mean.at( i ) = x;
	var.at( i ) = sqrt( xx - x*x );
      }
      // middle
      for ( int i=5; i<(int)wfm.size()-5; i++) {
	double x=0;
	double xx=0;
	for (int j=-5;j<5; j++) {
	  x += wfm.at( i+j );
	  xx += wfm.at( i+j )*wfm.at( i+j );
	}
	x /= 10.0;
	xx /= 10.0;
	mean.at( i ) = x;
	var.at( i ) = sqrt( xx-x*x );
      }
      // end
      for ( int i=wfm.size()-5; i<(int)wfm.size(); i++ ) {
	double x=0;
	double xx=0;
	for (int j=-10; j<0; j++) {
	  x += wfm.at( i+j );
	  xx += wfm.at( i+j )*wfm.at( i+j );
	}
	x /= 10.0;
	xx /= 10.0;
	mean.at( i ) = x;
	var.at( i ) = sqrt( xx - x*x );
      }

      // now we chuck any spots that have high variance connecting them linearly between good baselines
      
      // get first good mean and replace
      //std::cout << "a" << std::endl;
      double varthreshold = 1.0;
      int j=0; 

      while ( j<(int)var.size() && var.at(j)>varthreshold )
	j++;
	
      if ( j<(int)mean.size() ) {
	for (int i=0; i<j; i++) {
	  mean.at(i) = mean.at(j);
	  var.at(i) = 0.0;
	}
      }
      
      //std::cout << "b" << std::endl;

      int lastmean = j;
      for ( int i=j; i<(int)mean.size(); i++ ) {
	if ( var.at(i)>1.0 ) {
	  // find the next good mean
	  double k=i;
	  while ( k<(int)mean.size() && var.at(k)>1.0 )
	    k++;
	  if ( k<(int)mean.size() ) {
	    double dydx = ( mean.at(k)-mean.at(i-1) )/double( k-(i-1) );
	    for (int l=i; l<k; l++) {
	      var.at(l) = 0.0;
	      mean.at(l) = mean.at(i-1) + dydx*( l-(i-1) );
	    }
	    lastmean = i;
	  }
	}
	else
	  lastmean = i;
      }

      //std::cout << "c" << std::endl;
      
      for ( int i=lastmean; i<(int)mean.size(); i++ ) {
	mean.at( i ) = mean.at( lastmean );
      }

      //std::cout << "d" << std::endl;

      rollingmean[ *ich ] = mean;
      rollingvariance[ *ich ] = var;
      
    }//end of channel loop
  }//end of calcbaseline
}
