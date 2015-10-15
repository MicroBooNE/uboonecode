#include "CosmicWindowSubEvents.hh"
#include <iostream>
#include "SubEvent.hh"
#include "Flash.hh"
#include "SubEventList.hh"
#include "FlashList.hh"
#include "scintresponse.hh"

namespace subevent {

  void formCosmicWindowSubEvents( CosmicWindowHolder& cosmicwindows, SubEventModConfig& config, SubEventList& subevents ) {
    cosmicwindows.sort(); // sort into (sample_time,channel) order
    int matchid = 0;
    for ( CosmicWindowIndexIter it=cosmicwindows.hgindexbegin(); it!=cosmicwindows.hgindexend(); it++ ) {
      if ( (*it).claimed )
	continue;
      CosmicWinIndex& first = (*it);
      first.claimed = true;
      // now collect all of sametime
      std::vector< CosmicWinIndex > matches;
      matches.push_back( first );
      for ( CosmicWindowIndexIter itb=it+1; itb!=cosmicwindows.hgindexend(); itb++ ) {
	if (  abs(first.sample-(*itb).sample) < config.flashgate ) {
	  matches.push_back( *itb );
	  (*itb).claimed = true;
	}
      }
      // make subevent from all of these
      SubEvent asubevent;
      asubevent.tstart_sample = (*matches.begin()).sample;
      asubevent.tend_sample = (*matches.begin()).sample;
      for ( CosmicWindowIndexIter itb=matches.begin(); itb!=matches.end(); itb++ ) {


	// get waveform
	std::vector< double >& wfm = cosmicwindows.highGainWfmMap[ (*itb) ];
	// get waveform quantities
	int winlen = (int)wfm.size();
	int maxt = 0;
	double maxamp = 0.;
	for ( int iadc=0; iadc<winlen; iadc++ ) {
	  if ( wfm.at(iadc)>maxamp ) {
	    maxamp = wfm.at(iadc);
	    maxt = iadc;
	  }
	}
	// if saturated, we use the low gain waveform!
	if ( maxamp>2000.0 ) {
	  //std::cout << "use low gain for ch=" << (*itb).channel << ", t=" << (*itb).sample << " maxamp=" << maxamp << std::endl;
	  CosmicWinIndex lgindex = cosmicwindows.getLGindexFromHG((*itb));
	  //std::cout << "  lowgain index=(" << lgindex.sample << ", " << lgindex.channel << ")" << std::endl;
	  wfm = cosmicwindows.lowGainWfmMap[ lgindex ];
	  winlen = (int)wfm.size();
	  maxamp = 0.;
	  for ( int iadc=0; iadc<winlen; iadc++ ) {
	    if ( wfm.at(iadc)>maxamp ) {
	      maxamp = wfm.at(iadc);
	      maxt = iadc;
	    }
	  }
	}
	
	//std::cout << "matches id=" << matchid << ": (" << (*itb).sample << ", " << (*itb).channel << ")  maxamp=" << maxamp << std::endl;
	if ( maxamp < 0.0 ) {
	  // weird, skip
	  (*itb).claimed = true;
	  //std::cout << "  undershoot waveform. skipping" << std::endl;
	  continue;
	}
	// calculate expectation
	std::vector< double > expectation;
	expectation.reserve(200);
	calcScintResponseCPP( expectation, 0, 10000, maxt, config.spe_sigma, maxamp, 
			      config.fastconst_ns, config.slowconst_ns, config.nspersample, config.fastfraction, config.slowfraction, config.noslowthreshold );
	// update start and end times of subevent
	if ( asubevent.tstart_sample > (*itb).sample ) asubevent.tstart_sample = (*itb).sample;
	if ( asubevent.tend_sample < (*itb).sample+(int)expectation.size() )
	  asubevent.tend_sample = (*itb).sample+(int)expectation.size();
	// make flash and add to subevent
	asubevent.flashes.add( Flash( (*itb).channel, (*itb).sample, (*itb).sample+(int)expectation.size(), (*itb).sample+maxt, maxamp, expectation, wfm ) );
      }
      subevents.add( std::move( asubevent ) );
      matchid++;
    }
    std::cout << "Formed " << subevents.size() << " cosmic discriminator subevents" << std::endl;
  }
  
}
