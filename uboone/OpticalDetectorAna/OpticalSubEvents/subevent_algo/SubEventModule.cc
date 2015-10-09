#include "SubEventModule.hh"
#include "Flash.hh"
#include "FlashList.hh"
#include "SubEventList.hh"
#include "WaveformData.hh"
#include "uboone/OpticalDetectorAna/OpticalSubEvents/cfdiscriminator_algo/cfdiscriminator.hh"
#include "scintresponse.hh"
#include <algorithm>
#include <iostream>
#include <cmath>

namespace subevent {

  int findChannelFlash( int channel, std::vector< double >& waveform, SubEventModConfig& config, Flash& opflash ) {
    // ---------------------------------------
    // input
    // int channel: channel id
    // waveform: ADCs
    // config: SubEventModule configuration
    // output
    // opflash: Flash object
    // ---------------------------------------

    std::vector< int > t_fire;
    std::vector< int > amp_fire;
    std::vector< int > maxt_fire;
    std::vector< int > diff_fire;
    //std::cout << "CFD config: " << config.cfdconfig.delay << " " <<  config.cfdconfig.threshold << " " <<  config.cfdconfig.deadtime << " " <<  config.cfdconfig.width << std::endl;
    
    cpysubevent::runCFdiscriminatorCPP( t_fire, amp_fire, maxt_fire, diff_fire, waveform.data(), 
					config.cfdconfig.delay, config.cfdconfig.threshold, config.cfdconfig.deadtime, config.cfdconfig.width, waveform.size() );

    // find largest
    int largestCFD = -1;
    double maxamp = 0;
    for (int n=0; n<(int)t_fire.size(); n++) {
      if ( amp_fire.at(n) > maxamp ) {
	maxamp = amp_fire.at(n);
	largestCFD = n;
      }
    }

    if ( largestCFD==-1 )
      return 0;

    // store info from discriminator
    opflash.ch = channel;
    opflash.tstart = std::max(t_fire.at( largestCFD )-config.npresamples, 0);
    opflash.tmax = maxt_fire.at( largestCFD );
    opflash.maxamp = amp_fire.at( largestCFD );

    // calc scint response
    std::vector< double > expectation;
    expectation.reserve( 200 );
    subevent::calcScintResponseCPP( expectation, 
				    opflash.tstart, (int)waveform.size(), opflash.tmax, 
				    config.spe_sigma, maxamp-waveform.at( opflash.tstart ), config.fastconst_ns, config.slowconst_ns, config.nspersample );

    opflash.tend = std::min( opflash.tstart+(int)expectation.size(), (int)waveform.size()-1 );
    
    std::vector< double > subwfm( waveform.begin()+opflash.tstart, waveform.begin()+opflash.tend );
    opflash.storeWaveform( subwfm );
    std::vector< double > subexp( expectation.begin(), expectation.begin()+(int)subwfm.size() );
    opflash.storeExpectation(  subexp );

    // for debug
    //std::cout << "return opflash: " << opflash.ch << " " << opflash.tstart << " " << opflash.tend << " " << opflash.tmax << " " << opflash.maxamp << std::endl;
    // std::cout << " expectation of flash: ";
    // for ( int i=0; i<(int)expectation.size(); i++)
    // std::cout << expectation.at(i) << " ";
    // std::cout << std::endl;

    return 1;
  }

  int getChannelFlashes( int channel, std::vector< double >& waveform, SubEventModConfig& config, FlashList& flashes, std::vector<double>& postwfm ) {
    // corresponds to cyRunSubEventDiscChannel
    // input
    // channel: FEMCH number
    // waveform: ADCs
    // config
    // output

    // make our working copy of the waveform
    postwfm.clear();
    postwfm.reserve( waveform.size() );
    std::copy( waveform.begin(), waveform.end(), back_inserter(postwfm) );

    //  find subevent
    int maxsubevents = config.maxchflashes;
    //std::cout << "  maxsubevents=" << maxsubevents << std::endl;
    int nsubevents = 0;
    double fx = 0.0;
    double sig = 0.0;
    double thresh = 0.0;
    double chped = 0.;
    // for (int i=0; i<5; i++)
    //   chped += waveform.at(i);
    // chped /= 5.0;
    
    //flashes.clear();

    while ( nsubevents<maxsubevents ) {
      // find one subevent (finds the largest flash of light)
      Flash opflash;
      int found = findChannelFlash( channel, postwfm, config, opflash );
      //std::cout << "[getChannelFlashes] Found  " << found << " flashes in channel " << channel << std::endl;
      if ( found==0 )
	break;
      
      // subtract waveform below subevent threshold
      //double amp_start = waveform.at( opflash.tstart );
      //double amp_end   = waveform.at( opflash.tend );
      //double slope = (amp_end-amp_start)/( opflash.tend-opflash.tstart );
      opflash.area30 = 0.0;
      opflash.area = 0.0;
      for (int tdc=0; tdc<(int)opflash.expectation.size(); tdc++) {
	fx = opflash.expectation.at(tdc);
	sig = sqrt( fx/20.0 );
	thresh = fx + 3.0*sig*20.0; // 3 sigma variance
	if ( postwfm.at( opflash.tstart )-chped < thresh ) {
	  postwfm.at( opflash.tstart + tdc ) = chped;
	}
	if ( tdc<600 && opflash.tstart+tdc+20<(int)waveform.size() ) {
	  if ( tdc<30 )
	    opflash.area30 += waveform.at( opflash.tstart+tdc )-chped;
	  opflash.area += waveform.at( opflash.tstart+tdc )-chped;
	}
      }
      opflash.fcomp_gausintegral = (opflash.maxamp-chped)*(config.spe_sigma/15.625)*sqrt(2.0)*3.14159;
      nsubevents += 1;
      flashes.add( std::move(opflash) );
    }//end of subflash search
      
    return nsubevents;

  }

  void formFlashes( WaveformData& wfms, SubEventModConfig& config, FlashList& flashes ) {

    for ( ChannelSetIter it=wfms.chbegin(); it!=wfms.chend(); it++ ) {
      int ch = *it;
      std::vector< double > postwfm;
      getChannelFlashes( ch, wfms.get( ch ), config, flashes, postwfm );
      std::cout << "search for flashes in channel=" << ch << ". found=" << flashes.size() << std::endl;
    }
  }

  void fillFlashAccumulators( FlashList& flashes, std::map< int, double >& pmtspemap, SubEventModConfig& config, std::vector< double >& peacc, std::vector< double >& hitacc ) {

    for ( FlashListIter iflash=flashes.begin(); iflash!=flashes.end(); iflash++ ) {
      if ( (*iflash).claimed )
	continue;

      int start = std::max( int( (*iflash).tstart-0.5*config.flashgate), 0 );
      int end = std::min( int( (*iflash).tstart+0.5*config.flashgate ), (int)peacc.size() );
      //std::cout << "add flash acc: ch=" << (*iflash).ch << " maxamp=" <<  ((*iflash).maxamp)/pmtspemap[(*iflash).ch] << " t=[" << start << ", " << end << "]" << std::endl;
      for ( int t=start; t<end; t++ ) {
	peacc.at(t) += ((*iflash).maxamp)/pmtspemap[(*iflash).ch];
	hitacc.at(t) += 1.0;
      }
      
    }

  }

  void formSubEvents( WaveformData& wfms, SubEventModConfig& config, std::map< int, double >& pmtspemap, SubEventList& subevents ) {

    std::cout << "FormSubEvents" << std::endl;

    FlashList flashes;
    formFlashes( wfms, config, flashes );
    std::cout << "  total flashes: " << flashes.size() << std::endl;

    int nloops = 0;
    ChannelSetIter itch=wfms.chbegin();
    int nsamples = wfms.get( *itch ).size();
    std::vector< double > peacc( nsamples, 0.0 );
    std::vector< double > hitacc( nsamples, 0.0 );

    while ( nloops < config.maxsubeventloops ) {
      std::cout << " start subevent search: loop#" << nloops << std::endl;

      // accumulators: the summed pulse height and the number of hits
      // use first entry to set size
      peacc.assign( nsamples, 0.0 );
      hitacc.assign( nsamples, 0.0 );
      
      fillFlashAccumulators( flashes, pmtspemap, config, peacc, hitacc );

      // find maximums
      //double  hit_tmax = 0;
      double pe_tmax = 0;
      double pemax = 0;
      double hitmax = 0;
      for ( int tick=0; tick<(int)peacc.size(); tick++ ) {
	if ( peacc.at(tick)>pemax ) {
	  pemax = peacc.at(tick);
	  pe_tmax = tick;
	}
	if ( hitacc.at(tick)>hitmax ) {
	  hitmax = hitacc.at(tick);
	  //hit_tmax = tick;
	}
      }
      std::cout << "  accumulator max: t=" << pe_tmax << " amp=" << pemax << " hits=" << hitmax << std::endl;

      // organize flashes within maxima
      if ( pemax>config.ampthresh || hitmax>config.hitthresh ) {
	// passed! 
	SubEvent newsubevent;
	
	//if ( !flashes.sortedByTime()  ) flashes.sortByTime();
	//std::cout << "   sorted flashes by time" << std::endl;
	
	// form subevent by grouping flashes around tmax
	int nclaimed = 0;
	for ( FlashListIter iflash=flashes.begin(); iflash!=flashes.end(); iflash++ ) {
	  
	  if ( (*iflash).claimed ) continue;
	  
	  if ( abs( (*iflash).tstart - pe_tmax )< config.flashgate ) {

	    newsubevent.tstart_sample = (int)pe_tmax;
	    if ( newsubevent.tend_sample < (int)(*iflash).tend )
	      newsubevent.tend_sample = (int)(*iflash).tend;
	    newsubevent.maxamp = pemax;
	    //newsubevent.totpe += (*iflash).area/pmtspemap[ (*iflash).ch ];
	    newsubevent.totpe += (*iflash).area; // HACK
	    newsubevent.sumflash30 += ((*iflash).area30); // HACK
	    newsubevent.sumfcomp_gausintegral += (*iflash).fcomp_gausintegral; // HACK
	    (*iflash).claimed = true;
	    Flash copyflash( (*iflash ) );
	    newsubevent.flashes.add( std::move( copyflash ) ); 
	    nclaimed++;
	  }

	} //end of flash loop
	
	// store new subevent
	std::cout << "  subevent " << subevents.size() << ": tstart=" << newsubevent.tstart_sample << "  tend=" << newsubevent.tend_sample << " nflashes=" << newsubevent.flashes.size() << std::endl;
	subevents.add( std::move( newsubevent ) );
	
      }
      else {
	std::cout << "  did not find additional subevent" << std::endl;
	break;
      }

      nloops += 1;
    }//end of while loop
    
    std::cout << " end of formsubevents. found " << subevents.size() << std::endl;
    //std::cin.get();
  }

}
