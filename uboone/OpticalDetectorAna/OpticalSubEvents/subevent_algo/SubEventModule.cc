#include "SubEventModule.hh"
#include "Flash.hh"
#include "FlashList.hh"
#include "SubEventList.hh"
#include "WaveformData.hh"
#include "uboone/OpticalDetectorAna/OpticalSubEvents/cfdiscriminator_algo/cfdiscriminator.hh"
//#include "cfdiscriminator.hh"
#include "scintresponse.hh"
#include <algorithm>
#include <iostream>
#include <cmath>

namespace subevent {

  // ============================================================================================================
  // findChannelFlash
  int findChannelFlash( int channel, std::vector< double >& waveform, SubEventModConfig& config, std::string discrname, Flash& opflash ) {
    // ---------------------------------------
    // input
    // int channel: channel id
    // waveform: ADCs
    // config: SubEventModule configuration
    // chooses which CFD setting to use
    // output
    // opflash: Flash object
    // ---------------------------------------

    std::vector< int > t_fire;
    std::vector< int > amp_fire;
    std::vector< int > maxt_fire;
    std::vector< int > diff_fire;
    //std::cout << "CFD config: " << config.cfdconfig.delay << " " <<  config.cfdconfig.threshold << " " <<  config.cfdconfig.deadtime << " " <<  config.cfdconfig.width << std::endl;
    
    if ( discrname!="pass2" )
      cpysubevent::runCFdiscriminatorCPP( t_fire, amp_fire, maxt_fire, diff_fire, waveform.data(), 
					  config.cfdconfig.delay, config.cfdconfig.threshold, config.cfdconfig.deadtime, config.cfdconfig.width, waveform.size() );
    else
      cpysubevent::runCFdiscriminatorCPP( t_fire, amp_fire, maxt_fire, diff_fire, waveform.data(), 
					  config.cfdconfig_pass2.delay, config.cfdconfig_pass2.threshold, config.cfdconfig_pass2.deadtime, config.cfdconfig_pass2.width, waveform.size() );
	

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
				    config.spe_sigma, maxamp-waveform.at( opflash.tstart ), config.fastconst_ns, config.slowconst_ns, 
				    config.nspersample, config.fastfraction, config.slowfraction, config.noslowthreshold );

    opflash.tend = std::min( opflash.tstart+(int)expectation.size(), (int)waveform.size() );
    
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

  // ============================================================================================================
  // getChannelFlashes
  int getChannelFlashes( int channel, std::vector< double >& waveform, SubEventModConfig& config, std::string discrname, FlashList& flashes, std::vector<double>& postwfm ) {
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
      int found = findChannelFlash( channel, postwfm, config, discrname, opflash );
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
	// ---------------------------------------------
	// This is all a little hacky
	// this variance threshold should be tunned better?
	// it really serves as a first pass attempt at preventing repeated flash formation
	// later, we will refine the subevent to look for subevents inside the late-light tail
	// but really, someone please put xenon into the detector
	sig = sqrt( fx/20.0 );
	thresh = fx + 3.0*sig*20.0; // 3 sigma variance
	if ( postwfm.at( opflash.tstart )-chped < thresh ) {
	  postwfm.at( opflash.tstart + tdc ) = chped;
	}
	// we cap the area calculation and keep a 20 sample buffer from the end
	if ( tdc<600 && opflash.tstart+tdc+20<(int)waveform.size() ) { 
	  if ( tdc<30 )
	    opflash.area30 += waveform.at( opflash.tstart+tdc )-chped;
	  opflash.area += waveform.at( opflash.tstart+tdc )-chped;
	}
	// ---------------------------------------------
      }
      opflash.fcomp_gausintegral = (opflash.maxamp-chped)*(config.spe_sigma/15.625)*sqrt(2.0)*3.14159;
      nsubevents += 1;
      flashes.add( std::move(opflash) );
    }//end of subflash search
      
    return nsubevents;

  }

  // ============================================================================================================
  // formFlashes
  void formFlashes( WaveformData& wfms, SubEventModConfig& config, std::string discrname, FlashList& flashes, WaveformData& postwfms ) {
    for ( ChannelSetIter it=wfms.chbegin(); it!=wfms.chend(); it++ ) {
      int ch = *it;
      std::vector< double > postwfm;
      getChannelFlashes( ch, wfms.get( ch ), config, discrname, flashes, postwfm );
      //std::cout << "search for flashes in channel=" << ch << ". found=" << flashes.size() << std::endl;
      postwfms.set( ch, postwfm, wfms.isLowGain(ch) );
    }
  }

  // ============================================================================================================
  // fillFlashAccumulators
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

  // ============================================================================================================
  // formSubEvents
  void formSubEvents( WaveformData& wfms, SubEventModConfig& config, std::map< int, double >& pmtspemap, SubEventList& subevents, FlashList& unclaimed_flashes ) {

    //std::cout << "FormSubEvents" << std::endl;

    WaveformData postwfms; // this will store "post" waveforms. We remove sections of the waveforms used to build hits (by returning them to baseline)

    // We find flashes of light on each channel. We do this in two passes: a high threshold pass to look for pulses above any background spe light.
    FlashList flashes;
    formFlashes( wfms, config, "pass1", flashes, postwfms );
    //std::cout << "  total pass1/high-trehsold flashes: " << flashes.size() << std::endl;

    WaveformData postpostwfms;
    FlashList flashes_pass2;
    formFlashes( postwfms, config, "pass2", flashes_pass2, postpostwfms );
    //std::cout << "  total pass2/low-threshold flashes: " << flashes_pass2.size() << std::endl;

    int nloops = 0;
    ChannelSetIter itch=wfms.chbegin();
    int nsamples = wfms.get( *itch ).size();
    std::vector< double > peacc( nsamples, 0.0 );
    std::vector< double > hitacc( nsamples, 0.0 );

    while ( nloops < config.maxsubeventloops ) {
      //std::cout << " start subevent search: loop#" << nloops << std::endl;

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
      //std::cout << "  accumulator max: t=" << pe_tmax << " amp=" << pemax << " hits=" << hitmax << std::endl;

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
	    newsubevent.totpe += (*iflash).area;
	    newsubevent.sumflash30 += ((*iflash).area30); 
	    newsubevent.sumfcomp_gausintegral += (*iflash).fcomp_gausintegral;
	    (*iflash).claimed = true;
	    Flash copyflash( (*iflash ) );
	    newsubevent.flashes.add( std::move( copyflash ) ); 
	    nclaimed++;
	  }

	} //end of flash loop
	
	// store new subevent
	//std::cout << "  subevent " << subevents.size() << ": tstart=" << newsubevent.tstart_sample << "  tend=" << newsubevent.tend_sample << " nflashes=" << newsubevent.flashes.size() << std::endl;
	subevents.add( std::move( newsubevent ) );
	
      }
      else {
	//std::cout << "  did not find additional subevent" << std::endl;
	break;
      }

      nloops += 1;
    }//end of while loop
    
    //std::cout << " end of formsubevents. found " << subevents.size() << std::endl;

    // now add in second pass flashes
    int nadditions = 0;
    if ( subevents.size() > 0 ) {
      flashes_pass2.sortByTime();
      for( FlashListIter iflash2=flashes_pass2.begin(); iflash2!=flashes_pass2.end(); iflash2++ ) {
	if ( (*iflash2).claimed )
	  continue;
	if ( wfms.isLowGain( (*iflash2).ch ) && (*iflash2).maxamp < 40.0 )
	  continue; // skip low pulse height flashes in waveforms coming from low-gain. It is in the noise.
	
	double mintstart_diff = 1.0e6;
	int best_subevent = -1;
	for ( int isubevent=0; isubevent<(int)subevents.size(); isubevent++ ) {
	  double tdiff = fabs( subevents.get(isubevent).tstart_sample-(*iflash2).tstart );
	  if ( (*iflash2).tstart >= subevents.get(isubevent).tstart_sample-config.flashgate && (*iflash2).tend <= (subevents.get(isubevent).tend_sample)+config.flashgate ) {
	    if ( tdiff < mintstart_diff ) {
	      best_subevent = isubevent;
	      mintstart_diff = tdiff;
	    }
	  }
	}
	
	// add to subevent if a match found
	if ( best_subevent>=0 ) {
	  (*iflash2).claimed = true;
	  subevents.get(best_subevent).flashes_pass2.add( std::move(*iflash2) );
	  //std::cout << " adding second pass flash to subevent #" << best_subevent << std::endl;
	  nadditions += 1;
	}
	else {
	  // moved into unclaimed flashes list
	  unclaimed_flashes.add( std::move( *iflash2 ) );
	}
      }
    }
    //std::cout << " added " << nadditions << " 2nd pass flashes to the subevents." << std::endl;
    //std::cin.get();

    nadditions = 0;
    for ( int iflash=0; iflash<(int)flashes.size(); iflash++ ) {
      if ( !flashes.get( iflash ).claimed ) {
	unclaimed_flashes.add( std::move( flashes.get( iflash ) ) );
	nadditions ++;
      }
    }
    //std::cout << " added " << nadditions << " unclaimed 1st pass flashes to the subevents." << std::endl;

    // finally analyze subevents
    //std::cout << " analyze subevents" << std::endl;
    AnalyzeSubEvents( subevents );

    //std::cout << "subevents formed." << std::endl;
    
  }

  // ============================================================================================================
  // Analyze SubEvents
  void AnalyzeSubEvents( SubEventList& subevents ) {

    // We analyze each subevent, calculating quantities such as area
    for ( SubEventListIter isubevent=subevents.begin(); isubevent!=subevents.end(); isubevent++ ) {

      // subevent and timing
      SubEvent& subevent = (*isubevent);
      int width = subevent.tend_sample-subevent.tstart_sample;
      int tstart = subevent.tstart_sample;

      // combine flashes into one waveform
      std::vector< double > combined_flash( width, 0.0 );
      for ( FlashListIter iflash=subevent.flashes.begin(); iflash!=subevent.flashes.end(); iflash++ ) {
	Flash& aflash = (*iflash);
	int offset = aflash.tstart-tstart;
	for (int tdc=0; tdc<(int)aflash.waveform.size(); tdc++) {
	  if ( tdc+offset>=0 && tdc+offset<width )
	    combined_flash.at( tdc+offset ) += aflash.waveform.at( tdc ); // assumed that the pedestal has been subtracted -- should relax this later
	}
      }
      std::vector< double > combined_flash2( width, 0.0 );
      for ( FlashListIter iflash=subevent.flashes_pass2.begin(); iflash!=subevent.flashes_pass2.end(); iflash++ ) {
	Flash& aflash = (*iflash);
	int offset = aflash.tstart-tstart;
	for (int tdc=0; tdc<(int)aflash.waveform.size(); tdc++) {
	  if ( tdc+offset>=0 && tdc+offset<width )
	    combined_flash2.at( tdc+offset ) += aflash.waveform.at( tdc ); // assumed that the pedestal has been subtracted -- should relax this later
	}
      }
      

      // find tmax, maxamp
      subevent.maxamp = -1.0e-6;
      for ( int tdc=0; tdc<(int)combined_flash.size(); tdc++ ) {
	if ( subevent.maxamp < combined_flash[tdc] ) {
	  subevent.maxamp = combined_flash[tdc];
	  subevent.tmax_sample = tdc+tstart;
	}
      }

//       // find the start using the CFD, using a threshold of 10% maxamp
//       std::vector< int > t_fire;
//       std::vector< int > amp_fire;
//       std::vector< int > maxt_fire;
//       std::vector< int > diff_fire;
//       cpysubevent::runCFdiscriminatorCPP( t_fire, amp_fire, maxt_fire, diff_fire, combined_flash.data(), 4, 40.0, 24, 24, combined_flash.size() );

      // calculate total pe and pe in first 30 samples
      subevent.totpe = 0.0;
      subevent.pe30  = 0.0;
      subevent.totpe_1 = 0.0; // first pass prompt fraction
      subevent.pe30_1  = 0.0; // first pass prompt fraction
      int flash_start = std::max( 0, (subevent.tmax_sample-tstart)-10 );
//       if ( t_fire.size()>0 )
// 	flash_start = std::max( 0, t_fire.at(0)-4 );
      for ( int tdc=flash_start; tdc<(int)combined_flash.size(); tdc++ ) {
	subevent.totpe += combined_flash[tdc] + combined_flash2[tdc];
	subevent.totpe_1 += combined_flash[tdc];
	if ( tdc-flash_start<30 ) {
	  subevent.pe30   += combined_flash[tdc] + combined_flash2[tdc];
	  subevent.pe30_1 += combined_flash[tdc];
	}
      }
      
    }
  }

}
