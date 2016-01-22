#include "FEMBasicTriggerAlgo.h"
#include <iostream>

namespace fememu {

  void basicTrigger( int BeamWinSize, int NChannels, const BasicTriggerConfig& config, 
		     const std::vector< std::vector<int> >& chwfms, std::vector<int>& vmaxdiff, std::vector<int>& vmulti ) {

    //std::cout << __PRETTY_FUNCTION__ << std::endl;

    // first calculate accumulators for each waveform
    std::vector<int> chdiff[NChannels];
    std::vector<int> chhit[NChannels];
    for (int ch=0; ch<NChannels; ch++) {
      
      const std::vector<int>& wfm = chwfms[ch];
      
      chdiff[ch].resize( wfm.size(), 0 );
      chhit[ch].resize( wfm.size(), 0 );
      
      // memory for diff vectors
      std::vector<int> diff0( (int)wfm.size(), 0 );
      std::vector<int> diff3( (int)wfm.size(), 0 );
      for (int tick=config.fDiscr0delay; tick<(int)wfm.size(); tick++)
	diff0.at(tick) = wfm.at(tick)-wfm.at(tick-config.fDiscr0delay);
      for (int tick=config.fDiscr3delay; tick<(int)wfm.size(); tick++)
	diff3.at(tick) = wfm.at(tick)-wfm.at(tick-config.fDiscr3delay);


      // determine triggers and fill accumulators
      std::vector<int> ttrig0;
      std::vector<int> ttrig3;
    
      for (int tick=0; tick<(int)wfm.size(); tick++) {
	// discr0 must fire first: acts as pre-trigger. won't fire again until all discs are no longer active
	if ( diff0.at(tick)>=config.fDiscr0threshold ) {
	  if ( (ttrig0.size()==0 || ttrig0.at( ttrig0.size()-1 )+config.fDiscr0precount<tick )
	       && ( ttrig3.size()==0 || ttrig3.at( ttrig3.size()-1 )+config.fDiscr3deadtime < tick ) ) {
	    // form discr0 trigger
	    ttrig0.push_back( tick );
	  }
	} // end of if discr0 fires

	// discr3 fire
	if ( diff3.at(tick)>=config.fDiscr3threshold ) {
	  // must be within discr0 prewindow and outside of past discr3 deadtime
	  if ( ( ttrig0.size()>0 && tick-ttrig0.at( ttrig0.size()-1 ) < config.fDiscr0deadtime )
	       && ( ttrig3.size()==0 || ttrig3.at( ttrig3.size()-1 ) + config.fDiscr3deadtime < tick ) ) {
	    ttrig3.push_back( tick );
	    // find maxdiff
	    int tmaxdiff = diff3.at(tick);
	    for (int t=tick; t<std::min( tick+config.fDiscr3width, (int)diff3.size() ); t++) {
	      if ( tmaxdiff<diff3.at(t) )
		tmaxdiff = diff3.at(t);
	    }
	    // fill the accumulators
	    int tend = std::min( tick+config.fDiscr3deadtime, (int)diff3.size() );
	    for (int t=tick; t<tend; t++) {
	      chdiff[ch].at( t ) = tmaxdiff;
	      chhit[ch].at( t ) = 1;
	    }
	  }
	}
      }//end of wfm loop for trigger and accumulators
    }//end of channel loop

    // break up waveform into windows and calculate trigger vars for each window
    int wfmsize = (int)chwfms.at(0).size();
    if ( wfmsize < config.fMinReadoutTicks ) {
      std::cout << "Beam readout window size is too small! (" << wfmsize << " < " << config.fMinReadoutTicks << ")" << std::endl;
      return;
    }

    int nwindows = (wfmsize-1-config.fFrontBuffer)/BeamWinSize;
    vmaxdiff.clear();
    vmulti.clear();

    for (int iwin=0; iwin<nwindows; iwin++) {
      int winstart = config.fFrontBuffer + iwin*BeamWinSize;
      int winend   = config.fFrontBuffer + (iwin+1)*BeamWinSize;
      //winid = iwin;
      int winmaxmulti = 0;
      int winmaxdiff = 0;
      for (int tick=winstart; tick<winend; tick++) {
	int maxdiff_ = 0;
	int nhit_ = 0;
	for (int ch=0; ch<NChannels; ch++) {
	  maxdiff_ += chdiff[ch].at(tick);
	  nhit_    += chhit[ch].at(tick);
	}
	if ( winmaxdiff < maxdiff_ )
	  winmaxdiff = maxdiff_;
	if ( winmaxmulti < nhit_ )
	  winmaxmulti = nhit_;
      }

      // store for the window
      vmaxdiff.push_back( winmaxdiff );
      vmulti.push_back( winmaxmulti );
    
    }
  
  }

}
