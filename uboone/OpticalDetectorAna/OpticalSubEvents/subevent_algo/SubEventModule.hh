#ifndef __SubEventModule__
#define __SubEventModule__

#include <vector>
#include <map>
#include "SubEventModConfig.hh"
#include "SubEvent.hh"
#include "Flash.hh"
#include "FlashList.hh"
#include "SubEventList.hh"
#include "WaveformData.hh"

namespace subevent {
  
  // Main Routine ------------------------------------------------------------
  void formSubEvents( WaveformData& wfms, SubEventModConfig& config, std::map< int, double >& pmtspemap, SubEventList& subevents );
  // -------------------------------------------------------------------------

  int findChannelFlash( int ch, std::vector<double>& waveform, SubEventModConfig& config, Flash& returned_flash );
  int getChannelFlashes( int channel, std::vector< double >& waveform, SubEventModConfig& config, FlashList& flashes, std::vector<double>& postwfm );
  
  void formFlashes( WaveformData& wfms, SubEventModConfig& config, FlashList& flashes );
  void fillFlashAccumulators( FlashList& flashes, std::map< int, double >& pmtspemap, SubEventModConfig& config, std::vector< double >& peacc, std::vector< double >& hitacc );
  
}

#endif
