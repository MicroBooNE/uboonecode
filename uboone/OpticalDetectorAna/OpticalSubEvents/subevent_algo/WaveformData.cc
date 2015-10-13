#include "WaveformData.hh"

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
  
}
