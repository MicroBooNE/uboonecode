#include "WaveformData.hh"

namespace subevent {

  WaveformData::WaveformData() {
    waveforms.clear();
  }

  WaveformData::~WaveformData() {}

  void WaveformData::set( int ch, std::vector< double >& wfm ) { 
    waveforms[ch] = std::vector< double >( wfm.begin(), wfm.end() );
    channels.insert( ch );
  }

}
