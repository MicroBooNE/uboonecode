#ifndef __WAVEFORMDATA__
#define __WAVEFORMDATA__

#include <vector>
#include <map>
#include <set>

namespace subevent {

  typedef std::set<int>::iterator ChannelSetIter;

  class WaveformData {

  public:
    WaveformData();
    ~WaveformData();
    
    std::map< int, std::vector<double> > waveforms;
    std::map< int, bool > is_lowgain_channel;
    std::set< int > channels;
    ChannelSetIter chbegin() { return channels.begin(); };
    ChannelSetIter chend() { return channels.end(); };
    int nchannels() { return channels.size(); };

    std::vector< double >& get( int ch ) { return waveforms[ch]; };
    void set( int ch, std::vector<double>& wfm, bool islowgain );
    void setLowGain( int ch, bool islowgain );
    bool isLowGain( int ch );
    
  };

}

#endif
