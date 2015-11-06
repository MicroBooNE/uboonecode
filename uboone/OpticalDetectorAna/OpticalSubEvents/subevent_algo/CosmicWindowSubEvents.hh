#ifndef __CosmicWindowSubEvents__
#define __CosmicWindowSubEvents__

#include "SubEventModConfig.hh"
#include "SubEventList.hh"

#include <map>
#include <vector>
#include <algorithm>
#include <cstdlib>

namespace subevent {

  class CosmicWinIndex {
  public:
    CosmicWinIndex( int sample_, int channel_ ) {
      sample = sample_;
      channel = channel_;
      claimed = false;
    };
    ~CosmicWinIndex() {};
    int sample;
    int channel;
    bool claimed;

    bool operator==(const CosmicWinIndex& b ) const {
      if ( sample==b.sample && channel==b.channel ) 
	return true;
      return false;
    };
    bool operator<(const CosmicWinIndex& b ) const {
      if ( sample < b.sample )
	return true;
      else if ( sample > b.sample )
	return false;
      else {
	// if at sample time, use channel
	if ( channel < b.channel )
	  return true;
	else if ( channel > b.channel )
	  return false;
      }

      return false;
    };
  };

  
  typedef std::vector< CosmicWinIndex >::iterator CosmicWindowIndexIter;
  typedef std::map< CosmicWinIndex, std::vector< double > >::iterator CosmicWfmMapIter;

  class CosmicWindowHolder {

  public:
    CosmicWindowHolder() {};
    ~CosmicWindowHolder() {};

    std::vector< CosmicWinIndex > indexHG;
    std::vector< CosmicWinIndex > indexLG;
    std::map< CosmicWinIndex, std::vector< double > > highGainWfmMap;
    std::map< CosmicWinIndex, std::vector< double > > lowGainWfmMap;
    void addHG( int ch, int t_sample, const std::vector< double >& wfm ) { 
      indexHG.push_back( CosmicWinIndex( t_sample, ch ) );
      highGainWfmMap[ CosmicWinIndex( t_sample, ch ) ] = std::vector< double >( wfm.begin(), wfm.end() );
    };
    void addLG( int ch, int t_sample, const std::vector< double >& wfm ) { 
      indexLG.push_back( CosmicWinIndex( t_sample, ch ) );
      lowGainWfmMap[ CosmicWinIndex( t_sample, ch ) ] = wfm;
    };
    void sort() {
      std::sort( indexHG.begin(), indexHG.end() );
      std::sort( indexLG.begin(), indexLG.end() );
    };
    CosmicWinIndex getLGindexFromHG( CosmicWinIndex hgindex ) {
      CosmicWinIndex lgindex(-1,-1);
      int closest = -10000;
      for ( CosmicWindowIndexIter it=lgindexbegin(); it!=lgindexend(); it++ ) {
	if ( (*it).channel!=hgindex.channel )
	  continue;
	if ( closest<0 || closest < abs( (*it).sample-hgindex.sample ) ) {
	  closest = abs( (*it).sample-hgindex.sample );
	  lgindex = (*it);
	}
      }
      return lgindex;
    };
    CosmicWindowIndexIter hgindexbegin() { return indexHG.begin(); };
    CosmicWindowIndexIter hgindexend() { return indexHG.end(); };
    CosmicWindowIndexIter lgindexbegin() { return indexLG.begin(); };
    CosmicWindowIndexIter lgindexend() { return indexLG.end(); };

  };

  void formCosmicWindowSubEvents( CosmicWindowHolder& cosmicwindows, SubEventModConfig& config, SubEventList& subevents );
}

#endif
