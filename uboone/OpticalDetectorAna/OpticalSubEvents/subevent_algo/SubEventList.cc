#include "SubEventList.hh"
#include <algorithm>

#ifdef __BUILD_ROOT_DICT__
ClassImp( subevent::SubEventList )
#endif

namespace subevent {

  SubEventList::SubEventList() {
    fSubEvents.reserve(10);
  }
  SubEventList::~SubEventList() {}

  int SubEventList::add( SubEvent&& opflash ) {
    fSubEvents.emplace_back( opflash );
    return fSubEvents.size();
  }

  SubEvent& SubEventList::get( int i ) {
    return fSubEvents.at(i);
  }

  SubEventListIter SubEventList::begin() {
    return fSubEvents.begin();
  }
  
  SubEventListIter SubEventList::end() {
    return fSubEvents.end();
  }

  void SubEventList::sortByTime() {
    std::sort( begin(), end(), SubEventList::compareTime );
    sortMethod = kByTime;
  }

  void SubEventList::sortByCharge() {
    std::sort( begin(), end(), SubEventList::compareArea );
    sortMethod = kByCharge;
  }

  void SubEventList::sortByAmp() {
    std::sort( begin(), end(), SubEventList::compareAmp );
    sortMethod = kByAmp;
  }

}
