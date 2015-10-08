#include "FlashList.hh"

ClassImp( subevent::FlashList );

namespace subevent {

  FlashList::FlashList() {
    fFlashes.reserve(10);
  }
  FlashList::~FlashList() {}

  int FlashList::add( Flash&& opflash ) {
    fFlashes.emplace_back( opflash );
    return fFlashes.size();
  }

  Flash& FlashList::get( int i ) {
    return fFlashes.at(i);
  }

  FlashListIter FlashList::begin() {
    return fFlashes.begin();
  }
  
  FlashListIter FlashList::end() {
    return fFlashes.end();
  }

  void FlashList::sortByTime() {
    std::sort( begin(), end(), FlashList::compareTime );
    sortMethod = kByTime;
  }

  void FlashList::sortByCharge() {
    std::sort( begin(), end(), FlashList::compareArea );
    sortMethod = kByCharge;
  }

  void FlashList::sortByAmp() {
    std::sort( begin(), end(), FlashList::compareAmp );
    sortMethod = kByAmp;
  }

}
