#ifndef __SubEventList__
#define __SubEventList__

#ifdef __BUILD_ROOT_DICT__
#include "TObject.h"
#endif

#include "SubEvent.hh"
#include <vector>

namespace subevent {

  typedef std::vector< SubEvent >::iterator SubEventListIter;

#ifdef __BUILD_ROOT_DICT__
  class SubEventList : public TObject { 
#else
  class SubEventList {
#endif

  public:
    SubEventList();
    ~SubEventList();

#ifndef __CINT__    
#ifndef __GCCXML__
    int add( SubEvent&& opflash );
#endif
#endif
    SubEvent& get( int i );
    SubEventListIter begin();
    SubEventListIter end();
    void sortByTime();
    void sortByCharge();
    void sortByAmp();
    int size() { return fSubEvents.size(); };
    void clear() { fSubEvents.clear(); fSubEvents.reserve(20); };
    bool sortedByTime() { if (sortMethod==kByTime) return true; else return false; }; 
    bool sortedByCharge() { if (sortMethod==kByCharge) return true; else return false; }; 
    bool sortedByAmp() { if (sortMethod==kByAmp) return true; else return false; }; 

  protected:
    std::vector< SubEvent > fSubEvents;

    typedef enum { kUnsorted=-1, kByTime, kByCharge, kByAmp } SortMethod_t;
    SortMethod_t sortMethod;

    static bool compareTime( SubEvent& t1, SubEvent& t2 ) {
      if (t1.tstart_ns<t2.tstart_ns )
	return true;
      else
	return false;
    };
    static bool compareArea( SubEvent& q1, SubEvent& q2 ) {
      if ( q1.totpe<q2.totpe ) return true;
      else return false;
    };
    static bool compareAmp( SubEvent& amp1, SubEvent& amp2 ) {
      if ( amp1.maxamp<amp2.maxamp ) return true;
      else return false;
    };
    
#ifdef __BUILD_ROOT_DICT__
    ClassDef( SubEventList, 1 )
#endif

  };

}




#endif
