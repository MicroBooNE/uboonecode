#ifndef _UBOONETYPES_CHANNELDATAPMT_H
#define _UBOONETYPES_CHANNELDATAPMT_H
#include <memory>
#include <algorithm>
#include <sys/types.h>
#include <inttypes.h>
#include "evttypes.h"

#include <boost/serialization/list.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/binary_object.hpp>
#include <boost/serialization/map.hpp>

#include "constants.h"
#include "windowHeaderPMT.h"
#include "windowDataPMT.h"

namespace gov {
namespace fnal {
namespace uboone {
namespace datatypes {

using namespace gov::fnal::uboone;
 
/***
 *  Note: this is the serialization class that handles the card data.
 ***/

struct compareWindowHeaderPMT {
  bool operator() ( windowHeaderPMT lhs, windowHeaderPMT rhs) const
  {  
    if(lhs.getFrame()==rhs.getFrame())
      return lhs.getSample() < rhs.getSample();
    
    return lhs.getFrame() < rhs.getFrame();
  }
};

class channelDataPMT {

 public:
  static const uint8_t DAQ_version_number = gov::fnal::uboone::datatypes::constants::VERSION;
  
  channelDataPMT()
    {window_map.clear();}
  
  void insertWindow(windowHeaderPMT wH, windowDataPMT wD)
  { window_map.insert(std::pair<windowHeaderPMT,windowDataPMT>(wH,wD)); }

  void clearWindows(){ window_map.clear(); }

  int getNumberOfWindows() { return window_map.size(); }
  std::map<windowHeaderPMT,windowDataPMT, compareWindowHeaderPMT> getWindowMap() { return window_map; }

 private:
  std::map<windowHeaderPMT,windowDataPMT,compareWindowHeaderPMT> window_map;

  friend class boost::serialization::access;
  
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      if(version>0) {
	ar & window_map;
      }
    }

};


}  // end of namespace datatypes
}  // end of namespace uboone
}  // end of namespace fnal
}  // end of namespace gov

// This MACRO must be outside any namespaces.
BOOST_CLASS_VERSION(gov::fnal::uboone::datatypes::channelDataPMT, gov::fnal::uboone::datatypes::constants::VERSION)    

#endif /* #ifndef BOONETYPES_H */



