#ifndef _UBOONETYPES_GPS_H
#define _UBOONETYPES_GPS_H
#include <sys/types.h>
#include <inttypes.h>
#include "evttypes.h"
#include "share/boonetypes.h"

#include <boost/serialization/list.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>

#include "constants.h"

namespace gov {
namespace fnal {
namespace uboone {
namespace datatypes {

using namespace gov::fnal::uboone;

/**
   Note: this is the serialization class for the otherwise hardcoded gps_t struct,
   located in share/boonetypes. IF changes are made to gps_t, the appropriate changes
   should be made here as well, and the version number should be increased.
 **/

class gps {
 
public:
  static const uint8_t DAQ_version_number = gov::fnal::uboone::datatypes::constants::VERSION;
  gps() {};
  gps(uint32_t _lower, uint32_t _upper){ bt_gps.lower=_lower; bt_gps.upper=_upper; }

  uint32_t getLower() { return bt_gps.lower; }
  uint32_t getUpper() { return bt_gps.upper; }

  void setLower(uint32_t lower) { bt_gps.lower = lower; }
  void setUpper(uint32_t upper) { bt_gps.upper = upper; }

 private:
  gps_t bt_gps;
  friend class boost::serialization::access;
  
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      if( version > 0 )
	ar & bt_gps.lower & bt_gps.upper;
    }    
}; 


}  // end of namespace sebs
}  // end of namespace uboone
}  // end of namespace fnal
}  // end of namespace gov

// This MACRO must be outside any namespaces.
BOOST_CLASS_VERSION(gov::fnal::uboone::datatypes::gps, gov::fnal::uboone::datatypes::constants::VERSION)    

#endif /* #ifndef BOONETYPES_H */



