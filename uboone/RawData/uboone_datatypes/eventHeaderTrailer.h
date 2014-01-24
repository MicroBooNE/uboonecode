#ifndef _UBOONETYPES_EVENTHEADERTRAILER_H
#define _UBOONETYPES_EVENTHEADERTRAILER_H
#include <sys/types.h>
#include <inttypes.h>
#include "evttypes.h"
#include "share/boonetypes.h"

#include <boost/serialization/list.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>

#include "constants.h"

/**
   Note: these are the serialization classes for the otherwise hardcoded event header and
   trailer structs, located in share/boonetypes. IF changes are made to gps_t, the 
   appropriate changes should be made here as well, and the version number should be 
   increased.
 **/

namespace gov {
namespace fnal {
namespace uboone {
namespace datatypes {

using namespace gov::fnal::uboone;

class eventHeader {
 
public:
  static const uint8_t DAQ_version_number = gov::fnal::uboone::datatypes::constants::VERSION;
  eventHeader() {};

  uint32_t getHeader() { return bt_event_header.header; }
  void setHeader(uint32_t header) { bt_event_header.header = header; }

  void setEventHeader(event_header_t EH) { bt_event_header = EH; }

 private:
  event_header_t bt_event_header;
  friend class boost::serialization::access;
  
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      if( version > 0 )
	ar & bt_event_header.header;
    }    
}; 


class eventTrailer {
 
public:
  static const uint8_t DAQ_version_number = gov::fnal::uboone::datatypes::constants::VERSION;
  eventTrailer() {};

  uint32_t getTrailer() { return bt_event_trailer.trailer; }
  void setTrailer(uint32_t trailer) { bt_event_trailer.trailer = trailer; }

  void setEventTrailer(event_trailer_t ET) { bt_event_trailer = ET; }

 private:
  event_trailer_t bt_event_trailer;
  friend class boost::serialization::access;
  
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      if( version > 0 )
	ar & bt_event_trailer.trailer;
    }    
}; 


}  // end of namespace sebs
}  // end of namespace uboone
}  // end of namespace fnal
}  // end of namespace gov

// This MACRO must be outside any namespaces.
BOOST_CLASS_VERSION(gov::fnal::uboone::datatypes::eventHeader, gov::fnal::uboone::datatypes::constants::VERSION)    
BOOST_CLASS_VERSION(gov::fnal::uboone::datatypes::eventTrailer, gov::fnal::uboone::datatypes::constants::VERSION)    
#endif /* #ifndef */



