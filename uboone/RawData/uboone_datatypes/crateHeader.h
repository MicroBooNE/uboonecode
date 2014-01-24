#ifndef _UBOONETYPES_CRATEHEADER_H
#define _UBOONETYPES_CRATEHEADER_H
#include <memory>
#include <algorithm>
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

/***
   Note: this is the serialization class for the hardcoded crate_header_t struct
   that is used by the SEBs to attach a header per crate. If there is a change 
   there, it needs to be made here too. If there is a change made here, it may do
   no good without a change there.
***/


class crateHeader {
 public:
  static const uint8_t DAQ_version_number = gov::fnal::uboone::datatypes::constants::VERSION;
  crateHeader();
  crateHeader(crate_header_t cH) {  bt_crate_header = cH; }
  
  bool getCrateComplete() { return bt_crate_header.complete; }
  uint32_t getCrateSize() { return bt_crate_header.size; }
  uint8_t getCrateNumber() { return bt_crate_header.crate_number; }
  uint8_t getCardCount() { return bt_crate_header.card_count; }
  uint32_t getCrateEventNumber() { return bt_crate_header.event_number; }
  uint32_t getCrateFrameNumber() { return bt_crate_header.frame_number; }
  uint8_t getCrateType() { return (uint8_t)((bt_crate_header.crateBits >> 4) & 0xf ); }

  void setCrateComplete(bool value) { bt_crate_header.complete = value; }
  void setCrateSize(uint32_t size) { bt_crate_header.size = size; }
  void setCrateNumber(uint8_t cnum) { bt_crate_header.crate_number = cnum; }
  void setCardCount(uint8_t ccnt) { bt_crate_header.card_count = ccnt; }
  void setCrateEventNumber(uint32_t event) { bt_crate_header.event_number = event; }
  void setCrateFrameNumber(uint32_t frame) { bt_crate_header.frame_number = frame; }
  void setCrateType(uint8_t type) { bt_crate_header.crateBits = (bt_crate_header.crateBits & 0xf) + (type << 4); }

  crate_header_t getCrateHeader() { return bt_crate_header; }
  void setCrateHeader(crate_header_t cH) { bt_crate_header = cH; }

 private:
  crate_header_t bt_crate_header;

  friend class boost::serialization::access;
  
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      if(version==1)
	ar & bt_crate_header.complete
	   & bt_crate_header.size
	   & bt_crate_header.crate_number
	   & bt_crate_header.card_count
	   & bt_crate_header.event_number
 	   & bt_crate_header.frame_number;
      else if(version>1)
	ar & bt_crate_header.complete
	   & bt_crate_header.crateBits
	   & bt_crate_header.size
	   & bt_crate_header.crate_number
	   & bt_crate_header.card_count
	   & bt_crate_header.event_number
 	   & bt_crate_header.frame_number;
    }

};
}  // end of namespace datatypes
}  // end of namespace uboone
}  // end of namespace fnal
}  // end of namespace gov

// This MACRO must be outside any namespaces.
BOOST_CLASS_VERSION(gov::fnal::uboone::datatypes::crateHeader, gov::fnal::uboone::datatypes::constants::VERSION)    
#endif /* #ifndef BOONETYPES_H */



