#ifndef _UBOONETYPES_TRIGGERDATA_H
#define _UBOONETYPES_TRIGGERDATA_H
#include <sys/types.h>
#include <inttypes.h>
#include "evttypes.h"
#include "share/boonetypes.h"

#include <boost/serialization/list.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/date_time/gregorian/greg_serialize.hpp"
#include "boost/date_time/posix_time/time_serialize.hpp"
#include "constants.h"

namespace gov {
namespace fnal {
namespace uboone {
namespace datatypes {

using namespace gov::fnal::uboone;

/**
   Note: this is the serialization class for the otherwise hardcoded trigger_data_t struct,
   located in share/boonetypes. IF changes are made to trigger_data_t, the appropriate changes
   should be made here as well, and the version number should be increased.
 **/

class triggerData {

 public:
  static const uint8_t DAQ_version_number = gov::fnal::uboone::datatypes::constants::VERSION;
  triggerData();
  triggerData(trigger_data_t bt) { bt_trigger_data = bt; }

  uint32_t getTrigEventNum() { return bt_trigger_data.trig_event_num; }
  uint16_t getTrigEventType() { return bt_trigger_data.trig_event_type; }
  uint16_t getFrame() { return bt_trigger_data.frame; }
  uint64_t getClock() { return bt_trigger_data.clock; }

  void setTrigEventNum(uint32_t event) {bt_trigger_data.trig_event_num = event;}
  void setTrigEventType(uint16_t type) {bt_trigger_data.trig_event_type = type;}
  void setFrame(uint16_t frame) {bt_trigger_data.frame = frame;}
  void setClock(uint64_t clock) {bt_trigger_data.clock = clock;}

 private:
  trigger_data_t bt_trigger_data;
  friend class boost::serialization::access;
  
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      if( version >0 )
	ar & bt_trigger_data.trig_event_num
	   & bt_trigger_data.trig_event_type
	   & bt_trigger_data.frame
	   & bt_trigger_data.clock;
    }
};


}  // end of namespace datatypes
}  // end of namespace uboone
}  // end of namespace fnal
}  // end of namespace gov

// This MACRO must be outside any namespaces.
BOOST_CLASS_VERSION(gov::fnal::uboone::datatypes::triggerData, gov::fnal::uboone::datatypes::constants::VERSION)    

#endif /* #ifndef BOONETYPES_H */



