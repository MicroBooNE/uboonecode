#ifndef _UBOONETYPES_DAQSWTRIGGERDATA_H
#define _UBOONETYPES_DAQSWTRIGGERDATA_H
#include <sys/types.h>
#include <inttypes.h>
#include "evttypes.h"
#include <string>
//#include "share/boonetypes.h"
//
//#include <boost/serialization/list.hpp>
//#include <boost/serialization/string.hpp>
//#include <boost/serialization/version.hpp>
//#include "boost/date_time/gregorian/gregorian.hpp"
//#include "boost/date_time/posix_time/posix_time.hpp"
//#include "boost/date_time/gregorian/greg_serialize.hpp"
//#include "boost/date_time/posix_time/time_serialize.hpp"
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

class daqSoftwareTriggerData {

 public:
//  static const uint8_t DAQ_version_number = gov::fnal::uboone::datatypes::constants::VERSION;
  daqSoftwareTriggerData(){};
  //~daqSoftwareTriggerData();
//  daqSoftwareTriggerData(daqSoftwareTriggerData &td) { pass= td.pass;  } //copy constructor?

  bool getPass() { return pass; }
  uint32_t getPhmax() { return PHMAX; }
  uint32_t getMultiplicity() { return multiplicity; }
  uint32_t getTriggerTick() { return triggerTick; }
  std::string getTriggerAlgorithm() { return algorithm; }

  void setPass(bool passFail) {pass = passFail;}
  void setPhmax(uint32_t phmax) {PHMAX = phmax;}
  void setMultiplicity(uint32_t mult) {multiplicity = mult;}
  void setTriggerTick(uint32_t tick) {triggerTick = tick;}
  void setAlgorithm(std::string algo) {algorithm = algo;}

 private:
  uint32_t trigger_event_num;

  bool pass;
  uint32_t PHMAX;
  uint32_t multiplicity;
  uint32_t triggerTick;
  std::string algorithm;
//  friend class boost::serialization::access;
//  
//  template<class Archive>
//    void serialize(Archive & ar, const unsigned int version)
//    {
//      if( version >0 )
//	ar & bt_trigger_data.trig_event_num
//	   & bt_trigger_data.trig_event_type
//	   & bt_trigger_data.frame
//	   & bt_trigger_data.clock;
//    }
};


}  // end of namespace datatypes
}  // end of namespace uboone
}  // end of namespace fnal
}  // end of namespace gov

// This MACRO must be outside any namespaces.
//BOOST_CLASS_VERSION(gov::fnal::uboone::datatypes::daqSoftwareTriggerData, gov::fnal::uboone::datatypes::constants::VERSION)    

#endif /* #ifndef BOONETYPES_H */



