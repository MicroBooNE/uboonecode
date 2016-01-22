#ifndef _UBOONETYPES_DAQSWTRIGGERDATA_H
#define _UBOONETYPES_DAQSWTRIGGERDATA_H
#include <sys/types.h>
#include <inttypes.h>
//#include "evttypes.h"
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
//#include "constants.h"

namespace raw{

//using namespace gov::fnal::uboone;

/**
   Note: this is the serialization class for the otherwise hardcoded trigger_data_t struct,
   located in share/boonetypes. IF changes are made to trigger_data_t, the appropriate changes
   should be made here as well, and the version number should be increased.
 **/

class ubdaqSoftwareTriggerData {

 public:
  ubdaqSoftwareTriggerData(); // standard constructor

  bool getPass() { return pass; }
  uint32_t getPhmax() { return PHMAX; } // max adc sum at (software) trigger firing time
  uint32_t getMultiplicity() { return multiplicity; } // multiplicity at (software) trigger firing time
  uint32_t getTriggerTick() { return triggerTick; }  // tick since the beam-gate opened
  double getTimeSinceTrigger() { return triggerTime; } // time since the event (hardware) trigger, in us
  std::string getTriggerAlgorithm() { return algorithm; }

  void setPass(bool passFail) {pass = passFail;}
  void setPhmax(uint32_t phmax) {PHMAX = phmax;}
  void setMultiplicity(uint32_t mult) {multiplicity = mult;}
  void setTriggerTick(uint32_t tick) {triggerTick = tick;} // Does not affect trigger time number - needs to be set independently
  void setTimeSinceTrigger(double time) { triggerTime = time; } // Does not affect trigger tick number - needs to be set independently
  void setAlgorithm(std::string algo) {algorithm = algo;}

 private:

  bool pass;
  uint32_t PHMAX; // max adc sum at (software) trigger firing time
  uint32_t multiplicity; // multiplicity at (software) trigger firing time
  uint32_t triggerTick; // tick since the beam-gate opened
  double triggerTime; // time since the event (hardware) trigger, in us
  std::string algorithm;

};


}  // end of namespace raw

// This MACRO must be outside any namespaces.
//BOOST_CLASS_VERSION(gov::fnal::uboone::datatypes::daqSoftwareTriggerData, gov::fnal::uboone::datatypes::constants::VERSION)    

#endif /* #ifndef BOONETYPES_H */



