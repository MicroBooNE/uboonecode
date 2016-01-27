#ifndef _UBOONETYPES_DAQSWTRIGGERDATA_H
#define _UBOONETYPES_DAQSWTRIGGERDATA_H
#include <sys/types.h>
#include <inttypes.h>
//#include "evttypes.h"
#include <string>
#include <vector>
#include <iostream>
#include "uboone/TriggerSim/UBTriggerTypes.h"
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

class ubdaqSoftwareTriggerData {

 public:
  ubdaqSoftwareTriggerData(); // standard constructor

  void addAlgorithm(std::string name, bool pass, uint32_t phmax, uint32_t multiplicity, uint32_t triggerTick, double triggerTime, float prescale);
  
  // pass/veto by name
  bool passedAlgo(std::string algo);
  bool vetoAlgo(std::string algo);
  bool passedAlgos(std::vector<std::string> algos);
  bool vetoAlgos(std::vector<std::string> algos);

  //getters by entry number
  bool getPass(int i);// { return pass.at(i); }
  uint32_t getPhmax(int i);// { return PHMAX.at(i); } // max adc sum at (software) trigger firing time
  uint32_t getMultiplicity(int i);// { return multiplicity.at(i); } // multiplicity at (software) trigger firing time
  uint32_t getTriggerTick(int i);// { return triggerTick.at(i); }  // tick since the beam-gate opened
  double getTimeSinceTrigger(int i);// { return triggerTime.at(i); } // time since the event (hardware) trigger, in us
  std::string getTriggerAlgorithm(int i);// { return algorithm.at(i); }
  float getPrescale(int i);//{return prescale_weight.at(i);}

  //getters by name
  uint32_t getPhmax(std::string algo);  // max adc sum at (software) trigger firing time
  uint32_t getMultiplicity(std::string algo);  // multiplicity at (software) trigger firing time
  uint32_t getTriggerTick(std::string algo);   // tick since the beam-gate opened
  double getTimeSinceTrigger(std::string algo);  // time since the event (hardware) trigger, in us
  int getID(std::string algo); 
  float getPrescale(std::string algo);
  


 private:

  std::vector<std::pair<std::string,bool>> passAlgo;
  std::vector<uint32_t> PHMAX; // max adc sum at (software) trigger firing time
  std::vector<uint32_t> multiplicity; // multiplicity at (software) trigger firing time
  std::vector<uint32_t> triggerTick; // tick since the beam-gate opened
  std::vector<double> triggerTime; // time since the event (hardware) trigger, in us
  std::vector<float> prescale_weight; // 1/prescale_weight gives the fraction of events that are let through

//  void addEmptyEntries(int i); // not sure I need this
  
//  // setters - kept private so all users can do is add a whole entry
//  void setPass(int i, bool passFail);// { addEmptyEntries(i); pass.at(i) = passFail;}
//  void setPhmax(int i,uint32_t phmax);// {addEmptyEntries(i); PHMAX.at(i) = phmax;}
//  void setMultiplicity(int i, uint32_t mult);// {addEmptyEntries(i); multiplicity.at(i) = mult;}
//  void setTriggerTick(int i, uint32_t tick);// {addEmptyEntries(i); triggerTick.at(i) = tick;} // Does not affect trigger time number - needs to be set independently
//  void setTimeSinceTrigger(int i, double time);// {addEmptyEntries(i);  triggerTime.at(i) = time; } // Does not affect trigger tick number - needs to be set independently
//  void setAlgorithm(int i, std::string algo);// {addEmptyEntries(i); algorithm.at(i) = algo;}
//  
////  void setTriggerBit(::trigger::_ub_daqsw_trigger_t trig, bool pass);
//  void setPrescale(int i, float weight);//{addEmptyEntries(i); prescale_weight.at(i) = weight;}

};


}  // end of namespace raw

#endif 



