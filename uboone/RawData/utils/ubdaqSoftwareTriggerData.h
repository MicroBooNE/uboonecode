#ifndef _UBOONETYPES_DAQSWTRIGGERDATA_H
#define _UBOONETYPES_DAQSWTRIGGERDATA_H
#include <sys/types.h>
#include <inttypes.h>
//#include "evttypes.h"
#include <string>
#include <vector>
#include <iostream>
#include "uboone/TriggerSim/UBTriggerTypes.h"

namespace raw{

class ubdaqSoftwareTriggerData {

 public:
  ubdaqSoftwareTriggerData(); // standard constructor

  void addAlgorithm(std::string name, bool pass, uint32_t phmax, uint32_t multiplicity, uint32_t triggerTick, double triggerTime, float prescale);

  int getNumberOfAlgorithms(void);
  std::vector<std::string> getListOfAlgorithms(void);
  
  // pass/veto by name
  bool passedAlgo(std::string algo);
  bool vetoAlgo(std::string algo);
  bool passedAlgos(std::vector<std::string> algos);
  bool vetoAlgos(std::vector<std::string> algos);

  //getters by entry index
  bool getPass(int i);
  uint32_t getPhmax(int i);
  uint32_t getMultiplicity(int i);
  uint32_t getTriggerTick(int i);
  double getTimeSinceTrigger(int i);
  std::string getTriggerAlgorithm(int i);
  float getPrescale(int i);

  //getters by name
  uint32_t getPhmax(std::string algo);  // max adc sum at (software) trigger firing time
  uint32_t getMultiplicity(std::string algo);  // multiplicity at (software) trigger firing time
  uint32_t getTriggerTick(std::string algo);   // tick since the beam-gate opened
  double getTimeSinceTrigger(std::string algo);  // time since the event (hardware) trigger, in us
  int getID(std::string algo); // get the index of a given algorithm
  float getPrescale(std::string algo); // returns prescale_weight (see below)

 private:

  std::vector<std::pair<std::string,bool> > passAlgo; // list of algorithms and corresponding bools for pass/fail
  std::vector<uint32_t> PHMAX; // max adc sum at (software) trigger firing time
  std::vector<uint32_t> multiplicity; // multiplicity at (software) trigger firing time
  std::vector<uint32_t> triggerTick; // tick since the beam-gate opened
  std::vector<double> triggerTime; // time since the event (hardware) trigger, in us
  std::vector<float> prescale_weight; // 1/prescale_weight gives the fraction of events that are let through

};


}  // end of namespace raw

#endif 



