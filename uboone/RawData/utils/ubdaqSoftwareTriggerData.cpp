
#include <time.h>
#include "ubdaqSoftwareTriggerData.h"

//-----------------------------------------------------------------------------------

raw::ubdaqSoftwareTriggerData::ubdaqSoftwareTriggerData() {

}

void raw::ubdaqSoftwareTriggerData::addAlgorithm(std::string name_, bool pass_, bool pass_prescale_, uint32_t phmax_, uint32_t multiplicity_, uint32_t triggerTick_, double triggerTime_, float prescale_){
  std::pair<std::string, bool> tmp_pair;
  std::pair<std::string, bool> tmp_pair_prescale;
  tmp_pair.first = name_;
  tmp_pair.second = pass_;
  
  passAlgo.push_back(tmp_pair);
  passPrescale.push_back(pass_prescale_);
  PHMAX.push_back(phmax_);
  multiplicity.push_back(multiplicity_);
  triggerTick.push_back(triggerTick_);
  triggerTime.push_back(triggerTime_);
  prescale_weight.push_back(prescale_);
}

//-----------------------------------------------------------------------------------

std::vector<std::string> raw::ubdaqSoftwareTriggerData::getListOfAlgorithms(void) const{
  std::vector<std::string> list;
  for (auto const algoPair: passAlgo){
    list.push_back(algoPair.first);
  }
  return list;
}

//-----------------------------------------------------------------------------------

int raw::ubdaqSoftwareTriggerData::getNumberOfAlgorithms(void) const{
  return passAlgo.size();
}


//-----------------------------------------------------------------------------------

bool raw::ubdaqSoftwareTriggerData::passedAlgo(std::string algo) const{
  int id = getID(algo);
  return getPass(id);
}

//-----------------------------------------------------------------------------------

bool raw::ubdaqSoftwareTriggerData::passedPrescaleAlgo(std::string algo) const{
  int id = getID(algo);
  return getPassPrescale(id);
}

//-----------------------------------------------------------------------------------

bool raw::ubdaqSoftwareTriggerData::vetoAlgo(std::string algo) const{
  return not passedAlgo(algo);
}

//-----------------------------------------------------------------------------------

bool raw::ubdaqSoftwareTriggerData::passedAlgos(std::vector<std::string> algos) const{
  for (auto const &algo: algos){
    if (passedAlgo(algo)){
      return true;
    }
  }
  return (algos.size() == 0);
}


//-----------------------------------------------------------------------------------

bool raw::ubdaqSoftwareTriggerData::passedPrescaleAlgos(std::vector<std::string> algos) const{
  for (auto const &algo: algos){
    if (passedPrescaleAlgo(algo)){
      return true;
    }
  }
  return (algos.size() == 0);
}


//-----------------------------------------------------------------------------------

bool raw::ubdaqSoftwareTriggerData::vetoAlgos(std::vector<std::string> algos) const{
  return not passedAlgos(algos);
}

//-----------------------------------------------------------------------------------

uint32_t raw::ubdaqSoftwareTriggerData::getPhmax(std::string algo) const{  // max adc sum at (software) trigger firing time
  int id = getID(algo);
  return getPhmax(id);
}

//-----------------------------------------------------------------------------------

uint32_t raw::ubdaqSoftwareTriggerData::getMultiplicity(std::string algo) const{  // multiplicity at (software) trigger firing time
  int id = getID(algo);
  return getMultiplicity(id);
}

//-----------------------------------------------------------------------------------

uint32_t raw::ubdaqSoftwareTriggerData::getTriggerTick(std::string algo) const{   // tick since the beam-gate opened
  int id = getID(algo);
  return getTriggerTick(id);
}

//-----------------------------------------------------------------------------------

double raw::ubdaqSoftwareTriggerData::getTimeSinceTrigger(std::string algo) const{  // time since the event (hardware) trigger, in us
  int id = getID(algo);
  return getTimeSinceTrigger(id);
}

//-----------------------------------------------------------------------------------

int raw::ubdaqSoftwareTriggerData::getID(std::string algo) const{
  for (unsigned int i(0); i < passAlgo.size(); ++i){
    if (passAlgo.at(i).first == algo){ // correct algorithm
      return i;
    }
  }
  std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl;
  return -999;
} 

//-----------------------------------------------------------------------------------

float raw::ubdaqSoftwareTriggerData::getPrescale(std::string algo) const{
  int id = getID(algo);
  return getPrescale(id);
}

//-----------------------------------------------------------------------------------

bool raw::ubdaqSoftwareTriggerData::getPass(int i) const { 
  if (i >= 0){
    if ((unsigned)i >= passAlgo.size()){
      std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl;
      return 0;
    }
    return passAlgo.at(i).second;
  }
  std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl; // negative index;
  return 0;
}
//-----------------------------------------------------------------------------------

bool raw::ubdaqSoftwareTriggerData::getPassPrescale(int i) const { 
  if (i >= 0){
    if ((unsigned)i >= passAlgo.size()){
      std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl;
      return 0;
    }
    return passPrescale.at(i);
  }
  std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl; // negative index;
  return 0;
}

//-----------------------------------------------------------------------------------


uint32_t raw::ubdaqSoftwareTriggerData::getPhmax(int i) const{
  if (i >= 0){
    if ((unsigned)i >= PHMAX.size()){
      std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl;
      return 0;
    }
    return PHMAX.at(i);
  }
  std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl; // negative index;
  return 0;

} 

//-----------------------------------------------------------------------------------

uint32_t raw::ubdaqSoftwareTriggerData::getMultiplicity(int i) const{
  if (i >= 0){
    if ((unsigned)i >= multiplicity.size()){
      std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl;
      return 0;
    }
    return multiplicity.at(i);
  }
  std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl; // negative index
  return 0;
} 

//-----------------------------------------------------------------------------------

uint32_t raw::ubdaqSoftwareTriggerData::getTriggerTick(int i) const{
  if (i >= 0){
    if ((unsigned)i >= triggerTick.size()){
      std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl;
      return 0;
    }
    return triggerTick.at(i);
  }
  std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl; // negative index
  return 0;
}  

//-----------------------------------------------------------------------------------

double raw::ubdaqSoftwareTriggerData::getTimeSinceTrigger(int i) const{
  if (i >= 0){
    if ((unsigned)i >= triggerTime.size()){
      std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl;
      return -999;
    }
    return triggerTime.at(i);
  }
  std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl; // negative index
  return -999;
}

//-----------------------------------------------------------------------------------

std::string raw::ubdaqSoftwareTriggerData::getTriggerAlgorithm(int i) const{
  if (i >= 0){
    if ((unsigned)i >= passAlgo.size()){
      std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl;
      return "";
    }
    return passAlgo.at(i).first;
  }
  std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl; // negative index
  return "";
}

//-----------------------------------------------------------------------------------

float raw::ubdaqSoftwareTriggerData::getPrescale(int i) const{
  if (i >= 0){
    if ((unsigned)i >= prescale_weight.size()){
      std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl;
      return 1;
    }
    return prescale_weight.at(i);
  }
  std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl; // negative index
  return 1;
}

