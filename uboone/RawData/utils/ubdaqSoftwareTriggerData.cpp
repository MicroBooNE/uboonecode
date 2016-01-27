
#include <time.h>
#include "ubdaqSoftwareTriggerData.h"

//-----------------------------------------------------------------------------------

raw::ubdaqSoftwareTriggerData::ubdaqSoftwareTriggerData() {

}

void raw::ubdaqSoftwareTriggerData::addAlgorithm(std::string name_, bool pass_, uint32_t phmax_, uint32_t multiplicity_, uint32_t triggerTick_, double triggerTime_, float prescale_){
  std::pair<std::string, bool> tmp_pair;
  tmp_pair.first = name_;
  tmp_pair.second = pass_;
  
  passAlgo.push_back(tmp_pair);
  PHMAX.push_back(phmax_);
  multiplicity.push_back(multiplicity_);
  triggerTick.push_back(triggerTick_);
  triggerTime.push_back(triggerTime_);
  prescale_weight.push_back(prescale_);
}

//-----------------------------------------------------------------------------------

bool raw::ubdaqSoftwareTriggerData::passedAlgo(std::string algo){
  int id = getID(algo);
  return getPass(id);
//  for (auto const curr_algo: passAlgo){
//    if (curr_algo.first() == algo){ // found the correct algorithm
//      return passAlgo.second();
//    }
//  }
//  std::cout << "WARNING - asked for the status of a trigger algorithm that isn't present" << std::endl;
//  return false;
}

//-----------------------------------------------------------------------------------

bool raw::ubdaqSoftwareTriggerData::vetoAlgo(std::string algo){
  return not passedAlgo(algo);
}

//-----------------------------------------------------------------------------------

bool raw::ubdaqSoftwareTriggerData::passedAlgos(std::vector<std::string> algos){
  for (auto const &algo: algos){
    if (passedAlgo(algo)){
      return true;
    }
  }
  return false;
}

//-----------------------------------------------------------------------------------

bool raw::ubdaqSoftwareTriggerData::vetoAlgos(std::vector<std::string> algos){
  return not passedAlgos(algos);
}

//-----------------------------------------------------------------------------------

uint32_t raw::ubdaqSoftwareTriggerData::getPhmax(std::string algo){  // max adc sum at (software) trigger firing time
  int id = getID(algo);
  return getPhmax(id);
//  for (int i(0); i < PHMAX.size(); ++i){
//    if (passAlgo.at(i).first() == algo){ // correct algorithm
//      return PHMAX.at(i);
//    }
//  }
//  std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!"
//  return 0;
}

//-----------------------------------------------------------------------------------

uint32_t raw::ubdaqSoftwareTriggerData::getMultiplicity(std::string algo){  // multiplicity at (software) trigger firing time
  int id = getID(algo);
  return getMultiplicity(id);
//  for (int i(0); i < multiplicity.size(); ++i){
//    if (passAlgo.at(i).first() == algo){ // correct algorithm
//      return multiplicity.at(i);
//    }
//  }
//  std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!"
//  return 0;

}

//-----------------------------------------------------------------------------------

uint32_t raw::ubdaqSoftwareTriggerData::getTriggerTick(std::string algo){   // tick since the beam-gate opened
  int id = getID(algo);
  return getTriggerTick(id);
//  for (int i(0); i < triggerTick.size(); ++i){
//    if (passAlgo.at(i).first() == algo){ // correct algorithm
//      return triggerTick.at(i);
//    }
//  }
//  std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!"
//  return 0;

}

//-----------------------------------------------------------------------------------

double raw::ubdaqSoftwareTriggerData::getTimeSinceTrigger(std::string algo){  // time since the event (hardware) trigger, in us
  int id = getID(algo);
  return getTimeSinceTrigger(id);
//  for (int i(0); i < passAlgo.size(); ++i){
//    if (passAlgo.at(i).first() == algo){ // correct algorithm
//      return triggerTime.at(i);
//    }
//  }
//  std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!"
//  return -999;
}

//-----------------------------------------------------------------------------------

int raw::ubdaqSoftwareTriggerData::getID(std::string algo){
  for (unsigned int i(0); i < passAlgo.size(); ++i){
    if (passAlgo.at(i).first == algo){ // correct algorithm
      return i;
    }
  }
  std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl;
  return -999;
} 

//-----------------------------------------------------------------------------------

float raw::ubdaqSoftwareTriggerData::getPrescale(std::string algo){
  int id = getID(algo);
  return getPrescale(id);
//  for (int i(0); i < passAlgo.size(); ++i){
//    if (passAlgo.at(i).first() == algo){ // correct algorithm
//      return prescale_weight.at(i);
//    }
//  }
//  std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!"
//  return 1;
}

//-----------------------------------------------------------------------------------

bool raw::ubdaqSoftwareTriggerData::getPass(int i) { 
  if (i > 0){
    if ((unsigned)i > passAlgo.size()){
      std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl;
      return 0;
    }
    return passAlgo.at(i).second;
  }
  std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl; // negative index;
  return 0;
}

//-----------------------------------------------------------------------------------


uint32_t raw::ubdaqSoftwareTriggerData::getPhmax(int i) {
  if (i > 0){
    if ((unsigned)i > PHMAX.size()){
      std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl;
      return 0;
    }
    return PHMAX.at(i);
  }
  std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl; // negative index;
  return 0;

} 

//-----------------------------------------------------------------------------------

uint32_t raw::ubdaqSoftwareTriggerData::getMultiplicity(int i) {
  if (i > 0){
    if ((unsigned)i > multiplicity.size()){
      std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl;
      return 0;
    }
    return multiplicity.at(i);
  }
  std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl; // negative index
  return 0;
} 

//-----------------------------------------------------------------------------------

uint32_t raw::ubdaqSoftwareTriggerData::getTriggerTick(int i) {
  if (i > 0){
    if ((unsigned)i > triggerTick.size()){
      std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl;
      return 0;
    }
    return triggerTick.at(i);
  }
  std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl; // negative index
  return 0;
}  

//-----------------------------------------------------------------------------------

double raw::ubdaqSoftwareTriggerData::getTimeSinceTrigger(int i) {
  if (i > 0){
    if ((unsigned)i > triggerTime.size()){
      std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl;
      return -999;
    }
    return triggerTime.at(i);
  }
  std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl; // negative index
  return -999;
}

//-----------------------------------------------------------------------------------

std::string raw::ubdaqSoftwareTriggerData::getTriggerAlgorithm(int i) {
  if (i > 0){
    if ((unsigned)i > passAlgo.size()){
      std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl;
      return "";
    }
    return passAlgo.at(i).first;
  }
  std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl; // negative index
  return "";
}

//-----------------------------------------------------------------------------------

float raw::ubdaqSoftwareTriggerData::getPrescale(int i){
  if (i > 0){
    if ((unsigned)i > prescale_weight.size()){
      std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl;
      return 1;
    }
    return prescale_weight.at(i);
  }
  std::cout << "WARNING - asked for information on a trigger algorithm that isn't present!" << std::endl; // negative index
  return 1;
}



//-----------------------------------------------------------------------------------

//void raw::ubdaqSoftwareTriggerData::addEmptyEntries(int new_size){ // function adds empty entries to vectors, so we can set values for higher entries than currently exist
//  if (passAlgo.size() > new_size){return;}
//
//  passAlgo.reserve(new_size);
//  PHMAX.reserve(new_size);
//  multiplicity.reserve(new_size);
//  triggerTick.reserve(new_size);
//  triggerTime.reserve(new_size);
//  prescale_weight.reserve(new_size);
//  
//  for (int i(0); i < new_size - passAlgo.size(); ++i){
//    //std::pair<std::string,int> tmp;
//    passAlgo.emplace_back("",0);
//    PHMAX.emplace_back(0);
//    multiplicity.emplace_back(0);
//    triggerTick.emplace_back(0);
//    triggerTime.emplace_back(0);
//    prescale_weight.emplace_back(0);
//  }
//  
//  return;
//}

//-----------------------------------------------------------------------------------
