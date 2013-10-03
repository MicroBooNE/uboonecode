#ifndef _UBOONETYPES_BEAMDATA_H
#define _UBOONETYPES_BEAMDATA_H

#include <sys/types.h>

#include <iostream>
#include <vector>

#include <boost/serialization/list.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>

#include "constants.h"

namespace gov {
namespace fnal {
namespace uboone {
namespace datatypes {

using namespace gov::fnal::uboone;

class beamData {

 public:

  beamData();   
  
  std::string getDeviceName() const {return device_name;};
  std::string getUnits() const {return units;};
  uint32_t getSeconds() const {return seconds;}; // GPS clock. Since Jan 1, 2012. 
  uint16_t getMilliSeconds() const {return milli_seconds;};
  const std::vector<double>& getData() const {return device_data;};

  void setDeviceName(std::string val) {device_name=val;};
  void setUnits(std::string val) {units=val;};
  void setSeconds(uint32_t val) {seconds=val;};
  void setMilliSeconds(uint16_t val) {milli_seconds=val;};
  void setData(std::vector<double> val) {device_data=val;};
  void pushData(double val) {device_data.push_back(val);};
 
 private:

  std::string device_name;
  std::string units;
  uint32_t seconds; // GPS clock. Since Jan 1, 2012. 
  uint16_t milli_seconds;
  std::vector<double> device_data;   
  
  friend class boost::serialization::access;
  
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      if(version > 0)
	ar & device_name & units & seconds & milli_seconds & device_data;
    }
  
};
}  // end of namespace datatypes
}  // end of namespace uboone
}  // end of namespace fnal
}  // end of namespace gov

std::ostream & operator<<(std::ostream &os, const gov::fnal::uboone::datatypes::beamData &bd);

// This MACRO must be outside any namespaces.

BOOST_CLASS_VERSION(gov::fnal::uboone::datatypes::beamData, gov::fnal::uboone::datatypes::constants::VERSION)    

#endif /* #ifndef BOONETYPES_H */
