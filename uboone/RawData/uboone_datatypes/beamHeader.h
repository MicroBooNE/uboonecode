#ifndef _UBOONETYPES_BEAMHEADER_H
#define _UBOONETYPES_BEAMHEADER_H

#include <sys/types.h>

#include <boost/serialization/list.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>

#include "constants.h"

#include <iostream>

namespace gov {
namespace fnal {
namespace uboone {
namespace datatypes {

using namespace gov::fnal::uboone;

class beamHeader {

 public:
   beamHeader();

   uint8_t getRecordType() const {return record_type;};
   std::string getEventSignal() const {return event_signal;};
   uint32_t getSeconds() const {return seconds;}; // GPS clock. Since Jan 1, 2012. 
   uint16_t getMilliSeconds() const {return milli_seconds;};
   uint16_t getNumberOfDevices() const {return number_of_devices;};
   uint32_t getNumberOfBytesInRecord() const {return number_of_bytes_in_record;};

   void setRecordType(uint8_t val) {record_type=val;};
   void setEventSignal(std::string val) {event_signal=val;};
   void setSeconds(uint32_t val) {seconds=val;};
   void setMilliSeconds(uint16_t val) {milli_seconds=val;};
   void setNumberOfDevices(uint16_t val) {number_of_devices=val;};
   void setNumberOfBytesInRecord(uint32_t val) {number_of_bytes_in_record=val;};

   bool operator<( const beamHeader & h ) const {
     return ((this->seconds < h.seconds) || (this->seconds == h.seconds && this->milli_seconds<h.milli_seconds))  ;
   }
   bool operator<=( const beamHeader & h ) const {
     return ((this->seconds < h.seconds) || (this->seconds == h.seconds && this->milli_seconds<=h.milli_seconds))  ;
   }

 private:

   uint8_t record_type;
   std::string event_signal;
   uint32_t seconds; // GPS clock. Since Jan 1, 2012. 
   uint16_t milli_seconds;
   uint16_t number_of_devices;
   uint32_t number_of_bytes_in_record;

   friend class boost::serialization::access;

   template<class Archive>
     void serialize(Archive & ar, const unsigned int version)
     {

       if(version > 0)
	 ar & record_type & event_signal & seconds & milli_seconds & number_of_devices & number_of_bytes_in_record;

     }
  
};
}  // end of namespace datatypes
}  // end of namespace uboone
}  // end of namespace fnal
}  // end of namespace gov

std::ostream & operator<<(std::ostream &os, const gov::fnal::uboone::datatypes::beamHeader &bh);

// This MACRO must be outside any namespaces.

BOOST_CLASS_VERSION(gov::fnal::uboone::datatypes::beamHeader, gov::fnal::uboone::datatypes::constants::VERSION)    


#endif /* #ifndef BOONETYPES_H */
