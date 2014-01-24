#include <boost/serialization/list.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>

#include "beamHeader.h"

using namespace gov::fnal::uboone::datatypes;

std::ostream & operator<<(std::ostream &os, const gov::fnal::uboone::datatypes::beamHeader &bh)
{
  return os <<"Record type: " << (int)bh.getRecordType() << std::endl
	    <<"Event signal: "<< bh.getEventSignal() << std::endl
	    <<"Seconds: "<< bh.getSeconds() << std::endl
	    <<"Milli seconds: "<< bh.getMilliSeconds() << std::endl
	    <<"Number of devices: "<< (int)bh.getNumberOfDevices() << std::endl
	    <<"Number of bytes in record: "<< bh.getNumberOfBytesInRecord() << std::endl;	 
}

beamHeader::beamHeader()
{
  record_type=0;
  event_signal="";
  seconds=0;
  milli_seconds=0;
  number_of_devices=0;
  number_of_bytes_in_record=0;
}
