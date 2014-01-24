#include "globalHeader.h"


using namespace gov::fnal::uboone::datatypes;

globalHeader::globalHeader() {

  record_type        = RESERVED;
  record_origin      = 0xff;
  event_type         = UNUSED_TYPE;
  run_number         = 0xffffffff;
  subrun_number      = 0xffffffff;
  event_number       = 0xffffffff;
  event_number_crate = 0xffffffff; 
  
  // Do we need to worry about Leap seconds?

  seconds       = 0xffffffff; // GPS clock. Since Jan 1, 2012. 
  milli_seconds = 0xffff;
  micro_seconds = 0xffff;
  nano_seconds  = 0xffff;

  numberOfBytesInRecord = 0;
  number_of_sebs = 0;

}

