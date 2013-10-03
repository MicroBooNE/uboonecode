#ifndef _UBOONETYPES_DATAHEADERTRAILERPMT_H
#define _UBOONETYPES_DATAHEADERTRAILERPMT_H
#include <sys/types.h>
#include <inttypes.h>
#include "evttypes.h"
#include "share/boonetypes.h"

#include <boost/serialization/list.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>

#include "constants.h"

/**
   Note: these are the serialization classes for the otherwise hardcoded event header and
   trailer structs, located in share/boonetypes. IF changes are made to gps_t, the 
   appropriate changes should be made here as well, and the version number should be 
   increased.
 **/

namespace gov {
namespace fnal {
namespace uboone {
namespace datatypes {

using namespace gov::fnal::uboone;

class dataHeaderPMT {
 
public:
  static const uint8_t DAQ_version_number = gov::fnal::uboone::datatypes::constants::VERSION;
  dataHeaderPMT() {};

  uint16_t getHeader() { return bt_pmt_data_header.header; }
  void setHeader(uint16_t header) { bt_pmt_data_header.header = header; }

  void setDataHeader(pmt_data_header_t DH) { bt_pmt_data_header = DH; }

 private:
  pmt_data_header_t bt_pmt_data_header;
  friend class boost::serialization::access;
  
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      if( version > 0 )
	ar & bt_pmt_data_header.header;
    }    
}; 


class dataTrailerPMT {
 
public:
  static const uint8_t DAQ_version_number = gov::fnal::uboone::datatypes::constants::VERSION;
  dataTrailerPMT() {};

  uint16_t getTrailer() { return bt_pmt_data_trailer.trailer; }
  void setTrailer(uint16_t trailer) { bt_pmt_data_trailer.trailer = trailer; }

  void setDataTrailer(pmt_data_trailer_t DT) { bt_pmt_data_trailer = DT; }

 private:
  pmt_data_trailer_t bt_pmt_data_trailer;
  friend class boost::serialization::access;
  
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      if( version > 0 )
	ar & bt_pmt_data_trailer.trailer;
    }    
}; 


}  // end of namespace sebs
}  // end of namespace uboone
}  // end of namespace fnal
}  // end of namespace gov

// This MACRO must be outside any namespaces.
BOOST_CLASS_VERSION(gov::fnal::uboone::datatypes::dataHeaderPMT, gov::fnal::uboone::datatypes::constants::VERSION)    
BOOST_CLASS_VERSION(gov::fnal::uboone::datatypes::dataTrailerPMT, gov::fnal::uboone::datatypes::constants::VERSION)    
#endif /* #ifndef */



