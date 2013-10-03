#ifndef _UBOONETYPES_WINDOWHEADERPMT_H
#define _UBOONETYPES_WINDOWHEADERPMT_H
#include <memory>
#include <algorithm>
#include <sys/types.h>
#include <inttypes.h>
#include "evttypes.h"
#include "share/boonetypes.h"

#include <boost/serialization/list.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>

#include "constants.h"

namespace gov {
namespace fnal {
namespace uboone {
namespace datatypes {

using namespace gov::fnal::uboone;

/***
   Note: this is the serialization class for the hardcoded pmt_window_header_t struct.
   If there is a change there, it needs to be made here too. If there is a 
   change made here, it may do no good without a change there.
***/


class windowHeaderPMT {
 public:
  static const uint8_t DAQ_version_number = gov::fnal::uboone::datatypes::constants::VERSION;
  windowHeaderPMT();
  windowHeaderPMT(pmt_window_header_t wH) { bt_pmt_window_header = wH; }

  uint16_t getChannelAndDiscWord() { return bt_pmt_window_header.ch_and_disc; }
  uint16_t getFrameAndSample1Word() { return bt_pmt_window_header.frame_and_sample1; }
  uint16_t getSample2Word() { return bt_pmt_window_header.sample2; }

  void setChannelAndDiscWord(uint16_t word) { bt_pmt_window_header.ch_and_disc = word; }
  void setFrameAndSample1Word(uint16_t word) { bt_pmt_window_header.frame_and_sample1 = word; }
  void setSample2Word(uint16_t word) { bt_pmt_window_header.sample2 = word; }
  
  uint8_t getChannelNumber();
  uint8_t getDiscriminant();
  uint32_t getSample();
  uint16_t getFrame();

  pmt_window_header_t getWindowHeader() { return bt_pmt_window_header; }
  void setWindowHeader(pmt_window_header_t wH) { bt_pmt_window_header = wH; }
 

 private:
  pmt_window_header_t bt_pmt_window_header;

  friend class boost::serialization::access;
  
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      if(version>0)
	ar & bt_pmt_window_header.ch_and_disc
	   & bt_pmt_window_header.frame_and_sample1
	   & bt_pmt_window_header.sample2;
    }

};
}  // end of namespace datatypes
}  // end of namespace uboone
}  // end of namespace fnal
}  // end of namespace gov

// This MACRO must be outside any namespaces.
BOOST_CLASS_VERSION(gov::fnal::uboone::datatypes::windowHeaderPMT, gov::fnal::uboone::datatypes::constants::VERSION)    
#endif /* #ifndef BOONETYPES_H */



