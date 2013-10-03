#ifndef _UBOONETYPES_CARDHEADER_H
#define _UBOONETYPES_CARDHEADER_H
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
   Note: this is the serialization class for the hardcoded card_header_t struct.
   If there is a change there, it needs to be made here too. If there is a 
   change made here, it may do no good without a change there.
***/


class cardHeader {
 public:
  static const uint8_t DAQ_version_number = gov::fnal::uboone::datatypes::constants::VERSION;
  cardHeader();
  cardHeader(card_header_t cardH) { bt_card_header = cardH; }

  uint32_t getIDAndModuleWord() { return bt_card_header.id_and_module; }
  uint32_t getWordCountWord() { return bt_card_header.word_count; }
  uint32_t getEventWord() { return bt_card_header.event_number; }
  uint32_t getFrameWord() { return bt_card_header.frame_number; }  
  uint32_t getChecksumWord() { return bt_card_header.checksum; }
  
  void setIDAndModuleWord(uint32_t word) { bt_card_header.id_and_module = word; }
  void setWordCountWord(uint32_t word) { bt_card_header.word_count = word; }
  void setEventNumberWord(uint32_t word) { bt_card_header.event_number = word; }
  void setFrameNumberWord(uint32_t word) { bt_card_header.frame_number = word; }  
  void setChecksumWord(uint32_t word) { bt_card_header.checksum = word; }

  card_header_t getCardHeader() { return bt_card_header; }
  void setCardHeader(card_header_t cardH) { bt_card_header = cardH; }

  uint32_t getID();
  uint32_t getModule();
  uint32_t getEvent();
  uint32_t getFrame();
  uint32_t getChecksum();
  uint32_t getWordCount();
 
  size_t getCardDataSize();

 private:
  card_header_t bt_card_header;

  friend class boost::serialization::access;
  
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      if(version>0)
	ar & bt_card_header.id_and_module
	   & bt_card_header.word_count
	   & bt_card_header.event_number
 	   & bt_card_header.frame_number
	   & bt_card_header.checksum;
    }

};
}  // end of namespace datatypes
}  // end of namespace uboone
}  // end of namespace fnal
}  // end of namespace gov

// This MACRO must be outside any namespaces.
BOOST_CLASS_VERSION(gov::fnal::uboone::datatypes::cardHeader, gov::fnal::uboone::datatypes::constants::VERSION)    
#endif /* #ifndef BOONETYPES_H */



