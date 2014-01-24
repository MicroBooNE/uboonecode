#include "cardHeader.h"

using namespace gov::fnal::uboone::datatypes;

cardHeader::cardHeader(){
  bt_card_header.id_and_module = 0;
  bt_card_header.word_count = 0;
  bt_card_header.event_number = 0;
  bt_card_header.frame_number = 0;
  bt_card_header.checksum = 0;
}

uint32_t cardHeader::getID(){
  uint32_t Crate_ID = (((bt_card_header.id_and_module>>16) & 0xfff)>>5) & 0x7f;  //was called mod_id before
  return Crate_ID;
}

uint32_t cardHeader::getModule(){
  uint32_t Module = (bt_card_header.id_and_module>>16) & 0x1f;  //was called mod_id before
  return Module;
}

uint32_t cardHeader::getEvent(){
  uint32_t Event = ((bt_card_header.event_number>>16) & 0xfff) + ((bt_card_header.event_number& 0xfff) <<12);
  return Event;
}

uint32_t cardHeader::getFrame(){
  uint32_t Frame = ((bt_card_header.frame_number>>16) & 0xfff) + ((bt_card_header.frame_number & 0xfff) <<12);
  return Frame;
}

uint32_t cardHeader::getWordCount(){
  uint32_t wc = (  ((bt_card_header.word_count>>16) & 0xfff)+((bt_card_header.word_count & 0xfff)<<12) );
  return wc;
}

size_t cardHeader::getCardDataSize(){
  size_t DataSize = (getWordCount()+1) * sizeof(uint16_t);
  return DataSize;
}

uint32_t cardHeader::getChecksum(){
  return bt_card_header.checksum;
}
