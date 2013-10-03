#include "cardHeaderPMT.h"

using namespace gov::fnal::uboone::datatypes;

cardHeaderPMT::cardHeaderPMT(){
  bt_pmt_card_header.id_and_module = 0;
  bt_pmt_card_header.word_count = 0;
  bt_pmt_card_header.event_number = 0;
  bt_pmt_card_header.frame_number = 0;
  bt_pmt_card_header.checksum = 0;
  bt_pmt_card_header.trig_frame_and_sample = 0;
}

uint32_t cardHeaderPMT::getID(){
  uint32_t Crate_ID = (((bt_pmt_card_header.id_and_module>>16) & 0xfff)>>5) & 0x7f;  //was called mod_id before
  return Crate_ID;
}

uint32_t cardHeaderPMT::getModule(){
  uint32_t Module = (bt_pmt_card_header.id_and_module>>16) & 0x1f;  //was called mod_id before
  return Module;
}

uint32_t cardHeaderPMT::getEvent(){
  uint32_t Event = ((bt_pmt_card_header.event_number>>16) & 0xfff) + ((bt_pmt_card_header.event_number& 0xfff) <<12);
  return Event;
}

uint32_t cardHeaderPMT::getFrame(){
  uint32_t Frame = ((bt_pmt_card_header.frame_number>>16) & 0xfff) + ((bt_pmt_card_header.frame_number & 0xfff) <<12);
  return Frame;
}

uint32_t cardHeaderPMT::getWordCount(){
  uint32_t wc = (  ((bt_pmt_card_header.word_count>>16) & 0xfff)+((bt_pmt_card_header.word_count & 0xfff)<<12) );
  return wc;
}

size_t cardHeaderPMT::getCardDataSize(){
  size_t DataSize = (getWordCount()+1) * sizeof(uint16_t);
  return DataSize;
}

uint32_t cardHeaderPMT::getChecksum(){
  return bt_pmt_card_header.checksum;
}

uint32_t cardHeaderPMT::getTrigFrame(){
  uint32_t Trig_Frame = (((bt_pmt_card_header.trig_frame_and_sample>>16) & 0xfff)>>4 & 0xf) + (((getFrame())>>4)<<4);
  return Trig_Frame;
}

uint32_t cardHeaderPMT::getTrigSample(){
  uint32_t Trig_Sample = (((bt_pmt_card_header.trig_frame_and_sample>>16) & 0xf)<<8) + (bt_pmt_card_header.trig_frame_and_sample & 0xff);
  return Trig_Sample;
}

