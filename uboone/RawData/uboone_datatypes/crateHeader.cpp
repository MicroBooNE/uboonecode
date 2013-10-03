#include "crateHeader.h"

using namespace gov::fnal::uboone::datatypes;

crateHeader::crateHeader(){

  bt_crate_header.complete = false;
  bt_crate_header.crateBits = 15;
  bt_crate_header.size = 0;
  bt_crate_header.crate_number = 0;
  bt_crate_header.card_count = 0;
  bt_crate_header.event_number = 0;
  bt_crate_header.frame_number = 0;

}
