
#include <time.h>
#include "triggerData.h"

using namespace gov::fnal::uboone::datatypes;

triggerData::triggerData() {

  bt_trigger_data.trig_event_num = 0;
  bt_trigger_data.trig_event_type = 0;
  bt_trigger_data.frame = 0;
  bt_trigger_data.clock = 0;

}


