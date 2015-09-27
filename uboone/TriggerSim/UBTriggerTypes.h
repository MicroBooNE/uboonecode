
#ifndef UBTRIGGERTYPES_H
#define UBTRIGGERTYPES_H

namespace trigger{

  // Trigger bits for microboone trigger data product
  typedef enum _ubtrigger_t {
    kPMTTriggerBeam   = 0,
    kPMTTriggerCosmic = 1,
    kPMTTrigger   = 7,
    kTriggerPC    = 8,
    kTriggerEXT   = 9,
    kActive       = 10,
    kTriggerBNB   = 11,
    kTriggerNuMI  = 12,
    kVeto         = 13,
    kTriggerCalib = 14,
    kFakeGate     = 17,
    kFakeBeam     = 18,
    kSpare        = 19
  } UBTrigger_t;


}
#endif
