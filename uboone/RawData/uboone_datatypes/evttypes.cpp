#include "evttypes.h"

/* event type variables here */
#ifndef __CXX__                 /* g++ doesn't like multiple definitions. - 10/12/04 - sm */
uint32_t ev_type_toggle[N_EVENT_TYPES_MAX];
uint32_t ev_type_prescale[N_EVENT_TYPES_MAX];
uint32_t ev_type_count[N_EVENT_TYPES_MAX];
#endif

const char * TranslateEventType(uint16_t event_type, const int abbrev) {
    /* if abbrev is set to non-zero, then the abbreviation is used in some cases. Use macro
       TRANSLATE_ABBREV(=1) or TRANSLATE_NOABBREV(=0) to specify.*/

    switch (event_type) {
    case UNUSED_TYPE:
        return "";
        break;
    case BEAM_TYPE:
        return "Beam";
        break;
    case STROBE_TYPE:
        return "Strobe";
        break;
    case CALIB_TYPE:
        return abbrev ? "" : "Calibration";
        break;
    case BNB_TYPE:
        return abbrev ? "BNB" : "BNB Neutrino";
        break;
    case NUMI_TYPE:
        return abbrev ? "NuMI" : "NuMI Neutrino";
        break;
    case SN_TYPE:
        return abbrev ? "SuperN" : "Supernova";
        break;
    case GRANITE_TYPE:
        return abbrev ? "Granite" : "Granite";
        break;
    case BEGIN_TYPE:
        return abbrev ? "Begin" : "Begin Run event";
        break;
    case RESUME_TYPE:
        return abbrev ? "Resume" : "Resume Run event";
        break;
    case PAUSE_TYPE:
        return abbrev ? "Pause" : "Pause Run event";
        break;
    case TPC_ERROR_TYPE:
        return abbrev ? "TPC err" : "TPC error event";
        break;
    case END_TYPE:
        return abbrev ? "End" : "End Run event";
        break;
    case TEST_TYPE:
        return abbrev ? "Test" : "Test event";
        break;
    default:
        return abbrev ? "" : "Illegal Event type";
        break;
    }
}
