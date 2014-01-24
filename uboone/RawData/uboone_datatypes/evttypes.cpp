#include "evttypes.h"

/* event type variables here */
#ifndef __CXX__                 /* g++ doesn't like multiple definitions. - 10/12/04 - sm */
uint32_t ev_type_toggle[N_EVENT_TYPES_MAX];
uint32_t ev_type_prescale[N_EVENT_TYPES_MAX];
uint32_t ev_type_count[N_EVENT_TYPES_MAX];
#endif

char * TranslateEventType(uint16_t event_type, const int abbrev) {
    /* if abbrev is set to non-zero, then the abbreviation is used in some cases. Use macro
       TRANSLATE_ABBREV(=1) or TRANSLATE_NOABBREV(=0) to specify.*/

    switch (event_type) {
    case UNUSED_TYPE:
        return((char*)"");
        break;
    case BEAM_TYPE:
        return((char*)"Beam");
        break;
    case STROBE_TYPE:
        return((char*)"Strobe");
        break;
    case CALIB_TYPE:
        return(abbrev ? (char*)"" : (char*)"Calibration");
        break;
    case BNB_TYPE:
        return(abbrev ? (char*)"BNB" : (char*)"BNB Neutrino");
        break;
    case NUMI_TYPE:
        return(abbrev ? (char*)"NuMI" : (char*)"NuMI Neutrino");
        break;
    case SN_TYPE:
        return(abbrev ? (char*)"SuperN" : (char*)"Supernova");
        break;
    case GRANITE_TYPE:
        return(abbrev ? (char*)"Granite" : (char*)"Granite");
        break;
    case BEGIN_TYPE:
        return(abbrev ? (char*)"Begin" : (char*)"Begin Run event");
        break;
    case RESUME_TYPE:
        return(abbrev ? (char*)"Resume" : (char*)"Resume Run event");
        break;
    case PAUSE_TYPE:
        return(abbrev ? (char*)"Pause" : (char*)"Pause Run event");
        break;
    case TPC_ERROR_TYPE:
        return(abbrev ? (char*)"TPC err" : (char*)"TPC error event");
        break;
    case END_TYPE:
        return(abbrev ? (char*)"End" : (char*)"End Run event");
        break;
    case TEST_TYPE:
        return(abbrev ? (char*)"Test" : (char*)"Test event");
        break;
    default:
        return(abbrev ? (char*)"" : (char*)"Illegal Event type");
        break;
    }
}
