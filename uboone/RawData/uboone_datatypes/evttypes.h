#ifndef _EVENT_TYPES_H
#define _EVENT_TYPES_H 1
#include <inttypes.h>

/*  global_header.record_type possible bit settings. */
#define RESERVED 1
#define BOONE_SLOW_MON 3
#define ACNET_SLOW_MAGNET 4
#define ACNET_FAST_MULTIWIRE 5
#define LITTLE_MUON_COUNTER 6
#define BEAM_WALL_MONITOR 7
#define ACNET_FAST_IRM 8
#define TRIGGERED 9
#define SUPERNOVA_STREAM 10
#define PMT_DATA 63
#define TPC_DATA 64
#define TPC_ERROR_DATA 65

/*  global_header.record_origin possible values. */
#define DATA 0
#define MC   1


#define N_EVENT_TYPES_MAX 0x40 /* (decimal= 64)  */
#define UNUSED_TYPE 0  /* Don't use this one, since we don't want all-zero bitfield to be valid */
#define BEAM_TYPE  1
#define STROBE_TYPE 2
#define CALIB_TYPE 3
#define SN_TYPE 5
#define TPC_TYPE 6
#define PMT_TYPE 7
#define GRANITE_TYPE 8
#define VETO_TYPE 9               /* For some potential future hardware */ 
#define NO_UPPER_COSMICS_TYPE 10  /* Nevis trigger on inactivity */
#define BNB_TYPE 18
#define NUMI_TYPE 19
#define TPC_ERROR_TYPE 25
#define MRT_TYPE 56
#define TEST_TYPE 57      /* 57 used for diagnostics */
#define BEGIN_TYPE 60     /* 60 */
#define RESUME_TYPE 61    /* 61 */
#define PAUSE_TYPE 62     /* 62 */
#define END_TYPE 63       /* 63 */

/* macro used for string translation of event type values */
#define TRANSLATE_ABBREV 1
#define TRANSLATE_NOABBREV 0

/*Crate header types*/
#define UNKNOWN_HEADER_TYPE 0
#define TPC_HEADER_TYPE 0
#define PMT_HEADER_TYPE 1

/* function declarations */
#ifdef __CXX__
extern "C" char * TranslateEventType(uint16_t event_type, const int abbrev);
#else
char * TranslateEventType(uint16_t event_type, const int abbrev);
#endif

#endif /* #ifdef EVENT_TYPES_H */
