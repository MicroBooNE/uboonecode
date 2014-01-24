#ifndef BOONETYPES_H
#define BOONETYPES_H
#include <sys/types.h>
#include <inttypes.h>
#include "event_types.h"
#include "../gps/trigBoardClock.h"

/**
   Note: these are hardcoded structs that are used by the sebs for processing.
   They do not have explicit versioning attached to them, and thus care should 
   be made in handling them. There are serialized versions of these structs used
   by the assembler, located in datatypes/. If you make a change to the structs 
   in this file, you need to make appropriate changes there.

   AGAIN: IF YOU MAKE A CHANGE TO THE STRUCTS IN THIS FILE, YOU NEED TO MAKE THE 
   APPROPRIATE CHANGES IN THE CORRESPONDING CLASS IN THE 'DATATYPES' FOLDER. If 
   you do not, we will very likely have a bug.
 **/


/*
  EC, Sep, 2012. Event i comes from sebs with a crate_header, then
  N * {card_header + M*(channel_header+channel_data)} structures.
  N is # of cards/crate. M is # of channels per card.
  Crate 10 gives us additionally a trigger_header read over its devoted fiber.

  Assembler will tack on a global header and a trigger_header.
  The trigger_header will come fully formed from crate 10.
 */


/* The gps structure is what is read off the card
 * in the trigger crate and represents the GPS time at the arrival of a PPS.
 */

typedef struct gps {
    uint32_t lower;
    uint32_t upper;
} gps_t;  /* 8 bytes */


/* the trigger_data structure includes everything the
 * trigger is going to send to assembler. Will come faithfully
 * from trigger card. Hence, below is only a guess. Need to
 * get this from Chi/Bill.
 */
#define N_ACTIVITY_HIST 4
typedef struct trigger_data {
  uint32_t     trig_event_num;  /* trigger_event_number */
  uint16_t     trig_event_type; /* trigger event type e.g. beam, calib */
  uint16_t frame;  /* frame # where trigger happened*/
  uint64_t clock;  /* Master Crate Clock value where trigger happened*/
} trigger_data_t;  /* xyz bytes */


/** EVENT Control Words. These should be at begin and end of each crate.
Header word format: (should be 0xffffffff)
Trailer word format: (should be 0xe0000000)
**/
typedef struct event_header
{
  uint32_t header;
}event_header_t;

typedef struct event_trailer
{
  uint32_t trailer;
}event_trailer_t;


/**Card/Module Header format: (each module sends a header followed by adc data)
first word  {16'hffff,               4'h7,id[6:0],module[4:0]}
2nd word    {4'h7,eventnumber[23:12], 4'h7,eventnumber[11:0]}
3d word     {4'h7,wordcount[23:12],  4'h7,wordcount[11:0]}
4th word    {4'h7,checksum[23:12],   4'h7,checksum[11:0]}
**/
typedef struct card_header 
{
  uint32_t id_and_module; 
  uint32_t word_count;  //this is number of 16-bit words.
  uint32_t event_number;  
  uint32_t frame_number;  
  uint32_t checksum;
}card_header_t;


/** CHANNEL Control Words.
    First 4 bits are 0100 (0101) for first (last) word of channel.
    Remaining 12 bits are the channel #.
    Use this for degugging on seb side.
**/
typedef struct channel_header
{
    uint16_t channel_begin; /* OR'd channel|"first"*/
} channel_header_t;

typedef struct channel_trailer
{
    uint16_t channel_end; /* OR'd channel|"first"*/
} channel_trailer_t;


/**
ADC 16 bit data format: (each daq channel has a first and last word)
bit 15  14  13  12   11 to 0
first         1   1   0   0   channel[5:0]
last          1   1   0   1   adc[11:0]
not compressed    1   0   0   0   adc[11:0]
compressed    0  (code,code,             )
**/



/**
   This replaces tmub ../gps/symm.h
 **/
typedef struct gps_time 
{
  // 2^32 = 4.E9 . Thus 32 bits allows for (2013-1970)*3.14e7 seconds and 
  // enough nanoseconds to span a second.
  uint32_t second; // seconds since the epoch. 
  uint32_t nano;  // Nanoseconds since the second. 
  gps_time(){};
  //  gps_time(&gps_time_t){}; // copy constructor
} gps_time_t;



/**
   This header is created at seb, with crate-level information to be sent to assembler.
 **/
typedef struct crate_header
{
  bool complete; //bit for knowing if sub-event is guaranteed complete
  uint16_t crateBits; // 4 bits for crate 0 through 9, and one for PMT/TPC
  uint32_t size; //bytes, needs to be uint32_t for large events
  uint8_t crate_number; // Crate #
  uint8_t card_count; // Card count
  uint32_t event_number; // Event #
  uint32_t frame_number; // Frame #
  gps_time_t gps_time; // Inserted for SEB-10 only in rawFragmentDMASource.cpp: PPS time
  tbclkub_t daqClock_time; // Inserted for SEB-10 only in rawFragmentDMASource.cpp: PPS frame/sample/div
  crate_header() {};
} crate_header_t;



/**PMT Card/Module Header format: (each module sends a header followed by data)**/
typedef struct pmt_card_header 
{
  uint32_t id_and_module; 
  uint32_t word_count;  //this is number of 16-bit words.
  uint32_t event_number;  
  uint32_t frame_number;  
  uint32_t checksum;
  uint32_t trig_frame_and_sample;  
} pmt_card_header_t;


/** PMT Control Words. These should be at begin and end of each block of pmt data.
Header word format: (should be 0x4000)
Trailer word format: (should be 0xc000)
**/
typedef struct pmt_data_header
{
  uint16_t header;
}pmt_data_header_t;

typedef struct pmt_data_trailer
{
  uint16_t trailer;
}pmt_data_trailer_t;

typedef struct pmt_window_header 
{
  uint16_t ch_and_disc; 
  uint16_t frame_and_sample1; //bits 6-8: last 3 readout frame bits; bits 1-5: upper 5 bits of readout sample
  uint16_t sample2; //lower 12 bits of readout sample
} pmt_window_header_t;

#endif /* #ifndef BOONETYPES_H */



