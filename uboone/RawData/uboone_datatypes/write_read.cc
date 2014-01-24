#include <cstddef> // NULL
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/time.h>
#include "eventRecord.h"
#include "crateData.h"
#include "channelData.h"
#include "share/boonetypes.h"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_member.hpp>

/***********************************************************************************
 
 Overview:

 This is an example of how to read from an existing event record, and also how to 
 write out a modification of that event record. Full documentation will be written 
 up shortly -- this is mainly meant to provide examples of how to do things.

 The event record is written using the boost::serialization classes. It's not 
 necessary to be an expert or even to know that much about them other than they 
 provide a way to read in and write out data in a standardized form that allows the 
 user to recover an existing event structure from a file, and/or write an event 
 structure to the file. A versioning system also exists in these classes, so upates 
 to the event data structure may be made with no impact on reading old structures.

 Most importantly, it is not necessary for end users to understand the full event 
 strucutre. In fact, end users need only worry about where the information they need
 is contained and how to access it. The data types contained in the event record are
 all defined in the 'datatypes/' directory.
 
 End users may:
  (1) Add member functions to existing objects in the event record, like an 
      accessor, should they need them.
  (2) Create additional classes that take data from the event record and format it 
      in a way best suited for a specific purpose. These could be serialized classes, 
      or not.

 End users should not do DAQ things. So please don't:
  (1) Modify datatype members of the serialized classes in the event record. You may
      have more data you'd like permanently saved in the event record. If you do, 
      talk to your friend in DAQ about it.
  (2) Modify the "serialize()", "save()", or "load()" member functions of serialized
      classes. Big no-no. Talk to your friend in DAQ.
  (3) Change the data in officially written event by reading in an event record, 
      modifying it, and then overwriting the event record. If you want to do this, 
      talk to your friend in DAQ.

 ***********************************************************************************/


using namespace gov::fnal::uboone::datatypes;

int main() {

  /*****************************************************************************
   --- STEP 1: Read an eventRecord object from a file.
   This is the easy part: you can read in an entire event with basically one 
   line of code! After this, you have all the data in the event stored in your
   eventRecord object. Now it's just a matter of doing what you'd like with it.
  ******************************************************************************/

  eventRecord event_record;  // Declare an eventRecord object. This is yous.
  std::ifstream ifs("input_event.dat", std::ios::binary); // open the file that has your data (note: binary)
  boost::archive::binary_iarchive ia(ifs); // declare a boost::archive::binary_iarchive object
  ia >> event_record;  // read in from the archive into your eventRecord object
  ifs.close();


  /*****************************************************************************
   --- STEP 2: Update the granularity on your eventRecord object.
   First, a little bit about the data. When we write events in the assembler, we
   will likely write them at a very course granularity. The data from each crate
   arrives at the assembler as a big block of binary data. End users likely need
   something more than that: data broken down card-by-card or even 
   channel-by-channel. The "eventRecord::updateIOMode()" function does this. 
   Granularity modes are definied in datatypes/constants.h.

   It is possible to write events in higher granularity---so, if the event you 
   have is already at or above the mode you want, updating will do nothing. Note,
   if you want a mode coarser than what is there, the code cannot handle that 
   right now.
  ******************************************************************************/

  std::cout << "Current readout mode is " << (int)event_record.getIOMode();
  event_record.updateIOMode(IO_GRANULARITY_CHANNEL);


  /*****************************************************************************
   --- STEP 3: Access your event.
   Now that you have everything you want in the format you want it, go forth, 
   and process your data.
  ******************************************************************************/


  //For example, here's how to get info from the global header in the event record.
  //Other header-type objects will be very, very similar.
  globalHeader *global_header = event_record.getGlobalHeaderPtr();
  std::cout << "\n\nFrom global header, run number is " << std::hex << global_header->getRunNumber() << std::endl;
  std::cout << "From global header, event number is " << std::hex << global_header->getEventNumber() << std::endl;
  std::cout << "From global header, event size is " << std::dec << (double)global_header->getNumberOfBytesInRecord()/1024/1024 << " MB" << std::endl;
  std::cout << "From global header, seconds is " << std::hex << global_header->getSeconds() << std::endl;
  std::cout << "From global header, number of SEBs is " << std::hex << "0x" << (unsigned int)global_header->getNumberOfSEBs() << std::endl;

  /*************************************************************************************
   --- About the data:
   The format of the data is in levels: crate, card, channel.
   Level 1: The event record contains a map of (crateHeader,crateData) pairs.
   Level 2: Each crateData object may contain a map of (cardHeader,cardData) pairs.
   Level 3: Each cardData object may contain a map of (int,channelData) pairs. The int
            is the channel number.
  **************************************************************************************/

  //get the seb map, and do a loop over all sebs/crates
  std::map<crateHeader,crateData,compareCrateHeader> seb_map = event_record.getSEBMap();
  std::map<crateHeader,crateData>::iterator seb_it;
  for( seb_it = seb_map.begin(); seb_it != seb_map.end(); seb_it++){

    //get the crateHeader/crateData objects
    crateHeader crate_header = seb_it->first;
    crateData crate_data = seb_it->second;

    //can check some things in the crate header
    std::cout << "\nFrom crate header, crate (number,event,frame) is ... (" << std::dec
	      << (unsigned int)crate_header.getCrateNumber() << ","
	      << crate_header.getCrateEventNumber() << ","
	      << crate_header.getCrateFrameNumber() << ")" <<std::endl;

    std::cout << "From crate data, data size is " << std::dec << (double)crate_data.getCrateDataSize()/1024/1024 << " MB" << std::endl;
    /*
    std::copy(crate_data.getCrateDataPtr(),
	      crate_data.getCrateDataPtr()+sizeof(event_header_t),
	      (char*)memblkEH.get());
    std::copy(crate_data.getCrateDataPtr()+sizeof(event_header_t),
	      crate_data.getCrateDataPtr()+sizeof(event_header_t)+sizeof(card_header_t),
	      (char*)memblkCardH.get());
    
    std::cout << "Event Header " << std::hex << memblkEH->header << std::endl;
    std::cout << "First Card Header "
      << std::hex << memblkCardH->id_and_module << " " << memblkCardH->word_count << " " 
      << memblkCardH->event_number << " " << memblkCardH->frame_number<< " " << memblkCardH->checksum << std::dec << std::endl;
    */

    //now get the card map (for the current crate), and do a loop over all cards
    std::map<cardHeader,cardData>::iterator card_it;
    std::map<cardHeader,cardData,compareCardHeader> card_map = crate_data.getCardMap();
    for(card_it = card_map.begin(); card_it != card_map.end(); card_it++){

      //get the cardHeader/cardData objects
      cardHeader card_header = card_it->first;
      cardData card_data = card_it->second;

      //can check some things in the card header
      std::cout << "From CARD header, card (IDandModuleWord,WordCountWord,EventWord,FrameWord,ChecksumWord) is ... ( " << std::hex
		<< card_header.getIDAndModuleWord() << " , " << card_header.getWordCountWord() << " , "
		<< card_header.getEventWord() << " , " << card_header.getFrameWord() << " , " << card_header.getChecksumWord() << " )" << std::endl;
      std::cout << "From CARD header, card (ID,Module,WordCount,Event,Frame) is ... ( " << std::dec
		<< card_header.getID() << " , " << card_header.getModule()  << " , " << card_header.getWordCount()
		<< " , " << card_header.getEvent() << " , " << card_header.getFrame() << " )" << std::endl;

      //now get the channel map (for the current card), and do a loop over all channels
      std::map<int,channelData> channel_map = card_data.getChannelMap();
      std::map<int,channelData>::iterator channel_it;
      for(channel_it = channel_map.begin(); channel_it != channel_map.end(); channel_it++){

	//get the channel number and channelData
	int ch_num = channel_it->first;
	channelData chD = channel_it->second;
	
	//can pull some info from the channelData (channelData objects include some header/trailer words)
	std::cout << "(" << std::dec << ch_num << "," << std::hex << chD.getChannelHeader() << "," << chD.getChannelTrailer() << ") " << std::endl;;

	//and for kicks, let's print the first 10 words in each channel
	std::unique_ptr<uint16_t> blk_chD(new uint16_t);
	for (int j=0; j<10; j++){

	  std::copy(chD.getChannelDataPtr()+j*sizeof(uint16_t),
		    chD.getChannelDataPtr()+(j+1)*sizeof(uint16_t),
		    (char*)blk_chD.get());
	  if(j==0) std::cout << "\t";
	  std::cout << std::hex << *blk_chD << " ";
	} //end loop over channel data
	std::cout << std::endl;
      } //end loop over channels

      std::cout << std::endl;
      
    } //end loop over cards/modules

  } // end loop over sebs/crates

  
  /*****************************************************************************
   --- STEP 4 (optional): Write out your event record
   If you want to save your event record here for later use (especially if you 
   modified something in it, and want to keep it locally), it's easy to do so.
  ******************************************************************************/

  std::ofstream ofs("test.dat", std::ios::binary); //open output file
  boost::archive::binary_oarchive oa(ofs); // create output archive object
  oa << event_record; // write the event_record to that archive
  ofs.close();

  return 0;
}
