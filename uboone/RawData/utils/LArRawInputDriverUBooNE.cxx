////////////////////////////////////////////////////////////////////////
/// \file  LArRawInputDriverUBooNE.cxx
/// \brief Source to convert raw binary files to root files
/// \Original Authors
/// \version $Id: LArRawInputDriver.h,v 1.7 2010/01/14 19:20:33 brebel Exp $
/// \author  brebel@fnal.gov, soderber@fnal.gov
/// \MicroBooNE Author: jasaadi@fnal.gov (with much help from Wes and Eric)
////////////////////////////////////////////////////////////////////////

#include "RawData/utils/LArRawInputDriverUBooNE.h"

#include "RawData/RawDigit.h"
#include "RawData/DAQHeader.h"
#include "RawData/OpDetPulse.h" //<---I assume I need to include the Optical Detector
#include "Geometry/Geometry.h"
#include "SummaryData/RunData.h"

#include "art/Framework/IO/Sources/put_product_in_principal.h"
#include "art/Utilities/Exception.h"


extern "C" {
#include <sys/types.h>
#include <iostream>
#include <dirent.h>
}

// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ====== PUTTING IN SWIZZLE DEPENDENCIES DIRECTLY TO SEE IF THIS IS THE BEST WAY TO BUILD THIS =========
// ====== Very first steps to "swizzle" is nominally working...lots more work needed on this... =========
// =====                                    not for general use yet (4/7/13)                    =========
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================


using namespace gov::fnal::uboone::datatypes;

// ##################################################################################
// ### Putting in beamData.cpp methods stolen from swizzle_dependencies/datatypes ###
// ##################################################################################
std::ostream & operator<<(std::ostream &os, const gov::fnal::uboone::datatypes::beamData &bd)
{
  os <<"Device name: " << bd.deviceName << std::endl
     <<"Units: "<< bd.units << std::endl
     <<"Value(s): "<<bd.val.at(0);
  for (size_t i=1;i<bd.val.size();i++) os <<", "<<bd.val[i];
  os <<std::endl;

  return os;
}

beamData::beamData()
{
  deviceName="";
  units="";
  val.resize(0);
}

// ####################################################################################
// ### Putting in beamHeader.cpp methods stolen from swizzle_dependencies/datatypes ###
// ####################################################################################
std::ostream & operator<<(std::ostream &os, const gov::fnal::uboone::datatypes::beamHeader &bh)
{
  return os <<"Record type: " << (int)bh.record_type << std::endl
	    <<"Event signal: "<< bh.event_signal << std::endl
	    <<"Seconds: "<< bh.seconds << std::endl
	    <<"Milli seconds: "<< bh.milli_seconds << std::endl
	    <<"Number of devices: "<< (int)bh.number_of_devices << std::endl
	    <<"Number of bytes om: "<< bh.number_of_bytes_in_record << std::endl;	 
}

// ##################################################################################
// ### Putting in cardData.cpp methods stolen from swizzle_dependencies/datatypes ###
// ##################################################################################
beamHeader::beamHeader()
{
  record_type=0;
  event_signal="";
  seconds=0;
  milli_seconds=0;
  number_of_devices=0;
  number_of_bytes_in_record=0;
}


char* cardData::getCardDataPtr(){
  
  if(cardData_IO_mode >= IO_GRANULARITY_CHANNEL){
    std::cout << "ERROR! Granularity is above card level." 
	      << "Cannot return pointer to card data!" << std::endl;
    return nullptr;
  }
  else{
    return card_data_ptr.get();
  }
}

void cardData::setCardDataPtr(char* ptr){

  if(cardData_IO_mode >= IO_GRANULARITY_CHANNEL){
    std::cout << "ERROR! Granularity is above card level." 
	      << "Cannot set pointer to card data!" << std::endl;
  }
  else{
    card_data_ptr.reset(ptr);
  }
}

void cardData::updateIOMode(uint8_t new_mode){

  //we are already at card granularity...so get out if that's the case
  if(new_mode <= IO_GRANULARITY_CARD)
    return;

  if(new_mode >= IO_GRANULARITY_CHANNEL && cardData_IO_mode < IO_GRANULARITY_CHANNEL){

    size_t total_data_read = 0;
    size_t size16 = sizeof(uint16_t);
    while(total_data_read < card_data_size){
      
      //get the channel header word
      std::unique_ptr<uint16_t> channel_header(new uint16_t);
      std::copy(getCardDataPtr() + total_data_read,
		getCardDataPtr() + total_data_read + size16,
		(char*)channel_header.get());
      total_data_read += size16;
      
      //std::cout << "Channel header " << std::hex << *channel_header << std::endl;
      
      //now loop until we find the channel trailer
      std::unique_ptr<uint16_t> channel_trailer(new uint16_t);
      size_t channel_data_size = 0;
      while(total_data_read < card_data_size){
	std::copy(getCardDataPtr() + total_data_read,
		  getCardDataPtr() + total_data_read + size16,
		  (char*)channel_trailer.get());
	//std::cout << "Next Word " << std::hex << *channel_trailer << std::endl;
      	total_data_read += size16;

	if(*channel_trailer & 0x5000)
	  break;
	else 
	  channel_data_size += size16;
      }//end while over channel data

      //pointer to the channel data
      std::shared_ptr<char> channel_data_ptr(new char[channel_data_size]);
      std::copy(getCardDataPtr() + total_data_read - channel_data_size - size16,
		getCardDataPtr() + total_data_read - size16,
		channel_data_ptr.get());

      //now initialise channelData object, and store in map
      channelData chD(channel_data_ptr,channel_data_size,*channel_header,*channel_trailer);
      insertChannel(chD.getChannelNumber(),chD);

    }//end while over card data size

    card_data_ptr.reset();
    cardData_IO_mode = IO_GRANULARITY_CHANNEL;
  }//endif channel granularity update


}

void cardData::insertChannel(int channel_number, channelData chD){
  channel_map.insert(std::pair<int,channelData>(channel_number,chD));
}

// ####################################################################################
// ### Putting in cardHeader.cpp methods stolen from swizzle_dependencies/datatypes ###
// ####################################################################################
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

// ##################################################################################
// ### Putting in crateData.cpp methods stolen from swizzle_dependencies/datatypes ###
// ##################################################################################
char* crateData::getCrateDataPtr(){
  
  if(crateData_IO_mode >= IO_GRANULARITY_CARD){
    std::cout << "ERROR! Granularity is above crate level." 
	      << "Cannot return pointer to crate data!" << std::endl;
    return nullptr;
  }
  else {
    return crate_data_ptr.get();
  }
}

void crateData::setCrateDataPtr(char* ptr){

  if(crateData_IO_mode >= IO_GRANULARITY_CARD){
    std::cout << "ERROR! Granularity is above crate level." 
	      << "Cannot set pointer to crate data!" << std::endl;
  }
  else {
    crate_data_ptr.reset(ptr);
  }
}

void crateData::updateIOMode(uint8_t new_mode){

  //we are already at crate granularity...so get out if that's the case
  if(new_mode <= IO_GRANULARITY_CRATE)
    return;

  if(new_mode >= IO_GRANULARITY_CARD && crateData_IO_mode < IO_GRANULARITY_CARD){

    size_t data_read = 0;
    std::unique_ptr<event_header_t> memblkEH(new event_header_t);
    std::unique_ptr<event_trailer_t> memblkET(new event_trailer_t);

    std::copy(getCrateDataPtr() + data_read,
	      getCrateDataPtr() + data_read + sizeof(event_header_t),
	      (char*)memblkEH.get());
    event_header.setEventHeader(*memblkEH);
    data_read += sizeof(event_header_t);
    
    while(1){
      std::unique_ptr<card_header_t> memblkCardH(new card_header_t);
      std::copy(getCrateDataPtr() + data_read,
		getCrateDataPtr() + data_read + sizeof(card_header_t),
		(char*)memblkCardH.get());
      data_read += sizeof(card_header_t);
      
      cardHeader cardH(*memblkCardH);
      size_t cardDataSize = cardH.getCardDataSize();

      /*
      std::cout << "Card header ...\n"
	<< std::hex << memblkCardH->id_and_module << " " << memblkCardH->word_count << " " 
	<< memblkCardH->event_number << " " << memblkCardH->frame_number<< " " << memblkCardH->checksum << std::dec << std::endl;

      std::cout << std::hex
		<< cardH.getIDAndModuleWord() << " " << cardH.getWordCountWord() << " "
		<< cardH.getEventWord() << " " << cardH.getFrameWord() << " " << cardH.getChecksumWord() << std::endl;
      */

      std::shared_ptr<char> card_data(new char[cardDataSize]);
      std::copy(getCrateDataPtr() + data_read,
		getCrateDataPtr() + data_read + cardDataSize,
		(char*)card_data.get());
      data_read += cardDataSize;
      
      cardData cardD(card_data,cardDataSize);
      if(new_mode == IO_GRANULARITY_CHANNEL) cardD.updateIOMode(new_mode);
      insertCard(cardH,cardD);
      
      std::copy(getCrateDataPtr() + data_read,
		getCrateDataPtr() + data_read + sizeof(event_trailer_t),
		(char*)memblkET.get());
      if(memblkET->trailer == 0xe0000000) break;
    }
    
    event_trailer.setEventTrailer(*memblkET);
    data_read += sizeof(event_trailer_t);
    crate_data_ptr.reset();

    crateData_IO_mode = new_mode;
  } //endif on IO_GRANULARITY_CARD update

  if(new_mode == IO_GRANULARITY_CHANNEL && crateData_IO_mode < IO_GRANULARITY_CHANNEL){
    
    std::map<cardHeader,cardData>::iterator card_it;
    for( card_it = card_map.begin(); card_it != card_map.end(); card_it++)
      (card_it->second).updateIOMode(new_mode);
    
    crateData_IO_mode = new_mode; //eventRecords io_mode

  }//endif on IO_GRANULARITY_CHANNEL update

}

void crateData::insertCard(cardHeader cH, cardData cD){
  card_map.insert(std::pair<cardHeader,cardData>(cH,cD));
}


// #####################################################################################
// ### Putting in crateHeader.cpp methods stolen from swizzle_dependencies/datatypes ###
// #####################################################################################
crateHeader::crateHeader(){

  bt_crate_header.size = 0;
  bt_crate_header.crate_number = 0;
  bt_crate_header.card_count = 0;
  bt_crate_header.event_number = 0;
  bt_crate_header.frame_number = 0;

}


// #####################################################################################
// ### Putting in eventRecord.cpp methods stolen from swizzle_dependencies/datatypes ###
// #####################################################################################
eventRecord::eventRecord() {

  //use default constructors here
  global_header = globalHeader();
  trigger_data = triggerData();
  gps_data = gps();

  beam_header = beamHeader();
  beam_data_vector.clear();
  
  seb_map.clear();

  er_IO_mode = IO_GRANULARITY_CRATE;

}

//this updates all the crates and cards if necessary
void eventRecord::updateIOMode(uint8_t mode) {

  std::map<crateHeader,crateData>::iterator seb_it;
  for( seb_it = seb_map.begin(); seb_it != seb_map.end(); seb_it++)
    (seb_it->second).updateIOMode(mode);

  er_IO_mode = mode; //eventRecords io_mode
}

//insert crateHeader,crateData pair
void eventRecord::insertSEB(crateHeader cH, crateData cD){ 
  seb_map.insert( std::pair<crateHeader,crateData>(cH,cD) ); 
}


// ######################################################################################
// ### Putting in globalHeader.cpp methods stolen from swizzle_dependencies/datatypes ###
// ######################################################################################
globalHeader::globalHeader() {

  //record_type = RESERVED;
  record_origin = 2;
  run_number = 0xdadadada;
  event_number = 0xdadadada;
  event_number_crate = 0xdadadada; 
  
  seconds = 0xdadadada; // GPS clock. Since Jan 1, 2012. 
  // Do we need to worry about Leap seconds?
  milli_seconds = 0xdada;
  micro_seconds = 0xdada;
  nano_seconds = 0xdada;
  numberOfBytesInRecord = 0xdadadada;

  number_of_sebs = 0xd0;

}

// ######################################################################################
// ### Putting in triggerData.cpp methods stolen from swizzle_dependencies/datatypes ###
// ######################################################################################
triggerData::triggerData() {

  bt_trigger_data.trig_event_num = 0;
  bt_trigger_data.trig_event_type = 0;
  bt_trigger_data.frame = 0;
  bt_trigger_data.clock = 0;

}


// ============================================================================
// = \todo: Might need to include evttypes.cpp? Need to think about this some =
// ============================================================================



// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================







namespace {
// ======================================================================
// UBooNE data interface, adapted from code by Rebel/Soderberg:
//  modified J. Asaadi March 20th, 2013


  //Define Structures corresponding to Binary data file.
  /// \todo: First using the structures in boonetypes.h provided from
  /// \	     the swizzle package

  
  // ===================================================================== 
  /// \todo: This is a method stolen from swizzle_dependencies/share/boonetypes.h
  /// \      Clearly need to use swizzling instead...just need to figure out how
  
   struct event_header 	// Event Control Words. These should be at begin and end 
   			// of each crate. Header word format: (should be 0xffffffff)
			// Trailer word format: (should be 0xe0000000) */
   
   {
     uint32_t header;
   
   
   };
  
/**
ADC 16 bit data format: (each daq channel has a first and last word)
bit 15  14  13  12   11 to 0
first             1   1   0   0   channel[5:0]
last              1   1   0   1   adc[11:0]
not compressed    1   0   0   0   adc[11:0]
compressed    0  (code,code,             )
**/

   
   // ======================================================================
  int run( std::string s1 )
  {
    size_t p1 = s1.find("R");
    size_t p2 = s1.find("_E");

    int run = atoi((s1.substr(p1+1,p2-p1-1)).c_str());
    return run;
  }
   
   // ======================================================================
  int event( std::string s1 )
  {
    size_t p1 = s1.find("E");
    size_t p2 = s1.find("_T");

    int event = atoi((s1.substr(p1+1,p2-p1-1)).c_str());
    return event;
  }
   
   // ======================================================================
   bool compare( std::string s1, std::string s2 )
   {
     int r1 = run(s1);
     int r2 = run(s2);
     int e1 = event(s1);
     int e2 = event(s2);

     return r1 == r2 ?  e1 < e2
       :  r1 < r2;
   }
   
   
   // ======================================================================
   std::vector<std::string> getsortedfiles( std::string dir )
   {
     if( dir == "" )
       throw art::Exception( art::errors::Configuration )
         << "Vacuous directory name" << std::endl;

     std::vector<std::string> files;

     DIR * dp = NULL;
     if( (dp = opendir(dir.c_str())) == NULL ) {
       throw art::Exception( art::errors::FileOpenError )
         << "Error opening directory " << dir << std::endl;
     }

     dirent * dirp = NULL;
     while( (dirp = readdir(dp)) != NULL ) {
       std::string filename( dirp->d_name );
       if( filename.find("bin") != std::string::npos ) {
         files.push_back(filename);
       }
     }
     closedir(dp);

     sort( files.begin(), files.end(), compare );

     return files;
  }  // getsortedfiles()

   // =====================================================================
   struct EventFileSentry  //Use RAII (Resource Acquisition Is Initialization)
     {
     explicit EventFileSentry(std::string const &filepath)
     : infile(filepath.c_str(), std::ios_base::in | std::ios_base::binary){ }
     ~EventFileSentry() { infile.close(); }

      std::ifstream infile;
     };
 
 
 
   // ======================================================================
   // (This seems to be where we put DAQ::Header information....so going to swizzle)
   void process_LAr_file(std::string dir,
                        std::string  const &  filename,
                        std::vector<raw::RawDigit>& digitList,
                        raw::DAQHeader& daqHeader)
     {
       // Prepare the input file. The sentry is responsible for making the
       // file stream object, and will *automatically* close it when it
       // goes out of scope *for any reason*, including normal function
       // exit or exception throw.
       
       EventFileSentry efs(dir+"/"+filename);
       std::ifstream &infile =  efs.infile;
    
       // Throwing an exception if the file fails to open
       if( !infile.is_open() ) {
         throw art::Exception( art::errors::FileReadError )
          << "failed to open input file " << filename << std::endl;
      }
      
      // JAsaadi: Just checking that the file opened correctly
      if( infile.is_open()){std::cout<<"YES! The file opened"<<std::endl;}
      
      // ##########################################################################
      // Trying something from the swizzle
      eventRecord event_record;  // Declare an eventRecord object. This is yours.
      
      //std::ifstream ifs("/uboone/app/users/jasaadi/uBoone_DataFormat/xmit-bin-NU-04022013-203237-0.dat", std::ios::binary);
      
      // declare a boost::archive::binary_iarchive object
      boost::archive::binary_iarchive ia(infile); 
      //boost::archive::binary_iarchive ia(ifs); 
      // read in from the archive into your eventRecord object
 
      ia >> event_record;
      
      //Here's how to get info from the global header in the event record.
      //Other header-type objects will be very, very similar.
      globalHeader *global_header = event_record.getGlobalHeaderPtr();
      
      // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      // Reading in the event header information in a similar way to how ArgoNeuT did
      // ...need to think more about this
      
      event_header h1;
      //channel c1;
      //footer f1;
      
      
      //read in header section of file
      infile.read((char *) &h1, sizeof h1);
      
      
      // above here is still using old ArgoNeut Methods...need to fix eventually
      // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
      // ########################################################
      // ### Back to swizzling: using global_header-> methods ###
      // ########################################################
      time_t mytime = global_header->getSeconds();
      
      // I suspect this is antiquated...but might need to confirm?
      
      mytime = mytime << 32;//Nov. 2, 2010 - "time_t" is a 64-bit word on many 64-bit machines
      //so we had to change types in header struct to read in the correct
      //number of bits.  Once we have the 32-bit timestamp from the binary
      //data, shift it up to the upper half of the 64-bit timestamp.  - Mitch
      
      daqHeader.SetStatus(1);
      daqHeader.SetFixedWord(h1.header);
      daqHeader.SetFileFormat(global_header->getRecordType());
      daqHeader.SetSoftwareVersion(global_header->DAQ_version_number);
      daqHeader.SetRun(global_header->getRunNumber());
      daqHeader.SetEvent(global_header->getEventNumber());
      daqHeader.SetTimeStamp(mytime);
      
      /// \todo: What is the "spare word" ? Leaving it unset for now
      //daqHeader.SetSpareWord(h1.spare);
      
      
      
      // ### Swizzling to get the number of channels...trying the method used in write_read.cpp
      // ### provided by Wes --- About the data:
      // ### The format of the data is in levels: crate, card, channel.
      // ### Level 1: The event record contains a map of (crateHeader,crateData) pairs.
      // ### Level 2: Each crateData object may contain a map of (cardHeader,cardData) pairs.
      // ### Level 3: Each cardData object may contain a map of (int,channelData) pairs. The int
      // ### is the channel number.
      
      int channel_number = 0;
      
      //get the seb map, and do a loop over all sebs/crates
      std::map<crateHeader,crateData,compareCrateHeader> seb_map = event_record.getSEBMap();
      std::map<crateHeader,crateData>::iterator seb_it;
      for( seb_it = seb_map.begin(); seb_it != seb_map.end(); seb_it++){
      
      	//get the crateHeader/crateData objects
        crateHeader crate_header = seb_it->first;
        crateData crate_data = seb_it->second;
	
	//now get the card map (for the current crate), and do a loop over all cards
    	std::map<cardHeader,cardData>::iterator card_it;
    	std::map<cardHeader,cardData,compareCardHeader> card_map = crate_data.getCardMap();
    	for(card_it = card_map.begin(); card_it != card_map.end(); card_it++){
	   
	   //get the cardHeader/cardData objects
           cardHeader card_header = card_it->first;
           cardData card_data = card_it->second;
	   
	   //now get the channel map (for the current card), and do a loop over all channels
      	   std::map<int,channelData> channel_map = card_data.getChannelMap();
      	   std::map<int,channelData>::iterator channel_it;
      	   for(channel_it = channel_map.begin(); channel_it != channel_map.end(); channel_it++){
	   
           //get the channel number and channelData
	   int channel_number = channel_it->first;
	   channelData chD = channel_it->second;
      
           //can pull some info from the channelData (channelData objects include some header/trailer words)
      std::cout << "(" << std::dec << channel_number << "," << std::hex << chD.getChannelHeader() << "," << chD.getChannelTrailer() << ") " << std::endl;;
           
      	   }//<--End channel_it for loop
      	}//<---End card_it for loop
      }//<---End seb_it for loop
      
      std::cout<<"channel_number = "<<channel_number<<std::endl;
      
      //daqHeader.SetNChannels(h1.nchan);
      
      
      // THIS IS JUST WRONG...BUT I AM PUTTING SOMETHING HERE FOR NOW
      //vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
/*      unsigned int wiresPerPlane = 240; //FIX ME
      unsigned int planes = 2; //FIX ME
      
      //one digit for every wire on each plane
      digitList.clear();
      digitList.resize(wiresPerPlane*planes);
      //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
      for( int i = 0; i != h1.nchan; ++i ) {
      infile.read((char *) &c1, sizeof c1);
      //Create vector for ADC data, with correct number of samples for this event
      std::vector<short> adclist(c1.samples);
      infile.read((char*)&adclist[0],sizeof(short)*c1.samples);
      // std::cout << "Channel = " << c1.ch ;
      // std::cout << " #Samples = " << c1.samples ;
      // std::cout << " ADC[0] = " << adclist[0] << " ADC[2047] = " << adclist[2047] << std::endl;

      digitList[i] = raw::RawDigit((c1.ch-1), c1.samples, adclist);//subtract one from ch. number...
                                                                   //hence offline channels will always be one lower
                                                                   //than the DAQ480 definition. - mitch 7/8/2009
      digitList[i].SetPedestal(400.); //carl b assures me this will never change. bjr 4/15/2009
      }
      //read in footer section of file...though it's currently empty.
      infile.read((char *) &f1, sizeof f1);*/
      
      
      // ### Consider editing above this
      // ##############################################################################
    }
     
}//<---End namespace


// +++ Blatant stealing from LongBo
namespace lris {
  // ======================================================================
  // class c'tor/d'tor:
  LArRawInputDriverUBooNE::LArRawInputDriverUBooNE(fhicl::ParameterSet const &, // Not used
                                       art::ProductRegistryHelper &helper,
                                       art::PrincipalMaker const &pm)
    :
    principalMaker_(pm)
    , currentDir_        ()
    , inputfiles_        ( )
    , nextfile_          ( inputfiles_.begin() )
    , filesdone_         ( inputfiles_.end() )
    , currentSubRunID_   ( )
  {
    helper.reconstitutes<raw::DAQHeader,              art::InEvent>("daq");
    helper.reconstitutes<std::vector<raw::RawDigit>,  art::InEvent>("daq");
    helper.reconstitutes<sumdata::RunData,            art::InRun>  ("daq");
  }
  
  
  // ======================================================================
  void LArRawInputDriverUBooNE::closeCurrentFile()
  {
    // Nothing to do (See EventFileSentry).
  }
  
  
  // ======================================================================
  void LArRawInputDriverUBooNE::readFile(std::string const &name,
                                   art::FileBlock* &fb)
  {
    // Get the list of event files for this directory.
    currentDir_ = name;
    inputfiles_ = getsortedfiles(currentDir_);
    nextfile_ = inputfiles_.begin();
    filesdone_ = inputfiles_.end();
    currentSubRunID_ = art::SubRunID();

    // Fill and return a new Fileblock.
    fb = new art::FileBlock(art::FileFormatVersion(1, "LArRawInput 2011a"),
                            currentDir_);
  }
  
  // =====================================================================
  bool LArRawInputDriverUBooNE::readNext(art::RunPrincipal* const &inR,
                                   art::SubRunPrincipal* const &inSR,
                                   art::RunPrincipal* &outR,
                                   art::SubRunPrincipal* &outSR,
                                   art::EventPrincipal* &outE)
  {
    if (inputfiles_.empty() || nextfile_ == filesdone_ ) return false;

    // Create empty result, then fill it from current filename:
    std::unique_ptr<std::vector<raw::RawDigit> >  rdcollb ( new std::vector<raw::RawDigit>  );

    raw::DAQHeader daqHeader;
    bool firstEventInRun = (nextfile_ == inputfiles_.begin());

    process_LAr_file( currentDir_, *nextfile_++, *rdcollb, daqHeader );
    std::unique_ptr<raw::DAQHeader>              daqcollb( new raw::DAQHeader(daqHeader) );

    art::RunNumber_t rn = daqHeader.GetRun();
    art::Timestamp tstamp = daqHeader.GetTimeStamp();
    
    if (firstEventInRun)
      {
	std::unique_ptr<sumdata::RunData> rundata(new sumdata::RunData(geo::kBo) );
        currentSubRunID_ = art::SubRunID(rn, 1);
        outR = principalMaker_.makeRunPrincipal(rn, tstamp);
        outSR = principalMaker_.makeSubRunPrincipal(rn,
                                                    currentSubRunID_.subRun(),
                                                    tstamp);
        art::put_product_in_principal(std::move(rundata), *outR, "daq");
      } else if (rn != currentSubRunID_.run())
      {
        throw cet::exception("InconsistentEventStream")
          << "Encountered run #" << rn
          << " while processing events from run #" << currentSubRunID_.run()
          << "\n";
      }

    outE = principalMaker_.makeEventPrincipal(currentSubRunID_.run(),
                                              currentSubRunID_.subRun(),
                                              daqHeader.GetEvent(),
                                              tstamp);

    // Put products in the event.
    art::put_product_in_principal(std::move(rdcollb),
                                  *outE,
                                  "daq"); // Module label
    art::put_product_in_principal(std::move(daqcollb),
                                  *outE,
                                  "daq"); // Module label

    return true;
  }
  
}//<---Endlris

