////////////////////////////////////////////////////////////////////////
/// \file  LArRawInputDriverUBooNE.h
/// \brief Source to convert raw binary files to root files for MicroBooNE
///
/// \Original Version from:
/// \version $Id: T962ConvertBinaryToROOT.h,v 1.7 2010/01/14 19:20:33 brebel Exp $
/// \author  brebel@fnal.gov, soderber@fnal.gov
/// \MicroBooNE author: jasaadi@fnal.gov (with much help from Eric and Wes)
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ProductRegistryHelper.h"
#include "art/Framework/Core/PrincipalMaker.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "art/Framework/Principal/SubRunPrincipal.h"
#include "art/Framework/Principal/EventPrincipal.h"
#include "art/Persistency/Provenance/SubRunID.h"


#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/binary_object.hpp>


#include "RawData/utils/boonetypes.h"
#include <fstream>
#include <string>
#include <vector>

///Conversion of binary data to root files
namespace lris {
  class LArRawInputDriverUBooNE;
}


//================================================================================================================================
//================================================================================================================================
//================================================================================================================================
//================================================================================================================================
namespace gov {
namespace fnal {
namespace uboone {
namespace datatypes {

using namespace gov::fnal::uboone;



enum {
    IO_GRANULARITY_CRATE,
    IO_GRANULARITY_CARD,
    IO_GRANULARITY_CHANNEL
  };
  
namespace constants
{
  const int VERSION = 1; // A dB query eventually.
    // ... other related constants

} // namespace constants



class channelData {

 public:
  static const uint8_t DAQ_version_number = gov::fnal::uboone::datatypes::constants::VERSION;
  
  channelData()
    { channel_data_ptr.reset(); channel_data_size=0; }

  channelData(std::shared_ptr<char> data_ptr, size_t size, uint16_t header, uint16_t trailer)
    { channel_data_ptr.swap(data_ptr); channel_data_size=size; 
      channel_data_header = header; channel_data_trailer = trailer; }

  char* getChannelDataPtr() { return channel_data_ptr.get(); }
  void setChannelDataPtr(char* ptr) {channel_data_ptr.reset(ptr);}

  size_t getChannelDataSize() {return channel_data_size;}
  void setChannelDataSize(size_t size) { channel_data_size = size; }

  void setChannelHeader(uint16_t header) { channel_data_header = header; }
  void setChannelTrailer(uint16_t trailer) { channel_data_trailer = trailer; }

  uint16_t getChannelHeader() { return channel_data_header; }
  uint16_t getChannelTrailer() { return channel_data_trailer; }

  int getChannelNumber() 
  { int number = channel_data_header & 0x3fff; return number; }

 private:
  uint16_t channel_data_header;
  std::shared_ptr<char> channel_data_ptr;
  uint16_t channel_data_trailer;

  size_t channel_data_size;

  friend class boost::serialization::access;
  
  /***
      Use different save and load techniques here so that on the load, 
      we first read the data size, and then we declare space large 
      enough to hold it. After that we copy the data into our buffer, 
      and then swap the pointer to that buffer with out channelData member.
   ***/

  template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
      if(version>0) {
	ar & channel_data_size;
	ar & channel_data_header;
	ar & boost::serialization::make_binary_object(channel_data_ptr.get(),channel_data_size);
	ar & channel_data_trailer;
      }
    }

  template<class Archive>
    void load(Archive & ar, const unsigned int version) 
    {
      if(version>0) { 
	ar & channel_data_size;
	ar & channel_data_header;

	std::shared_ptr<char> data_ptr(new char[channel_data_size]);
	ar & boost::serialization::make_binary_object(data_ptr.get(),channel_data_size);
	channel_data_ptr.swap(data_ptr);

	ar & channel_data_trailer;
      }
    }
  BOOST_SERIALIZATION_SPLIT_MEMBER()
   
};

class beamData {

 public:
   std::string deviceName;
   std::string units;
   std::vector<double> val;
   
   beamData();   
   
 private:

   friend class boost::serialization::access;

   template<class Archive>
     void serialize(Archive & ar, const unsigned int version)
     {
       if(version > 0)
         ar & deviceName & units & val;
     }
   
};


class beamHeader {

 public:
   uint8_t record_type;
   std::string event_signal;
   uint32_t seconds; // GPS clock. Since Jan 1, 2012. 
   uint16_t milli_seconds;
   uint8_t  number_of_devices;
   uint32_t number_of_bytes_in_record;

   beamHeader();

 private:

   friend class boost::serialization::access;

   template<class Archive>
     void serialize(Archive & ar, const unsigned int version)
     {

       if(version > 0)
	 ar & record_type & event_signal & seconds & milli_seconds & number_of_devices & number_of_bytes_in_record;

     }
  
};



class cardData {

 public:
  static const uint8_t DAQ_version_number = gov::fnal::uboone::datatypes::constants::VERSION;
  
  cardData()
    { card_data_ptr.reset(); card_data_size=0; cardData_IO_mode = IO_GRANULARITY_CARD;}

  cardData(std::shared_ptr<char> data_ptr, size_t size)
    { card_data_ptr.swap(data_ptr); card_data_size=size; cardData_IO_mode = IO_GRANULARITY_CARD;}

  char* getCardDataPtr();
  void setCardDataPtr(char*);

  size_t getCardDataSize() {return card_data_size;}
  void setCardDataSize(size_t size) { card_data_size = size; }

  void updateIOMode(uint8_t);
  uint8_t getIOMode() { return cardData_IO_mode; }

  std::map<int,channelData> getChannelMap() { return channel_map; }
  void insertChannel(int,channelData);

 private:
  std::shared_ptr<char> card_data_ptr;
  size_t card_data_size;

  uint8_t cardData_IO_mode;

  std::map<int,channelData> channel_map;

  friend class boost::serialization::access;
  
  /***
      Use different save and load techniques here so that on the load, 
      we first read the data size, and then we declare space large 
      enough to hold it. After that we copy the data into our buffer, 
      and then swap the pointer to that buffer with out cardData member.
   ***/

  template<class Archive> void save(Archive & ar, const unsigned int version) const
    {
      if(version>0) {
	ar & card_data_size;
	ar & cardData_IO_mode;

	if(cardData_IO_mode==IO_GRANULARITY_CARD)
	  ar & boost::serialization::make_binary_object(card_data_ptr.get(),card_data_size);

	else if(cardData_IO_mode >=IO_GRANULARITY_CHANNEL)
	  ar & channel_map;
      }
    }

  template<class Archive> void load(Archive & ar, const unsigned int version) 
    {
      if(version>0) { 
	ar & card_data_size;
	ar & cardData_IO_mode;
	
	if(cardData_IO_mode==IO_GRANULARITY_CARD){
	  std::shared_ptr<char> data_ptr(new char[card_data_size]);
	  ar & boost::serialization::make_binary_object(data_ptr.get(),card_data_size);
	  card_data_ptr.swap(data_ptr);
	}
	else if(cardData_IO_mode>=IO_GRANULARITY_CHANNEL)
	  ar & channel_map;

      }//endif version 0
    }
  BOOST_SERIALIZATION_SPLIT_MEMBER()
   
};





class cardHeader {
 public:
  static const uint8_t DAQ_version_number = gov::fnal::uboone::datatypes::constants::VERSION;
  cardHeader();
  cardHeader(card_header_t cardH) { bt_card_header = cardH; }

  uint32_t getIDAndModuleWord() { return bt_card_header.id_and_module; }
  uint32_t getWordCountWord() { return bt_card_header.word_count; }
  uint32_t getEventWord() { return bt_card_header.event_number; }
  uint32_t getFrameWord() { return bt_card_header.frame_number; }  
  uint32_t getChecksumWord() { return bt_card_header.checksum; }
  
  void setIDAndModuleWord(uint32_t word) { bt_card_header.id_and_module = word; }
  void setWordCountWord(uint32_t word) { bt_card_header.word_count = word; }
  void setEventNumberWord(uint32_t word) { bt_card_header.event_number = word; }
  void setFrameNumberWord(uint32_t word) { bt_card_header.frame_number = word; }  
  void setChecksumWord(uint32_t word) { bt_card_header.checksum = word; }

  card_header_t getCardHeader() { return bt_card_header; }
  void setCardHeader(card_header_t cardH) { bt_card_header = cardH; }

  uint32_t getID();
  uint32_t getModule();
  uint32_t getEvent();
  uint32_t getFrame();
  uint32_t getChecksum();
  uint32_t getWordCount();
 
  size_t getCardDataSize();

 private:
  card_header_t bt_card_header;

  friend class boost::serialization::access;
  
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      if(version>0)
	ar & bt_card_header.id_and_module
	   & bt_card_header.word_count
	   & bt_card_header.event_number
 	   & bt_card_header.frame_number
	   & bt_card_header.checksum;
    }

};
class eventHeader {
 
public:
  static const uint8_t DAQ_version_number = gov::fnal::uboone::datatypes::constants::VERSION;
  eventHeader() {};

  uint32_t getHeader() { return bt_event_header.header; }
  void setHeader(uint32_t header) { bt_event_header.header = header; }

  void setEventHeader(event_header_t EH) { bt_event_header = EH; }

 private:
  event_header_t bt_event_header;
  friend class boost::serialization::access;
  
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      if( version > 0 )
	ar & bt_event_header.header;
    }    
}; 


class eventTrailer {
 
public:
  static const uint8_t DAQ_version_number = gov::fnal::uboone::datatypes::constants::VERSION;
  eventTrailer() {};

  uint32_t getTrailer() { return bt_event_trailer.trailer; }
  void setTrailer(uint32_t trailer) { bt_event_trailer.trailer = trailer; }

  void setEventTrailer(event_trailer_t ET) { bt_event_trailer = ET; }

 private:
  event_trailer_t bt_event_trailer;
  friend class boost::serialization::access;
  
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      if( version > 0 )
	ar & bt_event_trailer.trailer;
    }    
}; 



struct compareCardHeader {
  bool operator() ( cardHeader lhs, cardHeader rhs)
  { return lhs.getModule() < rhs.getModule(); }
};

class crateData {

 public:
  static const uint8_t DAQ_version_number = gov::fnal::uboone::datatypes::constants::VERSION;
  
  crateData()
    { crate_data_ptr.reset(); crate_data_size=0; crateData_IO_mode = IO_GRANULARITY_CRATE;}

  crateData(std::shared_ptr<char> data_ptr, size_t size)
    { crate_data_ptr.swap(data_ptr); crate_data_size=size; crateData_IO_mode = IO_GRANULARITY_CRATE; }

  size_t getCrateDataSize() {return crate_data_size;}
  void setCrateDataSize(size_t size) { crate_data_size = size; }

  char* getCrateDataPtr();// { return crate_data_ptr.get(); }
  void setCrateDataPtr(char*);// {crate_data_ptr.reset(ptr);}

  void updateIOMode(uint8_t);
  uint8_t getIOMode() { return crateData_IO_mode; }

  void insertCard(cardHeader,cardData);

  std::map<cardHeader,cardData,compareCardHeader> getCardMap() { return card_map;}

 private:
  uint8_t crateData_IO_mode;
  
  std::shared_ptr<char> crate_data_ptr;
  size_t crate_data_size;
  
  eventHeader event_header;
  std::map<cardHeader,cardData,compareCardHeader> card_map;
  eventTrailer event_trailer;

  friend class boost::serialization::access;
  
  /***
      Use different save and load techniques here so that on the load, 
      we first read the data size, and then we declare space large 
      enough to hold it. After that we copy the data into our buffer, 
      and then swap the pointer to that buffer with out crateData member.
   ***/

  template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
      
      if(version>0) {
	ar & crate_data_size;
	ar & crateData_IO_mode;
	
	if(crateData_IO_mode == IO_GRANULARITY_CRATE){
	  ar & boost::serialization::make_binary_object(crate_data_ptr.get(),crate_data_size);
	}
	else if(crateData_IO_mode >= IO_GRANULARITY_CARD){
	  ar & event_header;
	  ar & card_map;
	  ar & event_trailer;
	}
	
      }//endif version
      
    }
  
  template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
      
      if(version>0) { 
	ar & crate_data_size;
	ar & crateData_IO_mode;
	
	if(crateData_IO_mode==IO_GRANULARITY_CRATE){
	  std::shared_ptr<char> data_ptr(new char[crate_data_size]);
	  ar & boost::serialization::make_binary_object(data_ptr.get(),crate_data_size);
	  crate_data_ptr.swap(data_ptr);
	}
	else if(crateData_IO_mode>=IO_GRANULARITY_CARD){
	  ar & event_header;
	  ar & card_map;
	  ar & event_trailer;
	}
   
      }
      
    }
  
  BOOST_SERIALIZATION_SPLIT_MEMBER()
    
};


class crateHeader {
 public:
  static const uint8_t DAQ_version_number = gov::fnal::uboone::datatypes::constants::VERSION;
  crateHeader();  
  crateHeader(crate_header_t cH) { bt_crate_header = cH; }

  uint32_t getCrateSize() { return bt_crate_header.size; }
  uint8_t getCrateNumber() { return bt_crate_header.crate_number; }
  uint8_t getCardCount() { return bt_crate_header.card_count; }
  uint32_t getCrateEventNumber() { return bt_crate_header.event_number; }
  uint32_t getCrateFrameNumber() { return bt_crate_header.frame_number; }  

  void setCrateSize(uint32_t size) { bt_crate_header.size = size; }
  void setCrateNumber(uint8_t cnum) { bt_crate_header.crate_number = cnum; }
  void setCardCount(uint8_t ccnt) { bt_crate_header.card_count = ccnt; }
  void setCrateEventNumber(uint32_t event) { bt_crate_header.event_number = event; }
  void setCrateFrameNumber(uint32_t frame) { bt_crate_header.frame_number = frame; }

  crate_header_t getCrateHeader() { return bt_crate_header; }
  void setCrateHeader(crate_header_t cH) { bt_crate_header = cH; }

 private:
  crate_header_t bt_crate_header;

  friend class boost::serialization::access;
  
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      if(version>0)
	ar & bt_crate_header.size
	   & bt_crate_header.crate_number
	   & bt_crate_header.card_count
	   & bt_crate_header.event_number
 	   & bt_crate_header.frame_number;
    }

};



class globalHeader {

 public:
  static const uint8_t DAQ_version_number = gov::fnal::uboone::datatypes::constants::VERSION;
  globalHeader();
  
  void setRecordType(uint8_t type) { record_type = type; }
  void setRecordOrigin(uint8_t origin) { record_origin = origin; }
  void setRunNumber(uint32_t run) { run_number = run; }
  void setEventNumber(uint32_t event) { event_number = event; }
  void setEventNumberCrate(uint32_t event) { event_number_crate = event; }

  void setSeconds(uint32_t s) { seconds = s; }
  void setMilliSeconds(uint16_t ms) { milli_seconds = ms; }
  void setMicroSeconds(uint16_t us) { micro_seconds = us; }
  void setNanoSeconds(uint16_t ns) { nano_seconds = ns; }
  void setNumberOfBytesInRecord(uint32_t size) { numberOfBytesInRecord = size; }

  uint8_t getRecordType() { return record_type; }
  uint8_t getRecordOrigin() { return record_origin; }
  uint32_t getRunNumber() { return run_number; }
  uint32_t getEventNumber() { return event_number; }
  uint32_t getEventNumberCrate() { return event_number_crate; }
  
  uint32_t getSeconds() { return seconds; }
  uint16_t getMilliSeconds() { return milli_seconds; }
  uint16_t getMicroSeconds() { return micro_seconds; }
  uint16_t getNanoSeconds() { return nano_seconds; }
  uint32_t getNumberOfBytesInRecord() { return numberOfBytesInRecord; }

  uint8_t getNumberOfSEBs() { return number_of_sebs;}
  void setNumberOfSEBs(uint8_t s) { number_of_sebs = s;}

 private:
  uint8_t record_type;   /* From event_types.h */
  uint8_t record_origin; /* DATA or MC */
  uint32_t run_number;
  uint32_t event_number;
  uint32_t event_number_crate; /* Crate's sense of the evt #. */
  
  uint32_t seconds; // GPS clock. Since Jan 1, 2012. 
  uint16_t milli_seconds;
  uint16_t micro_seconds;
  uint16_t nano_seconds;
  uint32_t numberOfBytesInRecord;

  uint8_t number_of_sebs; //put this in just to test versioning

  friend class boost::serialization::access;
  
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      if(version>0)
	ar & record_type & record_origin
 	   & run_number & event_number & event_number_crate
	   & seconds & milli_seconds & micro_seconds & nano_seconds
	   & numberOfBytesInRecord & number_of_sebs;

    }
  
};


class gps {
 
public:
  static const uint8_t DAQ_version_number = gov::fnal::uboone::datatypes::constants::VERSION;
  gps() {};
  gps(uint32_t _lower, uint32_t _upper){ bt_gps.lower=_lower; bt_gps.upper=_upper; }

  uint32_t getLower() { return bt_gps.lower; }
  uint32_t getUpper() { return bt_gps.upper; }

  void setLower(uint32_t lower) { bt_gps.lower = lower; }
  void setUpper(uint32_t upper) { bt_gps.upper = upper; }

 private:
  gps_t bt_gps;
  friend class boost::serialization::access;
  
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      if( version > 0 )
	ar & bt_gps.lower & bt_gps.upper;
    }    
}; 


class triggerData {

 public:
  static const uint8_t DAQ_version_number = gov::fnal::uboone::datatypes::constants::VERSION;
  triggerData();
  triggerData(trigger_data_t bt) { bt_trigger_data = bt; }

  uint32_t getTrigEventNum() { return bt_trigger_data.trig_event_num; }
  uint16_t getTrigEventType() { return bt_trigger_data.trig_event_type; }
  uint16_t getFrame() { return bt_trigger_data.frame; }
  uint64_t getClock() { return bt_trigger_data.clock; }

  void setTrigEventNum(uint32_t event) {bt_trigger_data.trig_event_num = event;}
  void setTrigEventType(uint16_t type) {bt_trigger_data.trig_event_type = type;}
  void setFrame(uint16_t frame) {bt_trigger_data.frame = frame;}
  void setClock(uint64_t clock) {bt_trigger_data.clock = clock;}

 private:
  trigger_data_t bt_trigger_data;
  friend class boost::serialization::access;
  
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      if( version >0 )
	ar & bt_trigger_data.trig_event_num
	   & bt_trigger_data.trig_event_type
	   & bt_trigger_data.frame
	   & bt_trigger_data.clock;
    }
};

//used for map
struct compareCrateHeader {
  bool operator() ( crateHeader lhs, crateHeader rhs)
  { return lhs.getCrateNumber() < rhs.getCrateNumber(); }
};

class eventRecord {

 public:
  static const uint8_t DAQ_version_number = gov::fnal::uboone::datatypes::constants::VERSION;
  eventRecord();  
  void setGlobalHeader (globalHeader gH) { global_header = gH; }
  void setTriggerData (triggerData tD) { trigger_data = tD; }
  void setGPS (gps g) { gps_data = g; }
  void setBeamHeader (beamHeader bH) { beam_header = bH; }

  void insertBeamData (beamData bD) { beam_data_vector.push_back(bD); }
  void insertSEB(crateHeader,crateData); //in .cpp file
  
  globalHeader getGlobalHeader() { return global_header; }
  triggerData getTriggerData() { return trigger_data; }
  gps getGPS() { return gps_data; }
  beamHeader getBeamHeader() { return beam_header; }

  std::vector<beamData> getBeamDataVector() { return beam_data_vector; }
  std::map<crateHeader,crateData,compareCrateHeader> getSEBMap() { return seb_map; }

  globalHeader* getGlobalHeaderPtr() { return &global_header; }
  triggerData* getTriggerDataPtr() { return &trigger_data; }
  gps* getGPSPtr() { return &gps_data; }
  beamHeader* getBeamHeaderPtr() { return &beam_header; }

  int getSEBMap_size() { return seb_map.size(); }
  void clearSEBMap() { seb_map.clear(); }

  int getBeamDataVecotr_size() { return beam_data_vector.size(); }
  void clearBeamDataVector() { beam_data_vector.clear(); }

  uint8_t getIOMode() { return er_IO_mode; }
  void updateIOMode(uint8_t); //in .cpp file

 private:

  globalHeader global_header;
  triggerData trigger_data;
  gps gps_data;
  std::map<crateHeader,crateData,compareCrateHeader> seb_map;
  beamHeader beam_header;
  std::vector<beamData> beam_data_vector;
  
  uint8_t er_IO_mode;

  friend class boost::serialization::access;
  
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      if(version>0)
	ar & er_IO_mode
	   & global_header
	   & trigger_data
	   & gps_data
	   & beam_header & beam_data_vector //beam stuff...empty at first, added in later
	   & seb_map;
    }

};

}
}
}
}

//================================================================================================================================
//================================================================================================================================
//================================================================================================================================
//================================================================================================================================
//================================================================================================================================
//================================================================================================================================



class lris::LArRawInputDriverUBooNE {
  /// Class to fill the constraints on a template argument to the class,
  /// FileReaderSource
 public:
  // Required constructor
  LArRawInputDriverUBooNE(fhicl::ParameterSet const &pset,
                    art::ProductRegistryHelper &helper,
                    art::PrincipalMaker const &pm);

  // Required by FileReaderSource:
  void closeCurrentFile();
  void readFile(std::string const &name,
                art::FileBlock* &fb);
  bool readNext(art::RunPrincipal* const &inR,
                art::SubRunPrincipal* const &inSR,
                art::RunPrincipal* &outR,
                art::SubRunPrincipal* &outSR,
                art::EventPrincipal* &outE);

 private:
  // --- data members:
  typedef  std::vector<std::string>  stringvec_t;

  art::PrincipalMaker            principalMaker_;
  std::string                    currentDir_;
  stringvec_t                    inputfiles_;
  stringvec_t::const_iterator    nextfile_;
  stringvec_t::const_iterator    filesdone_;
  art::SubRunID                  currentSubRunID_; 
};  // LArRawInputDriverUBooNE
