#ifndef _UBOONETYPES_CARDDATA_H
#define _UBOONETYPES_CARDDATA_H
#include <memory>
#include <map>
#include <algorithm>
#include <sys/types.h>
#include <inttypes.h>
#include "evttypes.h"

#include <boost/serialization/list.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/binary_object.hpp>

#include "constants.h"
#include "channelData.h"

namespace gov {
namespace fnal {
namespace uboone {
namespace datatypes {

using namespace gov::fnal::uboone;
 
/***
 *  Note: this is the serialization class that handles the card data.
 ***/

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

  void updateIOMode(uint8_t,int);
  uint8_t getIOMode() { return cardData_IO_mode; }

  std::map<int,channelData> getChannelMap() { return channel_map; }
  int getNumberOfChannels() { return channel_map.size(); }
  void insertChannel(int,channelData);

  void decompress();

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


}  // end of namespace datatypes
}  // end of namespace uboone
}  // end of namespace fnal
}  // end of namespace gov

// This MACRO must be outside any namespaces.
BOOST_CLASS_VERSION(gov::fnal::uboone::datatypes::cardData, gov::fnal::uboone::datatypes::constants::VERSION)    

#endif /* #ifndef BOONETYPES_H */



