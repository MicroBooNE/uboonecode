#ifndef _UBOONETYPES_CHANNELDATA_H
#define _UBOONETYPES_CHANNELDATA_H
#include <memory>
#include <algorithm>
#include <sys/types.h>
#include <inttypes.h>
#include "evttypes.h"

#include <boost/serialization/list.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/binary_object.hpp>

#include "constants.h"

namespace gov {
namespace fnal {
namespace uboone {
namespace datatypes {

using namespace gov::fnal::uboone;
 
/***
 *  Note: this is the serialization class that handles the card data.
 ***/

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

  void decompress();

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


}  // end of namespace datatypes
}  // end of namespace uboone
}  // end of namespace fnal
}  // end of namespace gov

// This MACRO must be outside any namespaces.
BOOST_CLASS_VERSION(gov::fnal::uboone::datatypes::channelData, gov::fnal::uboone::datatypes::constants::VERSION)    

#endif /* #ifndef BOONETYPES_H */



