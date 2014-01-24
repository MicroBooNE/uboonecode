#ifndef _UBOONETYPES_CRATEDATA_H
#define _UBOONETYPES_CRATEDATA_H
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
#include "share/boonetypes.h"
#include "eventHeaderTrailer.h"
#include "cardHeader.h"
#include "cardData.h"

namespace gov {
namespace fnal {
namespace uboone {
namespace datatypes {

using namespace gov::fnal::uboone;
 
/***
 *  Note: this is the serialization class that handles the data.
 ***/

struct compareCardHeader {
  bool operator() ( cardHeader lhs, cardHeader rhs) const
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
  
  void decompress();

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


}  // end of namespace datatypes
}  // end of namespace uboone
}  // end of namespace fnal
}  // end of namespace gov

// This MACRO must be outside any namespaces.
BOOST_CLASS_VERSION(gov::fnal::uboone::datatypes::crateData, gov::fnal::uboone::datatypes::constants::VERSION)    

#endif /* #ifndef BOONETYPES_H */



