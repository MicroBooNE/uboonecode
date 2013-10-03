#ifndef _UBOONETYPES_WINDOWDATAPMT_H
#define _UBOONETYPES_WINDOWDATAPMT_H
#include <memory>
#include <algorithm>
#include <vector>
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

class windowDataPMT {

 public:
  static const uint8_t DAQ_version_number = gov::fnal::uboone::datatypes::constants::VERSION;
  
  windowDataPMT()
    { window_data_ptr.reset(); window_data_size=0;}

  ~windowDataPMT()
    { window_data_ptr.reset(); }

  windowDataPMT(size_t wd_size, char* wd_ptr)
    { window_data_size = wd_size;
      std::shared_ptr<char> data_ptr(new char[window_data_size]);
      std::copy(wd_ptr,wd_ptr+wd_size,data_ptr.get());
      window_data_ptr.swap(data_ptr); }

  char* getWindowDataPtr() { return window_data_ptr.get(); }
  void setWindowDataPtr(char* ptr) {window_data_ptr.reset(ptr);}

  size_t getWindowDataSize() {return window_data_size;}
  void setWindowDataSize(size_t size) { window_data_size = size; }

 private:
  std::shared_ptr<char> window_data_ptr;
  size_t window_data_size;

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
	ar & window_data_size;
	ar & boost::serialization::make_binary_object(window_data_ptr.get(),window_data_size);
      }
    }

  template<class Archive>
    void load(Archive & ar, const unsigned int version) 
    {
      if(version>0) { 
	ar & window_data_size;

	std::shared_ptr<char> data_ptr(new char[window_data_size]);
	ar & boost::serialization::make_binary_object(data_ptr.get(),window_data_size);
	window_data_ptr.swap(data_ptr);

      }
    }
  BOOST_SERIALIZATION_SPLIT_MEMBER()
   
};


}  // end of namespace datatypes
}  // end of namespace uboone
}  // end of namespace fnal
}  // end of namespace gov

// This MACRO must be outside any namespaces.
BOOST_CLASS_VERSION(gov::fnal::uboone::datatypes::windowDataPMT, gov::fnal::uboone::datatypes::constants::VERSION)    

#endif /* #ifndef BOONETYPES_H */



