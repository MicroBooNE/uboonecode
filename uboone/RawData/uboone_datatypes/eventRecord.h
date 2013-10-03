#ifndef _UBOONETYPES_EVENTRECORD_H
#define _UBOONETYPES_EVENTRECORD_H
#include <sys/types.h>
#include <inttypes.h>
#include <vector>
#include "evttypes.h"
#include "globalHeader.h"
#include "triggerData.h"
#include "gps.h"
#include "crateHeader.h"
#include "crateData.h"
#include "crateDataPMT.h"
#include "beamHeader.h"
#include "beamData.h"

#include <boost/serialization/list.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/binary_object.hpp>

#include "constants.h"

/***
    The eventRecord is meant to house all of the components of the 
    final data format as it leaves the assembler and is written to 
    disk. The data will be written as a boost binary_archive, so 
    we may version different pieces accordingly. Along with all the
    (independently) serialized headers, we have a map that pairs 
    crate headers and crate data.
 ***/


namespace gov {
namespace fnal {
namespace uboone {
namespace datatypes {

using namespace gov::fnal::uboone;

//used for map
struct compareCrateHeader {
  bool operator() ( crateHeader lhs, crateHeader rhs) const
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
  void insertSEB(crateHeader,crateDataPMT); //in .cpp file
  //void insertSEB_PMT(crateHeader,crateDataPMT); //in .cpp file
  //void insertSEB_TPC(crateHeader,crateData); //in .cpp file
  
  globalHeader getGlobalHeader() { return global_header; }
  triggerData getTriggerData() { return trigger_data; }
  gps getGPS() { return gps_data; }
  beamHeader getBeamHeader() { return beam_header; }

  std::vector<beamData> getBeamDataVector() { return beam_data_vector; }
  std::map<crateHeader,crateData,compareCrateHeader> getSEBMap() { return seb_map; }
  std::map<crateHeader,crateDataPMT,compareCrateHeader> getSEBPMTMap() { return seb_pmt_map; }

  globalHeader* getGlobalHeaderPtr() { return &global_header; }
  triggerData* getTriggerDataPtr() { return &trigger_data; }
  gps* getGPSPtr() { return &gps_data; }
  beamHeader* getBeamHeaderPtr() { return &beam_header; }

  int getSEBMap_size() { return seb_map.size(); }
  void clearSEBMap() { seb_map.clear(); }

  int getSEBPMTMap_size() { return seb_pmt_map.size(); }
  void clearSEBPMTMap() { seb_pmt_map.clear(); }

  int getBeamDataVecotr_size() { return beam_data_vector.size(); }
  void clearBeamDataVector() { beam_data_vector.clear(); }

  uint8_t getIOMode() { return er_IO_mode; }
  void updateIOMode(uint8_t); //in .cpp file

  void decompress();

 private:

  globalHeader global_header;
  triggerData trigger_data;
  gps gps_data;
  std::map<crateHeader,crateData,compareCrateHeader> seb_map;
  std::map<crateHeader,crateDataPMT,compareCrateHeader> seb_pmt_map;
  beamHeader beam_header;
  std::vector<beamData> beam_data_vector;
  
  uint8_t er_IO_mode;

  friend class boost::serialization::access;
  
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      if(version>1)
	ar & er_IO_mode
	   & global_header
	   & trigger_data
	   & gps_data
	   & beam_header & beam_data_vector //beam stuff...empty at first, added in later
	   & seb_map
	   & seb_pmt_map;

      else if(version>0)
	ar & er_IO_mode
	   & global_header
	   & seb_map;
    }

};

}  // end of namespace datatypes
}  // end of namespace uboone
}  // end of namespace fnal
}  // end of namespace gov

// This MACRO must be outside any namespaces.
BOOST_CLASS_VERSION(gov::fnal::uboone::datatypes::eventRecord, gov::fnal::uboone::datatypes::constants::VERSION)    

#endif /* #ifndef BOONETYPES_H */
