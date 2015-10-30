#ifndef _BEAMRUN_H
#define _BEAMRUN_H

#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include "boost/date_time/posix_time/posix_time.hpp"
#include <boost/archive/binary_oarchive.hpp>

#include "datatypes/ub_BeamHeader.h"
#include "datatypes/ub_BeamData.h"
#include "beamIFDBInterface.h"
#include "beamRunHeader.h"

namespace gov {
namespace fnal {
namespace uboone {
namespace beam {

class beamRun {

  typedef std::map<gov::fnal::uboone::datatypes::ub_BeamHeader, std::vector<gov::fnal::uboone::datatypes::ub_BeamData> > beamdatamap_t;

public:
  beamRun();
  ~beamRun();

  void StartRun(beamRunHeader& rh, boost::posix_time::ptime tstart);
  void Update(boost::posix_time::ptime tend=boost::posix_time::microsec_clock::local_time());
  void EndRun(boost::posix_time::ptime tend);

 private:
  
  void WriteData(std::string beamline, beamdatamap_t &data_map);
  void ProcessResponse(httpResponse* response, beamdatamap_t &data_map, std::string beamline,boost::posix_time::ptime tend);
  boost::posix_time::ptime ToPtime(gov::fnal::uboone::datatypes::ub_BeamHeader bh);

  std::map<std::string, boost::posix_time::ptime > fLastQueryTime;
    
   beamIFDBInterface fIFDB;

   //run info
   gov::fnal::uboone::beam::beamRunHeader fRunHeader;

   //output archive
   std::ofstream* fOut; 
   std::map<std::string, boost::archive::binary_oarchive* > fOA;

   //config parameters
   std::string fOutputDirData;
   std::string fOutputDirInfo;
   std::vector<std::string> fBeamLine;
   std::map<std::string, int> fEventTypeMap;
   std::map<std::string, int> fTimeWindowMap;
   std::map<std::string, float> fTimeOffsetMap;
   std::map<std::string, float> fTimePaddingMap;
   std::map<std::string, std::vector<std::string> > fBundle;
   std::string fIFDBURL;
   int fIFDBLatency;
   long fMWRTimeOffset;

   std::map<std::string, gov::fnal::uboone::datatypes::ub_BeamHeader> fLastBeamHeader;

   boost::posix_time::time_duration fZoneOffset;
};
}  // end of namespace beam
}  // end of namespace uboone
}  // end of namespace fnal
}  // end of namespace gov

#endif /* #ifndef BEAMRUN_H */
