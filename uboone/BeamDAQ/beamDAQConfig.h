#ifndef _BEAMDAQCONFIG_H
#define _BEAMDAQCONFIG_H

#include <string>
#include <vector>
#include <map>

#include "fhiclcpp/ParameterSet.h"

namespace gov {
namespace fnal {
namespace uboone {
namespace beam {

class beamDAQConfig {
  
 public:
  static beamDAQConfig* GetInstance();

  std::string GetDataOutputDir() {return fOutputDirData;};
  std::string GetInfoOutputDir() {return fOutputDirInfo;};
  std::vector<std::string> GetBeamLineList() {return fBeamLine;};
  std::map<std::string, int> GetEventTypeMap() {return fEventTypeMap;};
  std::map<std::string, int> GetTimeWindowMap() {return fTimeWindowMap;};
  std::map<std::string, float> GetTimeOffsetMap() {return fTimeOffsetMap;};
  std::map<std::string, float> GetTimePaddingMap() {return fTimePaddingMap;};
  std::map<std::string, std::vector<std::string> > GetBundles() {return fBundle;};
  std::string GetIFDBURL() {return fIFDBURL;};
  int GetMaxRunLength() {return fMaxRunLength;};
  int GetIFDBLatency() {return fIFDBLatency;};
  fhicl::ParameterSet GetParameterSet() {return fPSet;};
  long GetMWRTimeOffset() {return fMWRTimeOffset;};

 private:
  beamDAQConfig();

  beamDAQConfig(beamDAQConfig const&);  // Don't Implement
  void operator=(beamDAQConfig const&); // Don't implement

  static beamDAQConfig* fInstance;

  fhicl::ParameterSet fPSet;
  std::string fOutputDirData;
  std::string fOutputDirInfo;
  std::vector<std::string> fBeamLine;
  std::map<std::string, std::vector<std::string> > fBundle;
  std::string fIFDBURL;
  std::map<std::string, int> fEventTypeMap;
  std::map<std::string, int> fTimeWindowMap;
  std::map<std::string, float> fTimeOffsetMap;
  std::map<std::string, float> fTimePaddingMap;
  int fIFDBLatency;
  int fMaxRunLength;
  long fMWRTimeOffset;

};
}  // end of namespace beam
}  // end of namespace uboone
}  // end of namespace fnal
}  // end of namespace gov

#endif /* #ifndef BEAMDAQCONFIG_H */
