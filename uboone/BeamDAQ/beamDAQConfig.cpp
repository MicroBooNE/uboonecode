#include "beamDAQConfig.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <exception>

#include "fhiclcpp/make_ParameterSet.h"
#include "cetlib/filepath_maker.h"

using namespace gov::fnal::uboone::beam;

beamDAQConfig* beamDAQConfig::fInstance = NULL;

beamDAQConfig::beamDAQConfig()
{
  //load config
  cet::filepath_lookup fpath(getenv("FHICL_FILE_PATH"));
  std::string fcfgname="beamdaq_config.fcl";
  if (getenv("BEAMDAQ_CONFIG_FILE"))
     fcfgname=getenv("BEAMDAQ_CONFIG_FILE");

   try {
    fhicl::make_ParameterSet(fcfgname, fpath, fPSet);
  } catch ( std::exception& e ) {
    std::cerr<<"Unable to locate "<<fcfgname<<" in "
	     <<"$FHICL_FILE_PATH"<<std::endl
	     <<e.what()<<std::endl;
    exit(0);
  }

  //set parameters
  try {
    fOutputDirData=fPSet.get<std::string>("output_dir_data");
    fOutputDirInfo=fPSet.get<std::string>("output_dir_info");
    fBeamLine=fPSet.get<std::vector<std::string> >("beamlines");
    for (unsigned int i=0;i<fBeamLine.size();i++) 
      fBundle[fBeamLine[i]]=fPSet.get<std::vector<std::string> >(fBeamLine[i]);
      
    fIFDBURL=fPSet.get<std::string>("ifdb_url");
    fMaxRunLength=fPSet.get<int>("max_run_length");
    fIFDBLatency=fPSet.get<int>("ifdb_latency");
    fMWRTimeOffset=fPSet.get<long>("mwr_time_offset");
    
    fhicl::ParameterSet eventType=fPSet.get<fhicl::ParameterSet>("event_type");
    
    std::vector<std::string> keys=eventType.get_keys();
    for (unsigned int i=0;i<keys.size();i++) {
      int eid=eventType.get<int>(keys[i]);
      std::pair<std::string,int> p(keys[i],eid);
      fEventTypeMap.insert(p);
    }
    fhicl::ParameterSet timeWindow=fPSet.get<fhicl::ParameterSet>("time_window");
    
    keys=timeWindow.get_keys();
    for (unsigned int i=0;i<keys.size();i++) {
      int tw=timeWindow.get<int>(keys[i]);
      std::pair<std::string,int> p(keys[i],tw);
      fTimeWindowMap.insert(p);
    }
  } catch ( std::exception& e ) {
    std::cerr<<"Bad "<<fcfgname<<" file. \n"<<e.what()<<std::endl;
    exit(0);
  }

}

beamDAQConfig* beamDAQConfig::GetInstance()
{
  if ( !fInstance ) fInstance=new beamDAQConfig;

  return fInstance;
}
