#include "bdaq_get.h"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

using namespace gov::fnal::uboone::datatypes;
using namespace gov::fnal::uboone::beam;
using namespace std;
using namespace boost::posix_time;
using namespace boost::program_options;
using namespace boost::filesystem;

int main(int ac, char* av[]) 
{
  /**********************************************************************
   *
   * Parse options
   *
   *********************************************************************/
  vector<uint32_t> t0;//begin time
  vector<uint32_t> t1;//end time

  options_description opt("Options");
  opt.add_options()
    ("help,h", "Print help message")
    ("begin-time,i",value< vector<uint32_t> >(&t0)->multitoken(),"Begin (sub)run time in seconds milliseconds (unix time); takes up to two arguments")
    ("end-time,e",value< vector<uint32_t> >(&t1)->multitoken(),"End (sub)run time in seconds milliseconds (unix time); takes up to two arguments")
    ("run-number,r",value<int>()->default_value(99999),"Run number.")
    ("subrun-number,s",value<int>()->default_value(9999),"Subrun number.")
    ("fhicl-file,f",value<string>()->default_value(""),"Configuration fhicl file.");

  variables_map vm;

  try {
    store(parse_command_line(ac,av,opt),vm);
    notify(vm);
    if (vm.count("help")) {
      cerr<<opt<<endl;
      return 1;
    }
    if (!vm.count("begin-time") || !vm.count("end-time") ) {
      cerr<<"ERROR: Need to specify begin and end time!"<<endl;
      cerr<<opt<<endl;
      return 1;
    }
    if (!vm.count("run-number") || !vm.count("subrun-number") ) {
      cerr<<"ERROR: Need to specify run/subrun number!"<<endl;
      cerr<<opt<<endl;
      return 1;
    }
    if (t0.size()>2 || t1.size()>2) {
      cerr<<"Begin time and end time take up to two arguments for seconds and milliseconds (unix time)"<<endl;
      cerr<<opt<<endl;
      return 1;
    }
  } catch (error& e) {
    cerr << "ERROR: " <<e.what()<<endl<<endl;
    cerr << opt <<endl;
    return 1;
  }
  t0.resize(2);
  t1.resize(2);

  if (vm.count("fhicl-file")) {
    boost::filesystem::path fpath(vm["fhicl-file"].as<string>());
    std::string fhiclPath("BEAMDAQ_CONFIG_PATH=");
    if (fpath.parent_path().string() !="") 
      fhiclPath.append(fpath.parent_path().string());
    else 
      fhiclPath.append("./");

    //    std::cout<<"Going with path= "<<fhiclPath<<std::endl;
    ::putenv(const_cast<char *>(fhiclPath.c_str()));
    
    std::string fhiclFile("BEAMDAQ_CONFIG_FILE=");
    fhiclFile.append(fpath.filename().string());
    ::putenv(const_cast<char *>(fhiclFile.c_str()));

  }
  beamDAQConfig* bdconfig=beamDAQConfig::GetInstance();

  //start message facility
  mf::StartMessageFacility( mf::MessageFacilityService::MultiThread,
			    bdconfig->GetParameterSet().get<fhicl::ParameterSet>("message_facility"));
  mf::SetApplicationName( "BeamDAQ" );

  int run=vm["run-number"].as<int>();
  int subrun=vm["subrun-number"].as<int>();
  uint32_t t0_s  = t0.at(0);
  uint32_t t0_ms = t0.at(1);
  uint32_t t1_s  = t1.at(0);
  uint32_t t1_ms = t1.at(1);

  ptime pt0=from_time_t(t0_s)+hours(fZoneOffset.hours())+microseconds(t0_ms);
  ptime pt1=from_time_t(t1_s)+hours(fZoneOffset.hours())+microseconds(t1_ms);

  if (pt1-pt0>hours(bdconfig->GetMaxRunLength())) {
    mf::LogWarning("")<<"Refusing to start a run longer than "
		      <<bdconfig->GetMaxRunLength()  
		      <<" hours!";
    return 1;
  }
  if (pt0-pt1>seconds(0)) {
    mf::LogWarning("")<<"Refusing to start a run with t1>t0";
    return 1;
  }

  beamRun brm;
  beamRunHeader rh;
  rh.fRun=run;
  rh.fSubRun=subrun;
  brm.StartRun(rh,pt0);
  brm.EndRun(pt1);

  return 0;
}
