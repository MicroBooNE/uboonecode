////////////////////////////////////////////////////////////////////////
// Name: TFileCatalogMetadataMicroBooNE.h
//
// A struct datatype to hold the metadata information as it is extracted
// from various places.
//
// Created: 15-Jan-2015,  S. Gollapinni
//
////////////////////////////////////////////////////////////////////////
#ifndef TFILEMETADATAMICROBOONE_H
#define TFILEMETADATAMICROBOONE_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "art/Framework/Principal/SubRunPrincipal.h"
#include "art/Framework/Principal/Event.h"
#include <iostream>
#include <fstream>

using namespace std;

namespace util{

  class TFileMetadataMicroBooNE
  {
  public:
    TFileMetadataMicroBooNE(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~TFileMetadataMicroBooNE();
      
    void reconfigure(fhicl::ParameterSet const& p);
    
    struct metadata {
      std::tuple<std::string, std::string, std::string> fapplication;
      //no crc information yet
      //std::vector<std::string> fcrcs;
      std::string fdata_tier;
      time_t fend_time;
      unsigned int fevent_count=0;
      std::string ffile_format;
      std::string ffile_type;
      art::EventNumber_t ffirst_event=0;
      std::string fgroup; 
      art::EventNumber_t flast_event=0;
      std::set<std::string> fParents;
      std::vector<std::tuple<art::RunNumber_t,art::SubRunNumber_t,std::string>> fruns;
      time_t fstart_time=0; 
      //don't include stream_name yet
      //std::string fstream_name; 
      std::string ffcl_name;
      std::string ffcl_version;
      std::string fproject_name;
      std::string fproject_stage;
      std::string fproject_version;
    };
    
    metadata md;           
    std::set<art::SubRunNumber_t> fSubRunNumbers;  
    std::map<std::string,std::string> mdmap;    

  private:
   
    // Callbacks.
    void postBeginJob();
    void postOpenFile(std::string const& fn);
    void postEvent(art::Event const& ev);
    void postEndJob();

    // Private member functions.

    // Data members.

    // Fcl parameters.
    bool fGenerateTFileMetadata;  
    std::string frunType;                     

  }; // class TFileMetadataMicroBooNE

} //namespace utils
	
DECLARE_ART_SERVICE(util::TFileMetadataMicroBooNE, LEGACY)
	
#endif

