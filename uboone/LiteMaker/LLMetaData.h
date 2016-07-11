
#ifndef __LLMETADATA_H__
#define __LLMETADATA_H__

#include <string>
#include <map>
#include <set>

namespace larlite {

  namespace sam {

    struct FCLMetaData_t {
      std::string name;
      std::string version;
      std::string data_tier;
      std::string data_stream;
      FCLMetaData_t() : name        (" ")
		      , version     (" ")
		      , data_tier   (" ")
		      , data_stream (" ")
      {}
    };

    struct UBMetaData_t {
      std::string project_name;
      std::string project_stage;
      std::string project_version;
      UBMetaData_t() : project_name    (" ")
		     , project_stage   (" ")
		     , project_version (" ")
      {}
    };

    struct FileCatalogMetaData_t {
      std::string file_format;
      std::string file_type;
      std::string group;
      std::string application_family;
      std::string application_name;
      std::string application_version;
      FileCatalogMetaData_t() : file_format ("larlite")
			      , file_type   (" ")
			      , group       (" ")
			      , application_family  (" ")
			      , application_name    (" ")
			      , application_version (" ")

      {}
    };

    struct RunMetaData_t {
      std::string run_type;
      std::set<size_t> subruns;
      RunMetaData_t() : run_type(" ")
		      , subruns()
      {}
    };
    
    struct SAMBuiltInMetaData_t {
      size_t first_event;
      size_t last_event;
      size_t event_count;
      std::string start_time;
      std::string end_time;
      std::set<std::string> parents;
      std::map<size_t,larlite::sam::RunMetaData_t> runs_m;
      SAMBuiltInMetaData_t() : first_event (0)
			     , last_event  (0)
			     , event_count (0)
			     , start_time  (" ")
			     , end_time    (" ")
			     , parents     ()
			     , runs_m      ()
      {}
    };

  }
}
#endif
