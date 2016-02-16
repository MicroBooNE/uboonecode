
#include "uboone/Utilities/FileCatalogMetadataMicroBooNE.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/System/FileCatalogMetadata.h"

#include "LLMetaMaker.h"
#include <TTimeStamp.h>
#include <sstream>


//-----------------------------------------------------------------------------------------
util::LLMetaMaker::LLMetaMaker(fhicl::ParameterSet const& pset, art::ActivityRegistry &reg)
//----------------------------------------------------------------------------------------
{

  reconfigure(pset);

  reg.sPostOpenFile.watch     (this, &LLMetaMaker::postOpenFile     );

  reg.sPreBeginRun.watch      (this, &LLMetaMaker::preBeginRun      );
  reg.sPostBeginRun.watch     (this, &LLMetaMaker::postBeginRun     );

  reg.sPreProcessEvent.watch  (this, &LLMetaMaker::preProcessEvent  );
  reg.sPostProcessEvent.watch (this, &LLMetaMaker::postProcessEvent );

  reg.sPostBeginJob.watch     (this, &LLMetaMaker::postBeginJob );
  reg.sPostEndJob.watch       (this, &LLMetaMaker::postEndJob   );

  /*
  ::larlite::sam::FCLMetaData_t         _fcl_meta;
  ::larlite::sam::UBMetaData_t          _ub_meta;
  ::larlite::sam::FileCatalogMetaData_t _fcat_meta;
  ::larlite::sam::RunMetaData_t         _run_meta;
  ::larlite::sam::SAMBuiltInMetaData_t  _sam_meta;
  */

  _fcl_meta.data_tier = "larlite";
  _fcat_meta.file_format = "larlite";

}

//------------------------------------------------------------------
void util::LLMetaMaker::reconfigure(fhicl::ParameterSet const& pset)
//------------------------------------------------------------------
{}

//------------------------------------------------------------------
void util::LLMetaMaker::postBeginJob()
//------------------------------------------------------------------
{
  _sam_meta.start_time = TTimeStamp().AsString("s"); // UTC time stamp in SQL format
  art::ServiceHandle<art::FileCatalogMetadata> larmeta_handle;

  art::FileCatalogMetadata::collection_type larmeta;
  larmeta_handle->getMetadata(larmeta);

  for(auto const& key_value : larmeta) {
    if      (key_value.first == "file_type"         ) _fcat_meta.file_type           = key_value.second;
    else if (key_value.first == "group"             ) _fcat_meta.group               = key_value.second;
    else if (key_value.first == "applicationFamily" ) _fcat_meta.application_family  = key_value.second;
    else if (key_value.first == "process_name"      ) _fcat_meta.application_name    = key_value.second;
    else if (key_value.first == "applicationVersion") _fcat_meta.application_version = key_value.second;
  }
  
}

//------------------------------------------------------------------
void util::LLMetaMaker::postEndJob()
//------------------------------------------------------------------
{
  _sam_meta.end_time = TTimeStamp().AsString("s"); // UTC time stamp in SQL format

  art::ServiceHandle<util::FileCatalogMetadataMicroBooNE> ubmeta_handle;

  _fcl_meta.name     = ubmeta_handle->FCLName();
  _fcl_meta.version  = ubmeta_handle->FCLVersion();
  
  _ub_meta.project_name    = ubmeta_handle->ProjectName();
  _ub_meta.project_stage   = ubmeta_handle->ProjectStage();
  _ub_meta.project_version = ubmeta_handle->ProjectVersion();
}

//------------------------------------------------------------
void util::LLMetaMaker::preProcessEvent(const art::Event& evt)
//------------------------------------------------------------
{
  std::cout << "This is before event process" << std::endl;
}

//------------------------------------------------------------
void util::LLMetaMaker::postProcessEvent(const art::Event& evt)
//------------------------------------------------------------
{
  size_t run    = evt.run();
  size_t subrun = evt.subRun();
  size_t event  = evt.event();

  auto iter = _sam_meta.runs_m.find(run);
  if(iter == _sam_meta.runs_m.end()) {
    // New run found. Insert RunMetaData_t
    iter = (_sam_meta.runs_m.emplace(run,::larlite::sam::RunMetaData_t())).first;
    (*iter).second.run_type = "\" \"";

    art::ServiceHandle<art::FileCatalogMetadata> larmeta_handle;
    art::FileCatalogMetadata::collection_type larmeta;
    larmeta_handle->getMetadata(larmeta);
    for(auto const& key_value : larmeta) {
      if(key_value.first != "run_type") continue;
      (*iter).second.run_type = key_value.second;
      break;
    }
  }
  
  auto& run_meta = (*iter).second;
  run_meta.subruns.insert(subrun);

  if(_sam_meta.event_count == 0) _sam_meta.first_event = event;
  _sam_meta.last_event = event;
  ++_sam_meta.event_count;
}

//------------------------------------------------------
void util::LLMetaMaker::preBeginRun(art::Run const& run)
//------------------------------------------------------
{}

//------------------------------------------------------
void util::LLMetaMaker::postBeginRun(art::Run const& run)
//------------------------------------------------------
{}

//---------------------------------------------------------------
void util::LLMetaMaker::postOpenFile(const std::string& filename)
//---------------------------------------------------------------
{
  _sam_meta.parents.insert(filename);
}

//----------------------------------------------------------------------
std::string util::LLMetaMaker::GetContent(std::string stream_name) const
//----------------------------------------------------------------------
{
  std::stringstream msg;
  msg << "{\n"

      << "  \"application\": {\n"
      << "    \"family\"  : " << _fcat_meta.application_family  << ",\n"
      << "    \"name\"    : " << _fcat_meta.application_name    << ",\n"
      << "    \"version\" : " << _fcat_meta.application_version << "\n"
      << "  }\n"
      << "  \"file_format\" : \"" << _fcat_meta.file_format << "\",\n"
      << "  \"file_type\"   : "   << _fcat_meta.file_type   << ",\n"
      << "  \"group\"       : "   << _fcat_meta.group       << ",\n"

      << "  \"fcl.name\"    : \"" << _fcl_meta.name        << "\",\n"
      << "  \"fcl.version\" : \"" << _fcl_meta.version     << "\",\n"
      << "  \"data_tier\"   : \"" << _fcl_meta.data_tier   << "\",\n"
      << "  \"data_stream\" : \"" << stream_name           << "\",\n"

      << "  \"ub_project.name\"    : \"" << _ub_meta.project_name    << "\",\n"
      << "  \"ub_project.stage\"   : \"" << _ub_meta.project_stage   << "\",\n"
      << "  \"ub_project.version\" : \"" << _ub_meta.project_version << "\",\n"

      << "  \"start_time\"  : \"" << _sam_meta.start_time  << "\",\n"
      << "  \"end_time\"    : \"" << _sam_meta.end_time    << "\",\n"
      << "  \"first_event\" : "   << _sam_meta.first_event << ",\n"
      << "  \"last_event\"  : "   << _sam_meta.last_event  << ",\n"
      << "  \"event_count\" : "   << _sam_meta.event_count << ",\n"
      << "  \"parents\"     : [\n";
  
  for(auto const& parent : _sam_meta.parents) {
    
    size_t name_start = parent.rfind("/");
    if(name_start > parent.length()) name_start = 0;
    else ++name_start;
    msg << "    {  \"file_name\" : \"" << parent.substr(name_start) << "\"  }\n";
  }
  
  msg << "  ]\n"
      << "  \"runs\" : [\n";
  for(auto const& run_meta : _sam_meta.runs_m) {

    auto const& run = run_meta.first;
    auto const& run_type = run_meta.second.run_type;
    for(auto const& subrun : run_meta.second.subruns)

      msg << "    [  " << run << ",  " << subrun << ",  " << run_type << "]\n";
  }
  msg << "  ]\n"
      << "}\n";
  return msg.str();
}

namespace util{

  DEFINE_ART_SERVICE(LLMetaMaker)

} // namespace util  

