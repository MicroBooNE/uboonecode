////////////////////////////////////////////////////////////////////////
// Name:  FileCatalogMetadataMicroBooNE_service.cc.  
//
// Purpose:  Implementation for FileCatalogMetadataMicroBooNE.
//
// Created:  28-Oct-2014,  H. Greenlee
//
////////////////////////////////////////////////////////////////////////

#include "uboone/Utilities/FileCatalogMetadataMicroBooNE.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/System/FileCatalogMetadata.h"

//--------------------------------------------------------------------
// Constructor.

util::FileCatalogMetadataMicroBooNE::
FileCatalogMetadataMicroBooNE(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
{
  // Get parameters.

  fFCLName = pset.get<std::string>("FCLName");
  fFCLVersion = pset.get<std::string>("FCLVersion");
  fProjectName = pset.get<std::string>("ProjectName");
  fProjectStage = pset.get<std::string>("ProjectStage");
  fProjectVersion = pset.get<std::string>("ProjectVersion");    

  // Register for callbacks.

  reg.sPostBeginJob.watch(this, &FileCatalogMetadataMicroBooNE::postBeginJob);
}

//--------------------------------------------------------------------
// PostBeginJob callback.
// Insert per-job metadata via FileCatalogMetadata service.
void util::FileCatalogMetadataMicroBooNE::postBeginJob()
{
  // Get art metadata service.

  art::ServiceHandle<art::FileCatalogMetadata> mds;

  // Add metadata.

  mds->addMetadata("fclName", fFCLName);
  mds->addMetadata("fclVersion", fFCLVersion);
  mds->addMetadata("ubProjectName", fProjectName);
  mds->addMetadata("ubProjectStage", fProjectStage);
  mds->addMetadata("ubProjectVersion", fProjectVersion);
}

DEFINE_ART_SERVICE(util::FileCatalogMetadataMicroBooNE)
