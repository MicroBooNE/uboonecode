////////////////////////////////////////////////////////////////////////
// Name:  FileCatalogMetadataMicroBooNE_service.cc.  
//
// Purpose:  Art service adds microboone-specific per-job sam metadata.
//
//           FCL parameters:
//
//           FCLName        - FCL file name.
//           FCLVersion     - FCL file version.
//           ProjectName    - Project name.
//           ProjectStage   - Project stage.
//           ProjectVersion - Project version.
//
//           Above values are recorded in internal sam metadata generated
//           by art program.
//
//           This service does not have user-callable methods.  Simply
//           add to an art configuration in services.user block of job
//           file.
//
// Created:  28-Oct-2014,  H. Greenlee
//
////////////////////////////////////////////////////////////////////////

#include <string>
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/System/FileCatalogMetadata.h"

namespace util {

  // Class declaration.

  class FileCatalogMetadataMicroBooNE
  {
  public:

    // Constructor, destructor.

    FileCatalogMetadataMicroBooNE(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~FileCatalogMetadataMicroBooNE() = default;

  private:

    // Callbacks.

    void postBeginJob();

    // Data members.

    std::string fFCLName;
    std::string fFCLVersion;
    std::string fProjectName;
    std::string fProjectStage;
    std::string fProjectVersion;
  };

  //--------------------------------------------------------------------
  // Constructor.

  FileCatalogMetadataMicroBooNE::
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

} // namespace util

DECLARE_ART_SERVICE(util::FileCatalogMetadataMicroBooNE, LEGACY)
DEFINE_ART_SERVICE(util::FileCatalogMetadataMicroBooNE)
