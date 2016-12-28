////////////////////////////////////////////////////////////////////////////////
/// \file UBOONEGeometryHelper_service.cc
///
/// \version $Id
/// \author  rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

// Migration note:
// Geometry --> uboone/Geometry
#include "uboone/Geometry/UBooNEGeometryHelper.h"

#include "larcore/Geometry/ChannelMapAlg.h"
#include "larcore/Geometry/GeometryCore.h" // larcore. geo::GeometryData_t

// Migration note:
// Geometry --> uboone/Geometry for the two below
#include "uboone/Geometry/ChannelMapUBooNEAlg.h"

#include "TString.h"


namespace uboone
{

  UBooNEGeometryHelper::UBooNEGeometryHelper( fhicl::ParameterSet const & pset, art::ActivityRegistry & reg )
  :  fPset( pset )
     //fReg( reg )
  {}

  UBooNEGeometryHelper::~UBooNEGeometryHelper() throw()
  {}  
  
  void UBooNEGeometryHelper::doConfigureChannelMapAlg( fhicl::ParameterSet const & sortingParameters, geo::GeometryCore* geom ) 
  {
    fChannelMap.reset();
    std::string const detectorName = geom->DetectorName();

    if ( detectorName.find("microboone") == std::string::npos ) {
      std::cout << __PRETTY_FUNCTION__ << ": WARNING USING CHANNEL MAP ALG WITH NON-MICROBOONE GEO!" << std::endl;
    }

    fChannelMap = std::make_shared<geo::ChannelMapUBooNEAlg>( fPset, sortingParameters );

    if ( fChannelMap )
      {
        geom->ApplyChannelMap(fChannelMap); // calls Initialize(fGeoData) for us
      }

  }
  
  std::shared_ptr<const geo::ChannelMapAlg> UBooNEGeometryHelper::doGetChannelMapAlg() const
  {
    return fChannelMap;
  }

}

DEFINE_ART_SERVICE_INTERFACE_IMPL(uboone::UBooNEGeometryHelper, geo::ExptGeoHelperInterface)
