////////////////////////////////////////////////////////////////////////////////
/// \file UBOONEGeometryHelper_service.cc
///
/// \version $Id
/// \author  rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

// Migration note:
// Geometry --> uboone/Geometry
#include "uboone/Geometry/UBooNEGeometryHelper.h"

#include "Geometry/ChannelMapAlg.h"

// Migration note:
// Geometry --> uboone/Geometry for the two below
#include "uboone/Geometry/ChannelMapUBooNEAlg.h"

#include "TString.h"


namespace uboone
{

  UBooNEGeometryHelper::UBooNEGeometryHelper( fhicl::ParameterSet const & pset, art::ActivityRegistry & reg )
  :  fPset( pset ),
     fReg( reg ),
     fChannelMap()
  {}

  UBooNEGeometryHelper::~UBooNEGeometryHelper() throw()
  {}  
  
  void UBooNEGeometryHelper::doConfigureChannelMapAlg( const TString & detectorName,
                                                     fhicl::ParameterSet const & sortingParam,
                                                     std::vector<geo::CryostatGeo*> & c,
						     std::vector<geo::AuxDetGeo*>   & ad )
  {
    fChannelMap = std::shared_ptr<geo::ChannelMapAlg>( new geo::ChannelMapUBooNEAlg( sortingParam, fPset ) );
  }
  
  std::shared_ptr<const geo::ChannelMapAlg> UBooNEGeometryHelper::doGetChannelMapAlg() const
  {
    return fChannelMap;
  }

}

DEFINE_ART_SERVICE_INTERFACE_IMPL(uboone::UBooNEGeometryHelper, geo::ExptGeoHelperInterface)
