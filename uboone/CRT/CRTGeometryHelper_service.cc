#include "uboone/CRT/CRTGeometryHelper.hh"
#include "uboone/CRT/CRTChannelMapAlg.hh"
#include "larcore/Geometry/AuxDetGeometryCore.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include <vector>
#include "larcore/Geometry/AuxDetGeo.h"


namespace crt
{

  CRTGeometryHelper::CRTGeometryHelper( fhicl::ParameterSet const & pset, art::ActivityRegistry & reg ) :  fPset( pset ){}

  CRTGeometryHelper::~CRTGeometryHelper(){}

  void CRTGeometryHelper::doConfigureAuxDetChannelMapAlg( fhicl::ParameterSet const & sortingParameters, geo::AuxDetGeometryCore* geom )
  {
    fChannelMap.reset();
    fChannelMap = std::make_shared<CRTChannelMapAlg>( fPset, sortingParameters );
    if ( fChannelMap )
      	geom->ApplyChannelMap(fChannelMap);

  }

  CRTGeometryHelper::AuxDetChannelMapAlgPtr_t CRTGeometryHelper::doGetAuxDetChannelMapAlg() const
  {
    return fChannelMap;
  }

}

DEFINE_ART_SERVICE_INTERFACE_IMPL(crt::CRTGeometryHelper, geo::AuxDetExptGeoHelperInterface)
