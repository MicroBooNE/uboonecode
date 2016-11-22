#ifndef CRTGeometryHelper_HH_
#define CRTGeometryHelper_HH_

// framework libraries
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/AuxDetExptGeoHelperInterface.h"
#include "larcore/Geometry/AuxDetChannelMapAlg.h"

#include <memory> //For std::shared_ptr

namespace geo{
  class AuxDetChannelMapAlg;
}

namespace crt
{

  class CRTGeometryHelper: public geo::AuxDetExptGeoHelperInterface {
  public:
    CRTGeometryHelper(fhicl::ParameterSet const& pset,
                      art::ActivityRegistry& reg);
    ~CRTGeometryHelper();
  private:

    virtual void doConfigureAuxDetChannelMapAlg(fhicl::ParameterSet const & sortingParameters,
                                        geo::AuxDetGeometryCore* geom) override;

    virtual AuxDetChannelMapAlgPtr_t doGetAuxDetChannelMapAlg() const override;

    std::shared_ptr<geo::AuxDetChannelMapAlg> fChannelMap;
    fhicl::ParameterSet fPset;
  };
}
DECLARE_ART_SERVICE_INTERFACE_IMPL(crt::CRTGeometryHelper, geo::AuxDetExptGeoHelperInterface, LEGACY)

#endif // define
