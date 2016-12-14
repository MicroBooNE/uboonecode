/**
 * \class CRTGeometryHelper
 *
 * \ingroup crt
 *
 * \brief Interface class for the crt channel map.
 *
 * See `geo::AuxDetExptGeoHelperInterface` for full explanation.
 *
 * \author $Author: Kevin Wierman<kevin.wierman@pnnl.gov> 
 *
 * \version $Revision: 1.0 $
 *
 * \date $Date: 2016/12/12 $
 *
 * Contact: kevin.wierman@pnnl.gov
 *
 * Created on: Tuesday, December 13, 2016
 *
 */

#ifndef CRTGeometryHelper_HH_
#define CRTGeometryHelper_HH_

// framework libraries
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/AuxDetExptGeoHelperInterface.h"
#include "larcore/Geometry/AuxDetChannelMapAlg.h"

#include <memory> //For std::shared_ptr


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
