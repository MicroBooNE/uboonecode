#ifndef CRTChannelMapAlg_hh_
#define CRTChannelMapAlg_hh_

#include "larcore/Geometry/GeoObjectSorterStandard.h"
#include "larcore/Geometry/AuxDetChannelMapAlg.h"
#include "larcore/Geometry/AuxDetGeo.h"
#include "fhiclcpp/ParameterSet.h"
#include "TVector3.h"
#include <iostream>
#include <vector>
#include <set>


namespace geo{
  struct AuxDetGeometryData_t;
  class AuxDetGeo;
}

namespace crt {

  class CRTChannelMapAlg : public geo::AuxDetChannelMapAlg {
    //TODO: Redefine the mapping scheme
  protected:
    uint32_t ChannelNumberFromModuleAndPanel(uint32_t const& module,
                                              uint32_t const& panel) const;
  public:

    CRTChannelMapAlg(fhicl::ParameterSet const& pset,
                        fhicl::ParameterSet const& sortingParameters );

    ~CRTChannelMapAlg();

    void Initialize(geo::AuxDetGeometryData_t& geodata);
    void Uninitialize();

    uint32_t PositionToAuxDetChannel(double const  worldLoc[3],
                                     std::vector<geo::AuxDetGeo*> const& auxDets,
                                     size_t& ad,
                                     size_t& sv) const;

    const TVector3 AuxDetChannelToPosition(uint32_t const& channel,
                                           std::string const& auxDetName,
                                           std::vector<geo::AuxDetGeo*> const& auxDets) const;
  };
}

#endif
