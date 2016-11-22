#ifndef GEO_CHANNELMAP_UBOONE_ALG_H
#define GEO_CHANNELMAP_UBOONE_ALG_H

#include <vector>
#include <set>
#include <iostream>

#include "fhiclcpp/ParameterSet.h"

#include "larcore/Geometry/AuxDetChannelMapAlg.h"
#include "larcore/Geometry/GeoObjectSorterStandard.h"
#include "larcore/Geometry/AuxDetGeo.h"

namespace geo{

  // forward-declaration from geometry
  struct AuxDetGeometryData_t;
  class AuxDetGeo;
}

namespace uboone {

  class CRTChannelMapAlg : public geo::AuxDetChannelMapAlg {
    uint32_t fNModules;
    uint32_t fNStrips;//Should really be panels per module.
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
#endif // GEO_CHANNELMAPSTANDARDALG_H
