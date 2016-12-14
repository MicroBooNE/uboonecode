/**
 * \class CRTChannelMapAlg
 *
 * \ingroup crt
 *
 * \brief Provides an interface between position, channel and <module,strip>
 *
 * At any point in time, a class may have access to either the navigation position,
 * the channel number XOR the module and strip number for the CRT. This is meant to
 * provide a conversion between the three.
 *
 * \note At the moment, positions is center of the panel. Therefore, in order to do
 * panel intersections, this will require some work.
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

#ifndef CRTChannelMapAlg_hh_
#define CRTChannelMapAlg_hh_

#include "uboone/CRT/CRTGeoObjectSorter.hh"
#include "larcore/Geometry/GeoObjectSorterStandard.h"
#include "larcore/Geometry/AuxDetChannelMapAlg.h"
#include "larcore/Geometry/AuxDetGeo.h"
#include "fhiclcpp/ParameterSet.h"
#include "TVector3.h"
#include <iostream>
#include <vector>
#include <set>


namespace crt {

  class CRTChannelMapAlg : public geo::AuxDetChannelMapAlg {

    /// Sorts the list into the order of Module:Strip
    CRTGeoObjectSorter fSorter;

  public:

    /// Art style constructor
    CRTChannelMapAlg(fhicl::ParameterSet const& pset,
                        fhicl::ParameterSet const& sortingParameters );

    /// Default d'tor
    ~CRTChannelMapAlg();

    /// Sets up channel sorting based off the auxiliary detector information
    void Initialize(geo::AuxDetGeometryData_t& geodata);

    /// Nothing to uniit
    void Uninitialize() {}
    
    /// Given a position, return the channel number aux det number and sens det number.
    uint32_t PositionToAuxDetChannel(double const  worldLoc[3],
                                     std::vector<geo::AuxDetGeo*> const& auxDets,
                                     size_t& ad,
                                     size_t& sv) const;

    /// Given an aux det, returns the center of the detector
    const TVector3 AuxDetChannelToPosition(uint32_t const& channel,
                                           std::string const& auxDetName,
                                           std::vector<geo::AuxDetGeo*> const& auxDets) const;
  };
}

#endif
