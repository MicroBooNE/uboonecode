#include "uboone/CRT/CRTChannelMapAlg.hh"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/AuxDetGeometryCore.h"

#include <sstream>

namespace crt {

  CRTChannelMapAlg::CRTChannelMapAlg( fhicl::ParameterSet const& pvals,
    fhicl::ParameterSet const& sortingParameters ) : 
    fSorter(sortingParameters)
  {
  
  }

  CRTChannelMapAlg::~CRTChannelMapAlg()
  {
    return;
  }

  void CRTChannelMapAlg::Initialize(geo::AuxDetGeometryData_t& geodata)
  {
    /**
    * Sets up the following variables:
    *
    * fSorter: Sorts the geodata by name. We're assuming that the geodata is
    *          introduced without any sorting.
    *
    * fADGeoToName: maps the channel number to name
    *
    * fNameToADGeo: maps the name to channel number. This seems degenerate.
    *
    * fADGeoToChannelAndSV: largely unused since the geometry has 1 sensitive volume per detector
    *
    * \todo: The two for loops can be exchanged for a single iterator pattern.
    *        This will make the code more extensible.
    **/
    this->Uninitialize();
    std::vector<geo::AuxDetGeo*>& adgeo = geodata.auxDets;
    fSorter.SortAuxDets(adgeo);

    fADGeoToName.clear();
    fNameToADGeo.clear();
    fADGeoToChannelAndSV.clear();

    uint32_t index=0;
    for(uint32_t mod=0; mod<this->fSorter.GetNModules(); ++mod)
    {
      for(uint32_t strip=0; strip<this->fSorter.GetNStripsPerModule();++strip)
      {
        std::ostringstream stream;
        stream<<"volAuxDet_Module_"<<mod<<"_strip_"<<strip;
        fADGeoToName[index] = stream.str();
        fNameToADGeo[stream.str()] = index;
        fADGeoToChannelAndSV[index].push_back(std::make_pair(index,0));
        ++index;
      }
    }
  }

  uint32_t CRTChannelMapAlg::PositionToAuxDetChannel(
    double const worldLoc[3],
    std::vector<geo::AuxDetGeo*> const& auxDets,
    size_t& ad,
    size_t& sv) const
  {
    /**
    * Given a three space location, finds the nearest detector and 
    * returns the channel.
    *
    * \throws cet::exception if the aux det cannot be found or if the position
    * is not inside any of the detectors.
    **/

    uint32_t channel = UINT_MAX;
    ad = 0;
    sv = this->NearestSensitiveAuxDet(worldLoc, auxDets, ad);
    double svOrigin[3] = {0, 0, 0};
    double localOrigin[3] = {0, 0, 0};
    auxDets[ad]->SensitiveVolume(sv).LocalToWorld(localOrigin, svOrigin);
    auto gnItr = fADGeoToName.find(ad);
    if (gnItr != fADGeoToName.end())
    {
      auto csvItr = fADGeoToChannelAndSV.find(ad);
      if (csvItr == fADGeoToChannelAndSV.end()) 
      {
        throw cet::exception("CRTChannelMapAlg")
        << "No entry in channel and sensitive volume map for AuxDet index "
        << ad;
      }
      channel = 2 * sv + 0;
    }
    if (channel == UINT_MAX) 
    {
      throw cet::exception("CRTChannelMapAlg") << "position ("
      << worldLoc[0] << "," << worldLoc[1] << "," << worldLoc[2]
      << ") does not correspond to any AuxDet";
    }
    return channel;
  }

  const TVector3 CRTChannelMapAlg::AuxDetChannelToPosition(
    uint32_t const& channel,
    std::string const& auxDetName,
    std::vector<geo::AuxDetGeo*> const& auxDets) const
  {
    /**
    * Given a channel number and auxdet name, gets the position of the 
    * corresponding detector
    *
    * \throws cet::exception if the name is not registered in Initialize()
    * or if the channel number does not correspond to a detector.
    **/
    double x = 0;
    double y = 0;
    double z = 0;
    size_t ad = UINT_MAX;
    if (fNameToADGeo.count(auxDetName) > 0) 
    {
      ad = fNameToADGeo.find(auxDetName)->second;
    }
    else 
    {
      throw cet::exception("CRTChannelMapAlg")
      << "No AuxDetGeo with name " << auxDetName;
    }
    auto csvItr = fADGeoToChannelAndSV.find(ad);
    if (csvItr == fADGeoToChannelAndSV.end()) 
    {
      throw cet::exception("CRTChannelMapAlg")
      << "No entry in channel and sensitive volume"
      << " map for AuxDet index " << ad << " bail";
    }
    double svOrigin[3] = {0, 0, 0};
    double localOrigin[3] = {0, 0, 0};
    for (auto csv : csvItr->second) 
    {
      if (csv.first == channel) 
      {
        auxDets[ad]->SensitiveVolume(csv.second).LocalToWorld(localOrigin,
                    svOrigin);
        x = svOrigin[0];
        y = svOrigin[1];
        z = svOrigin[2];
        break;
      }
    }
    return TVector3(x, y, z);
  }
}
