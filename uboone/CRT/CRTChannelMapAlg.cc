#include "uboone/Geometry/CRTChannelMapAlg.hh"
#include "larcore/Geometry/AuxDetChannelMapAlg.h"
#include "larcore/Geometry/AuxDetGeo.h"
#include "larcore/Geometry/AuxDetGeometryCore.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <string>
#include <vector>

namespace uboone {

  //----------------------------------------------------------------------------
  CRTChannelMapAlg::CRTChannelMapAlg( fhicl::ParameterSet const& pvals,
    fhicl::ParameterSet const& sortingParameters ):fNModules(72), fNStrips(17)
  {
  }

  //----------------------------------------------------------------------------
  CRTChannelMapAlg::~CRTChannelMapAlg(){}

  void CRTChannelMapAlg::Initialize(geo::AuxDetGeometryData_t& geodata){
    unsigned nDets=0;
    for(std::vector<geo::AuxDetGeo*>::iterator it = geodata.auxDets.begin();
      it!= geodata.auxDets.end(); ++it)
    {
      geo::AuxDetGeo* geo = *it;
      for (uint32_t module=0;module<=this->fNModules;module++){
        for (uint32_t strip=0;strip<=this->fNStrips; strip++){
          char entryname[50];
          sprintf(entryname,"Module_%d_strip_%d",module, strip);
          if(geo->Name().find(entryname) != std::string::npos){
            mf::LogInfo("CRTChannelMapAlg")<<"Found CRT Panel in Geometry: "<<entryname;
            nDets++;
          }
        }
      }
    }
    mf::LogInfo("CRTChannelMapAlg")<<"Number of Detectors: "<<nDets;
  }

  void CRTChannelMapAlg::Uninitialize(){
  }

  uint32_t CRTChannelMapAlg::PositionToAuxDetChannel(
   double const worldLoc[3],
   std::vector<geo::AuxDetGeo*> const& auxDets,
   size_t& ad,
   size_t& sv) const{
     uint32_t channel = UINT_MAX;
     ad = this->NearestAuxDet(worldLoc, auxDets);
     sv = this->NearestSensitiveAuxDet(worldLoc, auxDets, ad);
     return 0;
  }

  const TVector3 CRTChannelMapAlg::AuxDetChannelToPosition(
    uint32_t const& channel,
    std::string const& auxDetName,
    std::vector<geo::AuxDetGeo*> const& auxDets) const{
      return TVector3(0., 0., 0.);
  }
} // namespace
