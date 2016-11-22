#include "uboone/CRT/CRTChannelMapAlg.hh"
//This next one is for debugging purposes only
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/AuxDetGeometryCore.h"

namespace crt {

  //---------------------------------------------------------------------------
  CRTChannelMapAlg::CRTChannelMapAlg( fhicl::ParameterSet const& pvals,
    fhicl::ParameterSet const& sortingParameters ) {
  }

  //---------------------------------------------------------------------------
  CRTChannelMapAlg::~CRTChannelMapAlg(){
    return;
  }

  void CRTChannelMapAlg::Initialize(geo::AuxDetGeometryData_t& geodata){
    for(auto it = geodata.auxDets.begin(); it!= geodata.auxdets.end(); ++it){
      mf::LogInfo("CRTChannelMapAlg")<<it->Name();
    }
  }

  void CRTChannelMapAlg::Uninitialize(){
    //TODO: Release Map
  }
  

  uint32_t CRTChannelMapAlg::PositionToAuxDetChannel(
    double const worldLoc[3],
    std::vector<geo::AuxDetGeo*> const& auxDets,
    size_t& ad,
    size_t& sv) const{
      ad = this->NearestAuxDet(worldLoc, auxDets);
      sv = this->NearestSensitiveAuxDet(worldLoc, auxDets, ad);
      return 0;
  }

  const TVector3 CRTChannelMapAlg::AuxDetChannelToPosition(
    uint32_t const& channel,
    std::string const& auxDetName,
    std::vector<geo::AuxDetGeo*> const& auxDets) const{
    mf::LogInfo("CRTChannelMapAlg")<<auxDetName;
    return TVector3(0.,0.,0.);
  }
  
}
