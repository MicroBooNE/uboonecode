#include "uboone/CRT/CRTChannelMapAlg.hh"

namespace crt {

  //---------------------------------------------------------------------------
  CRTChannelMapAlg::CRTChannelMapAlg( fhicl::ParameterSet const& pvals,
    fhicl::ParameterSet const& sortingParameters ) : 
    fNModules(72), fNStrips(17){
  }

  //---------------------------------------------------------------------------
  CRTChannelMapAlg::~CRTChannelMapAlg(){}

  void CRTChannelMapAlg::Initialize(geo::AuxDetGeometryData_t& geodata){
    // TODO: Redo this section with a better map
  }

  void CRTChannelMapAlg::Uninitialize(){
    //TODO: Release Map
  }
  uint32_t ChannelNumberFromModuleAndPanel(uint32_t const& module,
                                              uint32_t const& panel){
    //TODO: Fill
    return 0;
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
