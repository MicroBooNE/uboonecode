#include "uboone/CRT/CRTChannelMapAlg.hh"

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
    // TODO: Redo this section with a better map
  }

  void CRTChannelMapAlg::Uninitialize(){
    //TODO: Release Map
  }

}
