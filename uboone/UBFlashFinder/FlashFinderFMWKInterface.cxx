#ifndef __FLASHFINDERFMWKINTERFACE_CXX__
#define __FLASHFINDERFMWKINTERFACE_CXX__

//#include "cetlib/exception.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/OpDetGeo.h"
#include "uboone/Geometry/UBOpReadoutMap.h"


#include "FlashFinderFMWKInterface.h"
namespace pmtana {

  size_t NOpDets() {
    ::art::ServiceHandle<geo::Geometry> geo;
    return geo->NOpDets();
  }

  size_t OpDetFromOpChannel(size_t opch) {
    ::art::ServiceHandle<geo::Geometry> geo;
    return geo->OpDetFromOpChannel(opch);
  }

  void OpDetCenterFromOpChannel(size_t opch, double *xyz) {
    ::art::ServiceHandle<geo::Geometry> geo;
    geo->OpDetGeoFromOpChannel(opch).GetCenter(xyz); 
  }

}
#endif
