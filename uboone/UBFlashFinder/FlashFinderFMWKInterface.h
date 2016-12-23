#ifndef __FLASHFINDERFMWKINTERFACE_H__
#define __FLASHFINDERFMWKINTERFACE_H__

//#include "FhiclLite/ConfigManager.h"
#include "cetlib/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/Geometry.h"
#include <stdlib.h>
namespace pmtana {

  //typedef ::fcllite::PSet Config_t;
  typedef fhicl::ParameterSet Config_t;

  size_t NOpDets();

  size_t OpDetFromOpChannel(size_t opch);

  void OpDetCenterFromOpChannel(size_t opch, double *xyz);
  
}
#endif
