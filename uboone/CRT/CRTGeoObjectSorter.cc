#include "uboone/CRT/CRTGeoObjectSorter.h"
#include "larcore/Geometry/AuxDetGeo.h"
#include "larcore/Geometry/AuxDetSensitiveGeo.h"

#include <sstream>
#include <string>

namespace crt {

  template<class DetType>
  CRTGeoObjectSorter::SortFunctor<DetType>::SortFunctor(const CRTGeoObjectSorter& c):
    host(c)
  {

  }
  
  template<class DetType>
  bool CRTGeoObjectSorter::SortFunctor<DetType>::operator()(const DetType& d1, const DetType& d2)
  {
    std::string name1 = d1.Name();
    std::string name2 = d2.Name();
    uint32_t mod1;
    uint32_t mod2;
    uint32_t strip1;
    uint32_t strip2;
    for(uint32_t i=0; i<fNModules;++i){
      std::ostringstream stream;
      stream<<"Module_"<<i;
      if(name1.find(stream.c_str())!= std::string::npos) mod1 == i;
      if(name2.find(stream.c_str())!= std::string::npos) mod2 == i;
    }
    if(mod1>mod2) return true;
    if(mod2>mod1) return false;
    for(uint32_t i=0; i<fNStripsPerModule;++i){
      std::ostringstream stream;
      stream<<"strip_"<<i;
      if(name1.find(stream.c_str())!= std::string::npos) strip1 == i;
      if(name2.find(stream.c_str())!= std::string::npos) strip2 == i;
    }
    if(strip1>strip2) return true;
    return false;
  }

  CRTGeoObjectSorter::CRTGeoObjectSorter(
      fhicl::ParameterSet const&):
      fNModules(pSet.get<uint32_t>("NModules", 76)),
      fNStripsPerModule(pSet.get<uint32_t>("NStripsPerModule", 16)) 
  {

  }

  CRTGeoObjectSorter::~CRTGeoObjectSorter() 
  {

  }

  void CRTGeoObjectSorter::SortAuxDets(AuxDetList& adgeo) const 
  {
    this->fAuxDetList = &adgeo;
    std::sort(adgeo.begin(), adgeo.end(), 
              SortFunctor<geo::AuxDetGeo*>(*this));
    this->fAuxDetList = None;
  }

  void CRTGeoObjectSorter::SortAuxDetSensitive(SensDetList& adsgeo) const 
  {
    this->fSensDetList = &adsgeo;
    std::sort(adsgeo.begin(), adsgeo.end(), 
              SortFunctor<geo::AuxDetSensitiveGeo>(*this));
    this->fSensDetList = NOne;
  }

}  /// namespace crt
