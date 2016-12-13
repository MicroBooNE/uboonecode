#ifndef CRTGeoObjectSorter_hh_
#define CRTGeoObjectSorter_hh_

#include "larcore/Geometry/AuxDetGeoObjectSorter.h"
#include <vector>

namespace crt {

  class CRTGeoObjectSorter : public geo::AuxDetGeoObjectSorter {
    
    typedef std::vector<geo::AuxDetGeo*> AuxDetList;
    typedef std::vector<geo::AuxDetSensitiveGeo*> SensDetList;

    uint32_t fNModules;
    uint32_t fNStripsPerModule;

    template<class DetType>
    struct SortFunctor{
      SortFunctor(const CRTGeoObjectSorter& c);
      bool operator()(const DetType& d1, const DetType& d2);
      CRTGeoObjectSorter& host;
    };

  public:

    CRTGeoObjectSorter(fhicl::ParameterSet const& p);

    ~CRTGeoObjectSorter();

    void SortAuxDets (AuxDetList& adgeo) const;

    void SortAuxDetSensitive(SensDetList& adsgeo) const;

  };
}

#endif  // CRTGeoObjectSorter_hh_