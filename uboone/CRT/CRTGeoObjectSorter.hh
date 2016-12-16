/**
 * \class CRTGeoObjectSorter
 *
 * \ingroup crt
 *
 * \brief Sorts AuDetGeos by module number and strip number
 *
 * In order to map the CRT panels, they need to be arranged by
 * module and strip number. This is accomplished by sorting on this
 * pattern.
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

#ifndef CRTGeoObjectSorter_hh_
#define CRTGeoObjectSorter_hh_

#include "larcore/Geometry/AuxDetGeoObjectSorter.h"
#include <vector>

namespace crt {

  class CRTGeoObjectSorter : public geo::AuxDetGeoObjectSorter {
    
    /// Convenience typedef for templating on  geometry type
    typedef std::vector<geo::AuxDetGeo*> AuxDetList;
    typedef std::vector<geo::AuxDetSensitiveGeo*> SensDetList;

    /// Number of modules in a geometry
    uint32_t fNModules;
    /// Number of Strips per module
    uint32_t fNStripsPerModule;

    /// Helper struct for callbacks during STL sorting
    template<class DetType>
    struct SortFunctor{
      SortFunctor(const CRTGeoObjectSorter& c);
      bool operator()(DetType* d1, DetType* d2);
      const CRTGeoObjectSorter& host;
    };

  public:
    /// fNmodules getter
    inline uint32_t GetNModules() const {return this->fNModules;}

    /// fNStripsPerModule getter
    inline uint32_t GetNStripsPerModule() const {return this->fNStripsPerModule;}

    /// Artistic c'tor
    CRTGeoObjectSorter(fhicl::ParameterSet const& p);

    /// Default d'tor
    ~CRTGeoObjectSorter();

    /// Default sorter
    void SortAuxDets (AuxDetList& adgeo) const;

    /// Blanked because there is 1 SV per AD
    void SortAuxDetSensitive(SensDetList& adsgeo) const {}

  };
}

#endif  // CRTGeoObjectSorter_hh_