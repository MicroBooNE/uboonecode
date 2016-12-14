<<<<<<< HEAD
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

=======
>>>>>>> 94a7a81016e00031d5dd937416fbb3da0a440206
#ifndef CRTGeoObjectSorter_hh_
#define CRTGeoObjectSorter_hh_

#include "larcore/Geometry/AuxDetGeoObjectSorter.h"
#include <vector>

namespace crt {

  class CRTGeoObjectSorter : public geo::AuxDetGeoObjectSorter {
    
<<<<<<< HEAD
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
    
  protected:

    /// fNmodules getter
    inline uint32_t GetNModules() const {return this->fNModules;}

    /// fNStripsPerModule getter
    inline uint32_t GetNStripsPerModule() const {return this->fNStripsPerModule;}

  public:

    /// Artistic c'tor
    CRTGeoObjectSorter(fhicl::ParameterSet const& p);

    /// Default d'tor
    ~CRTGeoObjectSorter();

    /// Default sorter
    void SortAuxDets (AuxDetList& adgeo) const;

    /// Blanked because there is 1 SV per AD
    void SortAuxDetSensitive(SensDetList& adsgeo) const {}
=======
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
>>>>>>> 94a7a81016e00031d5dd937416fbb3da0a440206

  };
}

#endif  // CRTGeoObjectSorter_hh_