////////////////////////////////////////////////////////////////////////////////
/// \file UBooNEGeometryHelper.h
/// \brief Geometry helper service for UBooNE geometries. 
/// 
/// Handles UBooNE-specific information for the generic Geometry service
/// within LArSoft. Derived from the ExptGeoHelperInterface class
///
/// \verion $Id
/// \author rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

#ifndef UBooNE_ExptGeoHelperInterface_h
#define UBooNE_ExptGeoHelperInterface_h

#include "larcore/Geometry/ExptGeoHelperInterface.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/AuxDetGeo.h"

#include <memory>
#include <vector>

// Forward declarations
//
class TString;

namespace geo
{
  class ChannelMapAlg;
}

// Declaration
//
namespace uboone
{
  class UBooNEGeometryHelper : public geo::ExptGeoHelperInterface
  {
  public:
  
    UBooNEGeometryHelper( fhicl::ParameterSet const & pset, art::ActivityRegistry &reg );
    ~UBooNEGeometryHelper() throw();

    // Public interface for ExptGeoHelperInterface (for reference purposes)
    //
    // Configure and initialize the channel map.
    //
    // void  ConfigureChannelMapAlg( const TString & detectorName, 
    //                               fhicl::ParameterSet const & sortingParam,
    //                               std::vector<geo::CryostatGeo*> & c,
    //				     std::vector<geo::AuxDetGeo*>   & ad );
    //
    // Returns null pointer if the initialization failed
    // NOTE:  the sub-class owns the ChannelMapAlg object
    //
    // std::shared_ptr<const geo::ChannelMapAlg> & GetChannelMapAlg() const;
  
  private:

    void doConfigureChannelMapAlg( fhicl::ParameterSet const & sortingParameters, geo::GeometryCore* geom ) override;
    ChannelMapAlgPtr_t doGetChannelMapAlg() const override;
    
    fhicl::ParameterSet fPset; ///< copy of configuration parameter set
    //art::ActivityRegistry fReg; ///< copy of activity registry
    std::shared_ptr<geo::ChannelMapAlg> fChannelMap;
  
  };

}
DECLARE_ART_SERVICE_INTERFACE_IMPL(uboone::UBooNEGeometryHelper, geo::ExptGeoHelperInterface, LEGACY)

#endif // UBooNE_ExptGeoHelperInterface_h
