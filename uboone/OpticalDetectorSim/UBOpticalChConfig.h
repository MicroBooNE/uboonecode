/**
 * \file UBOpticalChConfig.h
 *
 * \ingroup OpticalDetectorSim
 * 
 * \brief Class def header for a class UBOpticalChConfig
 *
 * @author kazuhiro
 */

/** \addtogroup OpticalDetectorSim

    @{*/
#ifndef UBOPTICALCHCONFIG_H
#define UBOPTICALCHCONFIG_H

#include "Geometry/Geometry.h"
#include "SimpleChConfig.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

namespace opdet {
  /**
     \class ChConfig
     User defined class ChConfig ... these comments are used to generate
     doxygen documentation!
  */
  class UBOpticalChConfig : public SimpleChConfig{
    
  public:
    
    /// Default constructor
    UBOpticalChConfig(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    
    /// Default destructor
    virtual ~UBOpticalChConfig(){};

    void reconfigure(fhicl::ParameterSet const& pset);
    void doInitialization();

  private:
    fhicl::ParameterSet _pset;

  };

}


DECLARE_ART_SERVICE(opdet::UBOpticalChConfig, LEGACY)

#endif
/** @} */ // end of doxygen group 

