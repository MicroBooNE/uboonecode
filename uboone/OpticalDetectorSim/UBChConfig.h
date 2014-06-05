/**
 * \file UBChConfig.h
 *
 * \ingroup OpticalDetectorSim
 * 
 * \brief Class def header for a class UBChConfig
 *
 * @author kazuhiro
 */

/** \addtogroup OpticalDetectorSim

    @{*/
#ifndef UBCHCONFIG_H
#define UBCHCONFIG_H

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
  class UBChConfig : public SimpleChConfig{
    
  public:
    
    /// Default constructor
    UBChConfig(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    
    /// Default destructor
    virtual ~UBChConfig(){};

    void reconfigure(fhicl::ParameterSet const& pset);

  };

}


DECLARE_ART_SERVICE(opdet::UBChConfig, LEGACY)

#endif
/** @} */ // end of doxygen group 

