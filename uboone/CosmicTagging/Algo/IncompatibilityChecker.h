/**
 * \file IncompatibilityChecker.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class IncompatibilityChecker
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Algorithms

    @{*/
#ifndef INCOMPATIBILITYCHECKER_H
#define INCOMPATIBILITYCHECKER_H

#include <iostream>
#include "fhiclcpp/ParameterSet.h"
#include "uboone/LLSelectionTool/OpT0Finder/Base/OpT0FinderTypes.h"
//#include "uboone/LLSelectionTool/OpT0Finder/Base/BaseAlgorithm.h"
//#include "uboone/LLSelectionTool/OpT0Finder/Base/CustomAlgoFactory.h"

namespace flashana{
  
  /**
   \class IncompatibilityChecker
   User defined class IncompatibilityChecker ... these comments are used to generate
   doxygen documentation!
 */

  class IncompatibilityChecker {
    
  public:
    
    /// Default constructor
    IncompatibilityChecker();
    
    /// Default destructor
    ~IncompatibilityChecker(){}

    /// Configure function parameters
    void Configure(fhicl::ParameterSet const& p);

    /// Printd the current configuration
    void PrintConfig();
    
    /// Check if a beam flash (1st arg.) is _NOT_ compatible with an hypothesis flash (2nd arg.)
    bool CheckIncompatibility(const Flash_t &flash, const Flash_t &flash_hypo);

  protected:

    double _sigmaThreshold;
    double _nBinsRequirement;
    bool   _useFlashPosition;
  };
}

#endif
/** @} */ // end of doxygen group 

