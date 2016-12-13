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
#include "uboone/LLSelectionTool/OpT0Finder/Base/BaseAlgorithm.h"
#include "uboone/LLSelectionTool/OpT0Finder/Base/CustomAlgoFactory.h"

namespace flashana{
/**
   \class IncompatibilityChecker
   User defined class IncompatibilityChecker ... these comments are used to generate
   doxygen documentation!
 */

  class IncompatibilityChecker : public flashana::BaseAlgorithm {
    
  public:
    
    /// Default constructor
    IncompatibilityChecker(const std::string name="IncompatibilityChecker");
    
    /// Default destructor
    ~IncompatibilityChecker(){}

    ///
    bool CheckIncompatibility(const Flash_t &flash, const Flash_t &flash_hypo);


  protected:

    void _Configure_(const Config_t &pset);
    
    double _gap;
    double _sigmaThreshold;
    double _nBinsRequirement;
    bool   _useFlashPosition;
  };
  
  /**
     \class flashana::IncompatibilityCheckerFactory
  */
  class IncompatibilityCheckerFactory : public CustomAlgoFactoryBase {
  public:
    /// ctor
    IncompatibilityCheckerFactory() { CustomAlgoFactory::get().add_factory("IncompatibilityChecker",this); }
    /// dtor
    ~IncompatibilityCheckerFactory() {}
    /// creation method
    BaseAlgorithm* create(const std::string instance_name) { return new IncompatibilityChecker(instance_name); }
  };
} 

#endif
/** @} */ // end of doxygen group 

