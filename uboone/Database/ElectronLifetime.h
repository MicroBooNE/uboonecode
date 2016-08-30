/**
 * \file ElectronLifetime.h
 *
 * \ingroup IOVData
 * 
 * \brief Class def header for a class ElectronLifetime
 *
 * @author eberly@slac.stanford.edu
 */

/** \addtogroup IOVData

    @{*/
#ifndef IOVDATA_ELECTRONLIFETIME_H
#define IOVDATA_ELECTRONLIFETIME_H 1

#include "larevt/CalibrationDBI/IOVData/ChData.h"

namespace lariov {
  /**
     \class ElectronLifetime
  */
  class ElectronLifetime : public ChData {
    
    public:
    
      /// Constructor
      ElectronLifetime(unsigned int ch) : ChData(ch) {}
      
      /// Default destructor
      ~ElectronLifetime() {}
            
      float ExpPar()    const { return fExpPar; }
      float ConstPar()     const { return fConstPar; }
      
      void SetExpPar(float expPar)       { fExpPar    = expPar; }
      void SetConstPar(float constPar)   { fConstPar  = constPar; }
      
    private:
    
      float fExpPar;
      float fConstPar;
      
  }; // end class
} // end namespace lariov

#endif
/** @} */ // end of doxygen group 
