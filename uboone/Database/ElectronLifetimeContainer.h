/**
 * \file ElectronLifetime.h
 *
 * \ingroup IOVData
 * 
 * \brief Class def header for a class ElectronLifetime
 *
 * @author eberly@slac.stanford.edu
 */

#ifndef ELECTRONLIFETIMECONTAINER_H
#define ELECTRONLIFETIMECONTAINER_H 1

#include "larevt/CalibrationDBI/IOVData/ChData.h"

namespace lariov {
  /**
     \class ElectronLifetime
  */
  class ElectronLifetimeContainer : public ChData {
    
    public:
    
      /// Constructor
      ElectronLifetimeContainer(unsigned int ch) : ChData(ch) {}
      
      /// Default destructor
      ~ElectronLifetimeContainer() {}
            
      float ExpOffset()       const { return fExpOffset; }
      float TimeConstant()    const { return fTimeConstant; }
      float ExpOffsetErr()    const { return fExpOffsetErr; }
      float TimeConstantErr() const { return fTimeConstantErr; }
      
      void SetExpOffset(float val)       { fExpOffset       = val; }
      void SetTimeConstant(float val)    { fTimeConstant    = val; }
      void SetExpOffsetErr(float val)    { fExpOffsetErr    = val; }
      void SetTimeConstantErr(float val) { fTimeConstantErr = val; }
      
    private:
    
      float fExpOffset;
      float fTimeConstant;
      float fExpOffsetErr;
      float fTimeConstantErr;
      
  }; // end class
} // end namespace lariov

#endif
 
