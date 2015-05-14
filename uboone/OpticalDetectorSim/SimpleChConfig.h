/**
 * \file SimpleChConfig.h
 *
 * \ingroup OpticalDetector
 * 
 * \brief Class def header for a class SimpleChConfig
 *
 * @author kazuhiro
 */

/** \addtogroup OpticalDetector

    @{*/
#ifndef SIMPLECHCONFIG_H
#define SIMPLECHCONFIG_H

#include <vector>
#include "UBOpticalException.h"
#include "UBOpticalConstants.h"
namespace opdet {
  /**
     \class SimpleChConfig
     User defined class SimpleChConfig ... these comments are used to generate
     doxygen documentation!
  */
  class SimpleChConfig{
    
  protected:
    
    /// Default constructor
    SimpleChConfig() : fParams(kChConfigTypeMax,std::vector<float>()),
		       fDefault(kChConfigTypeMax,0)
    {
      size_t n_channels = kLogicStartChannel;
      
      for(size_t i=0; i<kChConfigTypeMax; ++i)
	fParams.at(i).resize(n_channels,0);
      
      fDefault.at(kPedestalMean)   = 2048;
      fDefault.at(kPedestalSpread) = 0.3;
      fDefault.at(kQE)             = 0.01;
      //fDefault.at(kHighGain)        = 20;
      //fDefault.at(kLowGain)        = 4;
      fDefault.at(kGain)           = 20.0; 
      fDefault.at(kPMTGain)        = 1.0;
      fDefault.at(kSplitterGain)   = 20.0;
      fDefault.at(kGainSpread)     = 0.05;
      fDefault.at(kT0)             = 0;
      fDefault.at(kT0Spread)       = 0;
      fDefault.at(kDarkRate)       = 1.e-5;
      
      for(size_t i=0; i<n_channels; ++i) {
	fParams.at(kPedestalMean).at(i)   = fDefault.at(kPedestalMean);
	fParams.at(kPedestalSpread).at(i) = fDefault.at(kPedestalSpread);
	fParams.at(kQE).at(i)             = fDefault.at(kQE);
	fParams.at(kGain).at(i)           = fDefault.at(kGain);
	fParams.at(kPMTGain).at(i)        = fDefault.at(kPMTGain);
	fParams.at(kSplitterGain).at(i)   = fDefault.at(kSplitterGain);
	fParams.at(kGainSpread).at(i)     = fDefault.at(kGainSpread);
	fParams.at(kT0).at(i)             = fDefault.at(kT0);
	fParams.at(kT0Spread).at(i)       = fDefault.at(kT0Spread);
	fParams.at(kDarkRate).at(i)       = fDefault.at(kDarkRate);
      }

    }
    
    /// Default destructor
    virtual ~SimpleChConfig(){};

  public:

    float GetParameter(const ChConfigType_t type, 
		       const unsigned short ch) const
    {
      if(ch == kINVALID_CHANNEL) return fDefault.at(type);
      if(ch >= fParams.at(type).size()) 
	throw UBOpticalException("Invalid channel number provided!");
      return fParams.at(type).at(ch);
    }

    const std::vector<float>& GetParameter(const ChConfigType_t type) const
    { return fParams.at(type); }


  protected:
    
    std::vector<std::vector<float> > fParams;

    std::vector<float>  fDefault;
    
  };
}

#endif
/** @} */ // end of doxygen group 

