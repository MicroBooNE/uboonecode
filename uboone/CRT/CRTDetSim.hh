/**
 * \class CRTDetSim
 *
 * \ingroup crt
 *
 * \brief Provides CRTData from simulations
 *
 * Converts IDEs from largeant (or whichever producer) to 
 * CRTData. This is meant to mimic the physical detector as much as
 * possible.
 *
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
**/

#ifndef CRTDetSim_HH_
#define CRTDetSim_HH_


#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"

#include <string>

#include <string>

namespace crt{
  class CRTDetSim :  public art:: EDProducer{
    // This is a basic ADC threshold.
    uint32_t fThreshold;
    /// This is the factor to go from IDE to ADC
    float fConversionFactor;
    /// Precision factor when converting to T1
    float fT1Precision;
    /// Name of the producer of the IDEs
    std::string fProducerName;
  public:

    /// Default ctor
    CRTDetSim(const fhicl::ParameterSet&);

    /// Default dtor
    ~CRTDetSim();

    /// art::EDProducer::produce implementation
    virtual void produce (art::Event&);

  };
}


#endif  //CRTDetSim_HH_
