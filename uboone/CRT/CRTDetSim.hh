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
    float fTDelayNorm;  //!< Time delay fit: Gaussian normalization
    float fTDelayShift;  //!< Time delay fit: Gaussian x shift
    float fTDelaySigma;  //!< Time delay fit: Gaussian width
    float fTDelayOffset;  //!< Time delay fit: Gaussian baseline offset
    float fTDelayRMSGausNorm;  //!< Time delay RMS fit: Gaussian normalization
    float fTDelayRMSGausShift;  //!< Time delay RMS fit: Gaussian x shift
    float fTDelayRMSGausSigma;  //!< Time delay RMS fit: Gaussian width
    float fTDelayRMSExpNorm;  //!< Time delay RMS fit: Exponential normalization
    float fTDelayRMSExpShift;  //!< Time delay RMS fit: Exponential x shift
    float fTDelayRMSExpScale;  //!< Time delay RMS fit: Exponential scale
    float fNpeScaleNorm;  //!< Npe vs. distance: 1/r^2 scale
    float fNpeScaleShift;  //!< Npe vs. distance: 1/r^2 x shift
    float fQ0;  // Average energy deposited for mips, for charge scaling [GeV]
    float fQPed;  // ADC offset for the single-peak peak mean [ADC]
    float fQSlope;  // Slope in mean ADC / Npe [ADC]
    float fQRMS;  // ADC single-pe spectrum width [ADC]
    float fTResInterpolator;  // Interpolator time resolution [ns]
    float fPropDelay;  // Delay in pulse arrival time [ns/m]
    float fPropDelayError;  // Delay in pulse arrival time, uncertainty [ns/m]
    float fAbsLenEff;  // Effective abs. length for transverse Npe scaling [cm]

    /// Name of the producer of the IDEs
    std::string fProducerName;

    /**
     * Get the channel trigger time relative to the start of the MC event.
     *
     * @param engine The random number generator engine
     * @param clock The clock to count ticks on
     * @param t0 The starting time (which delay is added to)
     * @param npe Number of observed photoelectrons
     * @param r Distance between the energy deposit and strip readout end [mm]
     * @return The channel trigger time [ns]
     */
    double getChannelTriggerTicks(CLHEP::HepRandomEngine* engine,
                                  detinfo::ElecClock& clock,
                                  float t0, float npeMean, float r);


  public:

    /// Default ctor
    CRTDetSim(const fhicl::ParameterSet&);

    /// Default dtor
    ~CRTDetSim();

    /// art::EDProducer::produce implementation
    virtual void produce (art::Event&);

    /// Set up the Configuration Parameters
    void reconfigure(fhicl::ParameterSet const & p) override;

  };
}


#endif  //CRTDetSim_HH_
