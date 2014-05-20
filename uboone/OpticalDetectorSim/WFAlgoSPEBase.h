/**
 * \file WFAlgoSPEBase.h
 *
 * \ingroup OpticalDetector
 * 
 * \brief Class def header for a class WFAlgoSPEBase
 *
 * @author kazuhiro
 */

/** \addtogroup OpticalDetectorSim

    @{*/
#ifndef WFALGOSPEBASE_H
#define WFALGOSPEBASE_H

#include "TString.h"
#include "WFAlgoBase.h"
#include "RandomServer.h"
namespace opdet {

  /**
     \class WFAlgoSPEBase
     User defined class WFAlgoSPEBase ... these comments are used to generate
     doxygen documentation!
  */
  class WFAlgoSPEBase : public WFAlgoBase{
    
  public:
    
    /// Default constructor
    WFAlgoSPEBase();
    
    /// Default destructor
    virtual ~WFAlgoSPEBase(){}

    /// Function to reset all attributes of the instance
    virtual void Reset();

    /// Function to set gain and fractional spread
    void SetGain(double gain, double fractional_spread)
    { fGain = gain; fGainSigma = fractional_spread; }

    /// Function to set T0 and absolute spread
    void SetT0(double t0, double t0_spread)
    { fT0 = t0; fT0Sigma = t0_spread; }

    /// Function to add G4 photon @ specific G4 time
    void AddPhoton(double g4time) { fPhotonTime.push_back(g4time); }

    /// Function to set G4 photons (in G4 time)
    void SetPhotons(const std::vector<double> &g4time);
    
    /// Function to enable spread of parameters
    void EnableSpread(bool doit=true) { fEnableSpread = doit; }

    const std::vector<double>& GetPhotonTime() const { return fPhotonTime;}

  protected:

    /// Boolean to enable spread
    bool fEnableSpread;

    /// Amplitude of SPE pulse
    double fGain;

    /// Spread of amplitude in fraction
    double fGainSigma;

    /// T0 of SPE pulse
    double fT0;

    /// T0 spread of SPE pulse
    double fT0Sigma;

    /// G4 photon times at which SPE will be generated/added
    std::vector<double> fPhotonTime;
    
  };
}

#endif
/** @} */ // end of doxygen group 

