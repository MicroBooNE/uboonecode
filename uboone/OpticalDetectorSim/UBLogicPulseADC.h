/**
 * \file UBLogicPulseADC.h
 *
 * \ingroup OpticalDetectorSim
 * 
 * \brief Class def header for a class UBLogicPulseADC
 *
 * @author kazuhiro
 */

/** \addtogroup OpticalDetectorSim

    @{*/
#ifndef UBLOGICPULSEADC_H
#define UBLOGICPULSEADC_H

#include "UBADCBase.h"

namespace opdet {
  /**
     \class UBLogicPulseADC
     A generator class that handles WF generation algorithm to create an waveform.
  */
  class UBLogicPulseADC : public UBADCBase {
    
  public:
    
    /// Default constructor
    UBLogicPulseADC();
    
    /// Default destructor
    virtual ~UBLogicPulseADC(){}

    /// Function to reset algorithm configuration
    virtual void Reset();

    /// Function to set Logic pulse amplitude
    void SetAmplitude(float amp) { fSPE.SetGain(amp,0); }

    /// Function to set Logic pulse pedestal
    void SetPedestal(float mean, float sigma) { fPED.SetPedestal(mean,sigma); }

    /// Function to add G4 photon in G4 time (ns as that is G4 natural unit)
    void AddPulse(double g4time){ fPulseTime.push_back(g4time); }

    /// Function to set G4 photons in G4 time (ns as that is G4 natural unit)
    void SetPulses(const std::vector<double>& g4time);

    /// Method to generate waveform for a specific channel
    void GenWaveform(std::vector<unsigned short> &logic_wf);

  protected:

    /// Photon time that is injected to the waveform (after QE applied)
    std::vector<double> fPulseTime;

    /// Algorithm to generate SPE waveform
    WFAlgoAnalyticalSPE fSPE;

    /// Algorithm to generate Pedestal
    WFAlgoPedestal      fPED;

  };
}

#endif
/** @} */ // end of doxygen group 

