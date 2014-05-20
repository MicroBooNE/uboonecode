/**
 * \file UBADCBase.h
 *
 * \ingroup OpticalDetectorSim
 * 
 * \brief Class def header for a class UBADCBase
 *
 * @author kazuhiro
 */

/** \addtogroup OpticalDetectorSim

    @{*/
#ifndef UBADCBase_H
#define UBADCBase_H

#include "UBChConfig.h"
#include "WFAlgoPedestal.h"
#include "WFAlgoAnalyticalSPE.h"
#include "RandomServer.h"

namespace opdet {
  /**
     \class UBADCBase
     A generator class that handles WF generation algorithm to create an waveform.
  */
  class UBADCBase{
    
  public:
    
    /// Default constructor
    UBADCBase();
    
    /// Default destructor
    virtual ~UBADCBase(){}

    /// Function to reset algorithm configuration
    virtual void Reset();

    /// Setter for waveform start time & duration
    void SetTimeInfo(const util::ElecClock &start_freq,
		     double duration);

  protected:

    /// Function to generate dark noise photon timings (called in GenWaveform)
    void GenDarkNoise(double dark_rate,double period);

    /// Method to digitize waveform
    void Digitize(const std::vector<float>& orig,
		  std::vector<unsigned short>& res) const;

  protected:

    /// Time information (digitization frequency & waveform start time)
    ::util::ElecClock fTimeInfo;

    /// Length of waveform
    double fDuration;

  };
}

#endif
/** @} */ // end of doxygen group 

