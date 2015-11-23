/**
 * \file WFAlgoDigitizedSPE.h
 *
 * \ingroup OpticalDetector
 * 
 * \brief Class def header for a class WFAlgoDigitizedSPE
 *
 * @author kazuhiro
 */

/** \addtogroup OpticalDetector

    @{*/
#ifndef WFALGODIGITIZEDSPE_H
#define WFALGODIGITIZEDSPE_H

#include "WFAlgoSPEBase.h"

namespace opdet {

  /**
     \class WFAlgoDigitizedSPE
     User defined class WFAlgoDigitizedSPE ... these comments are used to generate
     doxygen documentation!
  */
  class WFAlgoDigitizedSPE : public WFAlgoSPEBase{
    
  public:
    
    /// Default constructor
    WFAlgoDigitizedSPE();
    
    /// Default destructor
    virtual ~WFAlgoDigitizedSPE(){}

    /// Function to reset all attributes of the instance
    virtual void Reset();

    /**
       Core function: the algorithm is supposed to add relevant singal/noise ADC
       values to the input vector std::vector<float> wf. The tick=0 timing is 
       provided by detinfo::ElecClock start_time input variable which can also used
       to retrieve G4 time offset to electronics clock counting.
     */
    virtual void Process(std::vector<float> &wf,
			 const ::detinfo::ElecClock &start_time);

    /// Function to get SPE vector the instance currently holds
    const std::vector<float>& GetSPE() const { return fSPE; }

    /**
       Function to set SPE waveform. The second argument is detinfo::ElecClock which
       period should specify waveform tick size and timing specfies photon time
       w.r.t. start of the waveform.
     */
    void SetSPE( const std::vector<float> &wf,
		 const ::detinfo::ElecClock &time_info);

  protected:

    /// SPE waveform
    std::vector<float> fSPE;

    /// SPE waveform timing information (tick period & signal timing)
    ::detinfo::ElecClock fSPETime;

  };
}

#endif
/** @} */ // end of doxygen group 

