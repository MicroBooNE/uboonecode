/**
 * \file WFAlgoAnalyticalSPE.h
 *
 * \ingroup OpticalDetector
 * 
 * \brief Class def header for a class WFAlgoAnalyticalSPE
 *
 * @author kazuhiro
 */

/** \addtogroup OpticalDetectorSim

    @{*/
#ifndef WFALGOANALYTICALSPE_H
#define WFALGOANALYTICALSPE_H

#include <cmath>
#include "WFAlgoSPEBase.h"

namespace opdet {

  /**
     \class WFAlgoAnalyticalSPE
     User defined class WFAlgoAnalyticalSPE ... these comments are used to generate
     doxygen documentation!
  */
  class WFAlgoAnalyticalSPE : public WFAlgoSPEBase{
    
  public:
    
    /// Default constructor
    WFAlgoAnalyticalSPE();
    
    /// Default destructor
    virtual ~WFAlgoAnalyticalSPE(){}

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

  protected:

    /// Function to evaluate SPE formula @ time = x
    double EvaluateSPE(const double x) const;

  };
}

#endif
/** @} */ // end of doxygen group 

