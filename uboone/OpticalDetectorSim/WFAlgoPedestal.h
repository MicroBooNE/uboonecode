/**
 * \file WFAlgoPedestal.h
 *
 * \ingroup OpticalDetector
 * 
 * \brief Class def header for a class WFAlgoPedestal
 *
 * @author kazuhiro
 */

/** \addtogroup OpticalDetector

    @{*/
#ifndef WFALGOPEDESTAL_H
#define WFALGOPEDESTAL_H

#include "RandomServer.h"
#include "WFAlgoBase.h"

namespace opdet {
  /**
     \class WFAlgoPedestal
     Pedestal class for WF generation algorithms for optical detector
  */
  class WFAlgoPedestal : public WFAlgoBase{
    
  public:
    
    /// Default constructor
    WFAlgoPedestal();
    
    /// Default destructor
    virtual ~WFAlgoPedestal(){}

    /**
       Core function: the algorithm is supposed to add relevant singal/noise ADC
       values to the input vector std::vector<float> wf. The tick=0 timing is 
       provided by detinfo::ElecClock start_time input variable which can also used
       to retrieve G4 time offset to electronics clock counting.
     */
    virtual void Process(std::vector<float> &wf,
			 const ::detinfo::ElecClock &start_time);

    void SetPedestal(double mean, double sigma);

  protected:

    double fPedMean;

    double fPedSigma;

  };
}

#endif
/** @} */ // end of doxygen group 

