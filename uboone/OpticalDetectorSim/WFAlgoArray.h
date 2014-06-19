/**
 * \file WFAlgoArray.h
 *
 * \ingroup OpticalDetector
 * 
 * \brief Class def header for a class WFAlgoArray
 *
 * @author kazuhiro
 */

/** \addtogroup OpticalDetector

    @{*/
#ifndef WFALGOARRAY_H
#define WFALGOARRAY_H

#include <TStopwatch.h>
#include "WFAlgoBase.h"

namespace opdet {
  /**
     \class WFAlgoArray
     Array of WFAlgoBase inherit waveform generation algorithms.
     Itself is an algorithm, and can be attached to OpticalWFGenerator
  */
  class WFAlgoArray : public WFAlgoBase{
    
  public:
    
    /// Default constructor
    WFAlgoArray();
    
    /// Default destructor
    virtual ~WFAlgoArray(){}

    /// Adder of an algorithm
    void AddAlgo(WFAlgoBase* algo) {fAlgoArray.push_back(algo);}

    /**
       Core function: the algorithm is supposed to add relevant singal/noise ADC
       values to the input vector std::vector<float> wf. The tick=0 timing is 
       provided by util::ElecClock start_time input variable which can also used
       to retrieve G4 time offset to electronics clock counting.
     */
    virtual void Process(std::vector<float> &wf,
			 const ::util::ElecClock &start_time);
      
    virtual void Reset() { for(auto &algo : fAlgoArray) {algo->Reset();} }

    const std::vector<double>& GetTimeProfile() const { return fTimeArray;}


  protected:

    /// Algorithms to run
    std::vector<opdet::WFAlgoBase*> fAlgoArray;

    /// Stopwatch
    TStopwatch fWatch;

    /// Time profile result
    std::vector<double> fTimeArray;

  };
}

#endif
/** @} */ // end of doxygen group 

