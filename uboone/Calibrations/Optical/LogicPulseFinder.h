/**
 * \file LogicPulseFinder.h
 *
 * \ingroup PulseReco
 * 
 * \brief Class def header for a class LogicPulseFinder
 *
 * @author jarrett
 */

/** \addtogroup PulseReco

    @{*/
#ifndef LOGICPULSEFINDER_H
#define LOGICPULSEFINDER_H

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>
#include <utility>

/**
   \class LogicPulseFinder
   User defined class LogicPulseFinder ... these comments are used to generate
   doxygen documentation!
 */

template <class T>
class LogicPulseFinder {

public:

  /// Default constructor
  LogicPulseFinder(){}

  /// Default destructor
  ~LogicPulseFinder(){}

//Simple TTL pulse finder, scans through a vector which may contain a series of TTL pulses.
//Uses a simple global maximum as a reference for all potential peaks. Works since for a good
//TTL waveform the peaks should all be very similar in structure. Simply scans through looking for
//rising points at 50% of max and calls them the pulse start time. Fills an output vector
//with these points as well as outputting the number of peaks found.
  std::vector<unsigned short> Get_TTL_Starts(std::vector<T> wfm, T threshold);

  std::map<std::string,float> TTL_Health(std::vector<T> wfm);
};

#endif
/** @} */ // end of doxygen group 

