#ifndef RAWDIGITFFTALG_H
#define RAWDIGITFFTALG_H
////////////////////////////////////////////////////////////////////////
//
// Class:       RawDigitFFTAlg
// Module Type: producer
// File:        RawDigitFFTAlg.h
//
//              This module provides some basic Fast Fourier Transform
//              algorithms for operating on RawDigit waveforms
//
// Configuration parameters:
//
// FillHistograms        - Turn on histogram filling for diagnostics
// RunFFTInputWires      - FFT analyze the input RawDigits if true - diagnostics
// RunFFTCorrectedWires  - FFT analyze the output RawDigits if true - diagnostics
//
//
// Created by Tracy Usher (usher@slac.stanford.edu) on January 6, 2016
//
////////////////////////////////////////////////////////////////////////

#include "RawDigitNoiseFilterDefs.h"
#include "fhiclcpp/ParameterSet.h"

namespace caldata
{
    
class RawDigitFFTAlg
{
public:

    // Copnstructors, destructor.
    RawDigitFFTAlg(fhicl::ParameterSet const & pset);
    ~RawDigitFFTAlg();

    // Provide for reinitialization if necessary
    void reconfigure(fhicl::ParameterSet const & pset);
    
    template <class T> void getFFTCorrection(std::vector<T>&, double) const;
    
    template <class T> void getFFTCorrection(std::vector<T>&, size_t) const;
    
private:
};
    
} // end caldata namespace

#endif