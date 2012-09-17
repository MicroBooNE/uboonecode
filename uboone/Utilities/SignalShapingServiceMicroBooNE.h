///////////////////////////////////////////////////////////////////////
///
/// \file   SignalShapingServiceMicroBooNE.h
///
/// \brief  Service to provide microboone-specific signal shaping for
///         simulation (convolution) and reconstruction (deconvolution).
///
/// \author H. Greenlee 
///
/// This service inherits from SignalShaping and supplies
/// microboone-specific configuration.  It is intended that SimWire and
/// CalWire modules will access this service.
///
/// FCL parameters:
///
/// FieldBins       - Number of bins of field response.
/// Col3DCorrection - 3D path length correction for collection plane.
/// Ind3DCorrection - 3D path length correction for induction plane.
/// ColFieldRespAmp - Collection field response amplitude.
/// IndFieldRespAmp - Induction field response amplitude.
/// ShapeTimeConst  - Time constants for exponential shaping.
/// ColFilter       - Root parameterized collection plane filter function.
/// ColFilterParams - Collection filter function parameters.
/// IndFilter       - Root parameterized induction plane filter function.
/// IndFilterParams - Induction filter function parameters.
///
////////////////////////////////////////////////////////////////////////

#ifndef SIGNALSHAPINGSERVICEMICROBOONE_H
#define SIGNALSHAPINGSERVICEMICROBOONE_H

#include <vector>
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "Utilities/SignalShaping.h"
#include "TF1.h"

namespace util {
  class SignalShapingServiceMicroBooNE {
  public:

    // Constructor, destructor.

    SignalShapingServiceMicroBooNE(const fhicl::ParameterSet& pset,
				   art::ActivityRegistry& reg);
    ~SignalShapingServiceMicroBooNE();

    // Update configuration parameters.

    void reconfigure(const fhicl::ParameterSet& pset);

    // Accessors.

    const util::SignalShaping& SignalShaping(int plane) const;

    // Do convolution calcution (for simulation).

    template <class T> void Convolute(int plane, std::vector<T>& func) const;

    // Do deconvolution calcution (for reconstruction).

    template <class T> void Deconvolute(int plane, std::vector<T>& func) const;

  private:

    // Private configuration methods.

    // Post-constructor initialization.

    void init() const{const_cast<SignalShapingServiceMicroBooNE*>(this)->init();}
    void init();

    // Calculate response functions.
    // Copied from SimWireMicroBooNE.

    void SetFieldResponse();
    void SetElectResponse();

    // Calculate filter functions.

    void SetFilters();

    // Attributes.

    bool fInit;               ///< Initialization flag.

    // Fcl parameters.

    int fNFieldBins;                      ///< number of bins for field response	      
    double fCol3DCorrection;   		  ///< correction factor to account for 3D path of
                              		  ///< electrons thru wires			      
    double fInd3DCorrection;   		  ///< correction factor to account for 3D path of
                              		  ///< electrons thru wires			      
    double fColFieldRespAmp;  		  ///< amplitude of response to field 	      
    double fIndFieldRespAmp;  		  ///< amplitude of response to field             
    std::vector<double> fShapeTimeConst;  ///< time constants for exponential shaping
    TF1* fColFilterFunc;                  ///< Parameterized collection filter function.
    TF1* fIndFilterFunc;      		  ///< Parameterized induction filter function. 

    // Following attributes hold the convolution and deconvolution kernels

    util::SignalShaping fColSignalShaping;
    util::SignalShaping fIndSignalShaping;

    // Field response.

    std::vector<double> fColFieldResponse;
    std::vector<double> fIndFieldResponse;

    // Electronics response.

    std::vector<double> fElectResponse;

    // Filters.

    std::vector<TComplex> fColFilter;
    std::vector<TComplex> fIndFilter;
  };
}

//----------------------------------------------------------------------
// Do convolution.
template <class T> inline void util::SignalShapingServiceMicroBooNE::Convolute(int plane, std::vector<T>& func) const
{
  SignalShaping(plane).Convolute(func);
}


//----------------------------------------------------------------------
// Do deconvolution.
template <class T> inline void util::SignalShapingServiceMicroBooNE::Deconvolute(int plane, std::vector<T>& func) const
{
  SignalShaping(plane).Deconvolute(func);
}


#endif
