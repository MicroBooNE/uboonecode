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
/// \update notes:  Yun-Tse Tsai (yuntse@slac.stanford.edu), July 17th, 2014
///                 1. Read in field responses from input histograms
///                 2. Allow different sampling rates in the input
///                    field response
///                    NOTE: The InputFieldRespSamplingRate parameter has
///                    NOT implemented for the field response input
///                    as a function (UseFunctionFieldShape)
///                 3. Allow different electron drift velocities from
///                    which the input field responses are obtained
///                 4. Convolute the field and electronic responses,
///                    and then sample the convoluted function with
///                    the nominal sampling rate (util::DetectorProperties).
///                    NOTE: Currently this doesn't include the filter 
///                    function and the deconvolution kernel.
///                    We may want to include them later?
///                 5. Disable fColSignalShaping.SetPeakResponseTime(0.),
///                    so that the peak time in the input field response
///                    is preserved.
///                 6. Somebody needs to unify the units of time (microsec
///                    or nanosec); I'm fainting!
///
/// New function:   void SetResponseSampling();
///
/// Modified functions: void init();
///                     void SetFieldResponse();
///
/// New FCL parameters:
/// DefaultDriftVelocity       - The electron drift velocity used to obtain
///                              the input field response waveforms
/// InputFieldRespSamplingRate - The sampling rate in the input field response
/// UseHistogramFieldShape     - Use the field response from an input histogram,
///                              if both UseFunctionFieldShape and 
///                              UseHistogramFieldShape are false, we will
///                              use the toy field responses (a bipolar square
///                              function for induction planes, a ramp function
///                              for collection planes.)
/// FieldResponseFname         - Name of the file containing the input field 
///                              response histograms
/// FieldResponseHistoName     - Name of the field response histograms,
///                              the format in the code will be 
///                              FieldResponseHistoName_U(V,Y)
///update notes: Jyoti Joshi (jjoshi@bnl.gov), Jan 13, 2015 
//               1. Modification to GetShapingTime function to read in different
//                  shaping time for different planes
//               2. Modification to GetASICGain fucntion to read in different gain 
//                  settings for different planes    
////////////////////////////////////////////////////////////////////////

#ifndef SIGNALSHAPINGSERVICEMICROBOONE_H
#define SIGNALSHAPINGSERVICEMICROBOONE_H

#include <vector>
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "Utilities/SignalShaping.h"
#include "TF1.h"
#include "TH1D.h"

// LArSoft include
#include "Geometry/Geometry.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Utilities/TimeService.h"


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

    //double GetASICGain() { return fASICGainInMVPerFC; }
    
    std::vector<double> GetNoiseFactInd() { return fNoiseFactInd; }
    std::vector<double> GetNoiseFactColl() { return fNoiseFactColl; }
    

    //double GetShapingTime() { return fShapeTimeConst.at(1); }

    double GetASICGain(unsigned int const channel) const;
    double GetShapingTime(unsigned int const channel) const; 

    const util::SignalShaping& SignalShaping(unsigned int channel) const;

    int FieldResponseTOffset(unsigned int const channel) const;

    // Do convolution calcution (for simulation).

    template <class T> void Convolute(unsigned int channel, std::vector<T>& func) const;

    // Do deconvolution calcution (for reconstruction).

    template <class T> void Deconvolute(unsigned int channel, std::vector<T>& func) const;
    
    void SetDecon(int fftsize);
  private:

    // Private configuration methods.

    // Post-constructor initialization.

    void init() const{const_cast<SignalShapingServiceMicroBooNE*>(this)->init();}
    void init();

    // Calculate response functions.
    // Copied from SimWireMicroBooNE.

    void SetFieldResponse();

    void SetElectResponse(double shapingtime, double gain);  //changed to read different peaking time for different planes

    // Calculate filter functions.

    void SetFilters();

    // Attributes.

    bool fInit;               ///< Initialization flag.

    // Sample the response function, including a configurable
    // drift velocity of electrons

    void SetResponseSampling();

    

    // Fcl parameters.
    double fADCPerPCAtLowestASICGain; ///< Pulse amplitude gain for a 1 pc charge impulse after convoluting it the with field and electronics response with the lowest ASIC gain setting of 4.7 mV/fC

    //double fASICGainInMVPerFC;                  ///< Cold electronics ASIC gain setting in mV/fC
    std::vector<double> fASICGainInMVPerFC;                  ///< Cold electronics ASIC gain setting in mV/fC
    std::vector<double> fNoiseFactColl;         ///< RMS Noise in ADCs for lowest Gain Setting
    std::vector<double> fNoiseFactInd;          ///< RMS Noise in ADCs for lowest Gain Setting

    std::vector<double> fDefaultDriftVelocity;  ///< Default drift velocity of electrons in cm/usec
    std::vector<double> fFieldResponseTOffset;  ///< Time offset for field response in ns
    int fNFieldBins;         			///< number of bins for field response
    double fInputFieldRespSamplingPeriod;       ///< Sampling period in the input field response.
    double fCol3DCorrection; 			///< correction factor to account for 3D path of 
						///< electrons thru wires
    double fInd3DCorrection;  			///< correction factor to account for 3D path of 
						///< electrons thru wires
    double fColFieldRespAmp;  			///< amplitude of response to field 
    double fIndUFieldRespAmp;  			///< amplitude of response to field in U plane
    double fIndVFieldRespAmp;  			///< amplitude of response to field in V plane
    std::vector<double> fShapeTimeConst;  	///< time constants for exponential shaping
    std::vector<int> fDeconvPol;                ///< switch for DeconvKernel normalization sign (+ -> max pos ADC, - -> max neg ADC). Entry 0,1,2 = U,V,Y plane settings
    TF1* fColFilterFunc;      			///< Parameterized collection filter function.
    TF1* fIndUFilterFunc;      			///< Parameterized induction filter function for U plane.
    TF1* fIndVFilterFunc;      			///< Parameterized induction filter function for V plane

    
    bool fUseFunctionFieldShape;   		///< Flag that allows to use a parameterized field response instead of the hardcoded version
    bool fUseHistogramFieldShape;               ///< Flag that turns on field response shapes from histograms
    bool fGetFilterFromHisto;   		///< Flag that allows to use a filter function from a histogram instead of the functional dependency
    TF1* fColFieldFunc;      			///< Parameterized collection field shape function.
    TF1* fIndUFieldFunc;      			///< Parameterized induction field shape function for U plane.
    TF1* fIndVFieldFunc;      			///< Parameterized induction field shape function for V plane.

    TH1F *fFieldResponseHist[3];                ///< Histogram used to hold the field response, hardcoded for the time being

    TH1D *fFilterHist[3];    			///< Histogram used to hold the collection filter, hardcoded for the time being
    
    // Following attributes hold the convolution and deconvolution kernels

    util::SignalShaping fIndUSignalShaping;
    util::SignalShaping fIndVSignalShaping;
    util::SignalShaping fColSignalShaping;
    // Field response.

    std::vector<double> fIndUFieldResponse;
    std::vector<double> fIndVFieldResponse;
    std::vector<double> fColFieldResponse;

    // Electronics response.

    std::vector<double> fElectResponse;

    // Filters.

    std::vector<TComplex> fIndUFilter;
    std::vector<TComplex> fIndVFilter;
    std::vector<TComplex> fColFilter;
  };
}


//----------------------------------------------------------------------
// Do convolution.
template <class T> void util::SignalShapingServiceMicroBooNE::Convolute(unsigned int channel, std::vector<T>& func) const
{
  SignalShaping(channel).Convolute(func);
  //negative number
  int time_offset = FieldResponseTOffset(channel);
  
  std::vector<T> temp;
  if (time_offset <=0){
    temp.assign(func.begin(),func.begin()-time_offset);
    func.erase(func.begin(),func.begin()-time_offset);
    func.insert(func.end(),temp.begin(),temp.end());
  }else{
    temp.assign(func.end()-time_offset,func.end());
    func.erase(func.end()-time_offset,func.end());
    func.insert(func.begin(),temp.begin(),temp.end());
  }
  
}


//----------------------------------------------------------------------
// Do deconvolution.
template <class T> void util::SignalShapingServiceMicroBooNE::Deconvolute(unsigned int channel, std::vector<T>& func) const
{
  SignalShaping(channel).Deconvolute(func);
  int time_offset = FieldResponseTOffset(channel);
  
  std::vector<T> temp;
  if (time_offset <=0){
    temp.assign(func.end()+time_offset,func.end());
    func.erase(func.end()+time_offset,func.end());
    func.insert(func.begin(),temp.begin(),temp.end());
  }else{
    temp.assign(func.begin(),func.begin()+time_offset);
    func.erase(func.begin(),func.begin()+time_offset);
    func.insert(func.end(),temp.begin(),temp.end());

    
  }

}

DECLARE_ART_SERVICE(util::SignalShapingServiceMicroBooNE, LEGACY)
#endif
