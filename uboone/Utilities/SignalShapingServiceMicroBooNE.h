///////////////////////////////////////////////////////////////////////
///
/// \file   SignalShapingServiceMicroBooNE.h
///
/// \brief  Service to provide microboone-specific signal shaping for
///         simulation (convolution) and reconstruction (deconvolution)./Users/Lsrea/newSim/SignalShapingServiceMicroBooNE.h
///
/// \author H. Greenlee, major mods by L. Rochester
///
/// This service inherits from SignalShaping and supplies
/// microboone-specific configuration.  It is intended that SimWire and
/// CalWire modules will access this service.
///
/// FCL parameters:
///
/// FieldBins       - Number of bins of field response (generated from the histogram).
/// Col3DCorrection - 3D path length correction for collection plane. (not invoked)
/// Ind3DCorrection - 3D path length correction for induction plane.  (not invoked)
/// FieldRespAmpVec - vector of response amplitudes, one for each view
/// ShapeTimeConst  - Time constants for exponential shaping.
/// FilterVec       - vector of filter function parameters, one for each view
/// FilterParamsVec - Vector of filter function parameters.
///
/// \update notes: Leon Rochester (lsrea@slac.stanford.edu, Jan 12, 2015
///                many changes, need to be documented better
///                 1. the three (or n) views are now represented by a vector of views
///                 2. There are separate SignalShaping objects for convolution and
///                    deconvolution
///
///                Yun-Tse Tsai (yuntse@slac.stanford.edu), July 17th, 2014
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

using DoubleVec = std::vector<double>;

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

    //  double GetASICGain()                                  { return fASICGainInMVPerFC; }
    std::vector<DoubleVec> GetNoiseFactVec()                { return fNoiseFactVec; }
    //double GetShapingTime()                                 { return fShapeTimeConst.at(1); }; 
    std::vector<std::vector<size_t> > GetNResponses()       { return fNResponses; }
    std::vector<std::vector<size_t> > GetNActiveResponses() { return fNActiveResponses; }

    std::vector<size_t> GetViewIndex()       { return fViewIndex; }
    //double GetASICGain() { return fASICGainInMVPerFC; }

    //double GetShapingTime() { return fShapeTimeConst.at(1); }

    double GetASICGain(unsigned int const channel) const;
    double GetShapingTime(unsigned int const channel) const; 

    double GetRawNoise(unsigned int const channel) const ;
    double GetDeconNoise(unsigned int const channel) const;

    const std::vector<TComplex>& GetConvKernel(unsigned int channel, unsigned int wire) const;  // M. Mooney
    double Get2DFilterVal(size_t planeNum, size_t freqDimension, double binFrac) const;  // M. Mooney
    double Get2DFilterNorm(size_t planeNum) const;  // M. Mooney

    const util::SignalShaping& SignalShaping(unsigned int channel, unsigned wire = 0, size_t ktype = 0) const;

    int FieldResponseTOffset(unsigned int const channel, size_t ktype) const;

    // Do convolution calcution (for simulation).

    template <class T> void Convolute(size_t channel, std::vector<T>& func) const;
    template <class T> void Convolute(size_t channel, size_t wire, std::vector<T>& func) const;

    // Do deconvolution calcution (for reconstruction).

    template <class T> void Deconvolute(size_t channel, std::vector<T>& func) const;
    template <class T> void Deconvolute(size_t channel, size_t wire, std::vector<T>& func) const;
    
    void SetDecon(int fftsize);
    double GetDeconNorm(){return fDeconNorm;};

  private:

    // Private configuration methods.

    // Post-constructor initialization.

    void init() const{const_cast<SignalShapingServiceMicroBooNE*>(this)->init();}
    void init();

    // Calculate response functions.
    // Copied from SimWireMicroBooNE.


    void SetFieldResponse(size_t ktype);
    // void SetElectResponse(size_t ktype);
   
    void SetElectResponse(size_t ktype, double shapingtime, double gain);  //changed to read different peaking time for different planes

    // Calculate filter functions.

    void SetFilters();

    // Attributes.

    bool fInit;               ///< Initialization flag.

    // Sample the response function, including a configurable
    // drift velocity of electrons

    void SetResponseSampling(size_t ktype);

    void SetFieldResponseTOffsets( const TH1F* resp, const size_t ktype);

    size_t fNPlanes;
    size_t fNViews;

    
    

    // Fcl parameters.
    std::vector<size_t>      fViewIndex;
    std::map<size_t, size_t> fViewMap;
    size_t                   fViewForNormalization;
    double fDeconNorm;
    double fADCPerPCAtLowestASICGain; ///< Pulse amplitude gain for a 1 pc charge impulse after convoluting it the with field and electronics response with the lowest ASIC gain setting of 4.7 mV/fC

	  //double fASICGainInMVPerFC;                ///< Cold electronics ASIC gain setting in mV/fC
    std::vector<DoubleVec> fNoiseFactVec;       ///< RMS noise in ADCs for lowest gain setting

    std::vector<std::vector<size_t> > fNResponses;
    std::vector<std::vector<size_t> > fNActiveResponses;
    //double fASICGainInMVPerFC;                  ///< Cold electronics ASIC gain setting in mV/fC
    std::vector<double> fASICGainInMVPerFC;       ///< Cold electronics ASIC gain setting in mV/fC
    //std::vector<double> fNoiseFactColl;         ///< RMS Noise in ADCs for lowest Gain Setting
    //std::vector<double> fNoiseFactInd;          ///< RMS Noise in ADCs for lowest Gain Setting

    DoubleVec fDefaultDriftVelocity;  ///< Default drift velocity of electrons in cm/usec
    std::vector<DoubleVec>  fFieldResponseTOffset;  ///< Time offset for field response in ns

    std::vector<double> fCalibResponseTOffset; // calibrated time offset to align U/V/Y Signals 

    int fNFieldBins[2];         		///< number of bins for field response
    int fFieldLowEdge[2];           ///< low edge of the field response histo (for test output)
    double fFieldBinWidth[2];       ///<  Bin with of the input field response.

    DoubleVec f3DCorrectionVec;  ///< correction factor to account for 3D path of electrons, 1 for each plane (default = 1.0)
    
    DoubleVec fFieldRespAmpVec;
    std::vector<double> fShapeTimeConst; ///< time constants for exponential shaping
    std::vector<int> fDeconvPol;         ///< switch for DeconvKernel normalization sign (+ -> max pos ADC, - -> max neg ADC). Entry 0,1,2 = U,V,Y plane settings
    std::vector<TF1*> fFilterTF1Vec;     ///< Vector of Parameterized filter functions
    std::vector<std::string> fFilterFuncVec;
    std::vector<std::vector<TComplex> > fFilterVec;

    // Induced charge deconvolution additions (M. Mooney)
    std::vector<TF1*> fFilterTF1VecICTime;
    std::vector<std::string> fFilterFuncVecICTime;
    std::vector<TF1*> fFilterTF1VecICWire;
    std::vector<std::string> fFilterFuncVecICWire;
    std::vector<double> fFilterScaleVecICTime;
    std::vector<double> fFilterScaleVecICWire;
    std::vector<double> fFilterNormVecIC;

    std::vector<double> fFilterICTimeMaxFreq;
    std::vector<double> fFilterICTimeMaxVal;

    std::vector<double> fFilterICWireMaxFreq;
    std::vector<double> fFilterICWireMaxVal;

    

    bool fGetFilterFromHisto;   		///< Flag that allows to use a filter function from a histogram instead of the functional dependency

    std::vector<std::vector<std::vector<TH1F*> > > fFieldResponseHistVec;

    double fDefaultEField;
    double fDefaultTemperature;

    DoubleVec fTimeScaleParams;
    
    std::vector<TH1D*> fFilterHistVec;
    
    // Following attributes hold the convolution and deconvolution kernels

    std::vector<std::vector<std::vector<util::SignalShaping> > >fSignalShapingVec;
    // Field response.

    std::vector<std::vector<std::vector<DoubleVec> > >fFieldResponseVec;

    // Electronics response.

    std::vector<DoubleVec> fElectResponse;

    // Filters.

    std::vector<std::vector<TComplex> > FilterVec;

    bool fPrintResponses;
    bool fManualInterpolation;
  };
}




//----------------------------------------------------------------------
// Do convolution.
template <class T> inline void util::SignalShapingServiceMicroBooNE::Convolute(size_t channel, std::vector<T>& func) const
{
  SignalShaping(channel, 0).Convolute(func);

  //negative number
  int time_offset = FieldResponseTOffset(channel,0);
  
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

// Do convolution.
template <class T> inline void util::SignalShapingServiceMicroBooNE::Convolute(size_t channel, size_t wire, std::vector<T>& func) const
{
  SignalShaping(channel, wire).Convolute(func);

  //negative number
  int time_offset = FieldResponseTOffset(channel,0);
  
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
template <class T> inline void util::SignalShapingServiceMicroBooNE::Deconvolute(size_t channel, std::vector<T>& func) const
{
  size_t ktype = 1;
  SignalShaping(channel, 0, ktype).Deconvolute(func);

  int time_offset = FieldResponseTOffset(channel,ktype);
  
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

//----------------------------------------------------------------------
// Do deconvolution.

template <class T> inline void util::SignalShapingServiceMicroBooNE::Deconvolute(size_t channel, size_t wire, std::vector<T>& func) const
{
  size_t ktype = 1;
  SignalShaping(channel, wire, ktype).Deconvolute(func);

  int time_offset = FieldResponseTOffset(channel,ktype);
  
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

