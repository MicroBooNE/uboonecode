
#include "RawDigitFFTAlg.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geometry/Geometry.h"
#include "Utilities/DetectorProperties.h"

#include <cmath>
#include <algorithm>

#include "TVirtualFFT.h"

namespace caldata
{

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
RawDigitFFTAlg::RawDigitFFTAlg(fhicl::ParameterSet const & pset)
{
    reconfigure(pset);

    // Report.
    mf::LogInfo("RawDigitFFTAlg") << "RawDigitFFTAlg configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
RawDigitFFTAlg::~RawDigitFFTAlg()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void RawDigitFFTAlg::reconfigure(fhicl::ParameterSet const & pset)
{
}

template <class T> void RawDigitFFTAlg::getFFTCorrection(std::vector<T>& corValVec, double minPowerThreshold) const
{
    // This version will take FFT of input waveform and then remove bins in the time domain with a power less
    // than the threshold input above.
    int fftDataSize = corValVec.size();
    
    TVirtualFFT* fftr2c = TVirtualFFT::FFT(1, &fftDataSize, "R2C");

    // In the first step we copy the input array into a container of doubles to pass to the FFT
    double fftInputArray[fftDataSize];
    
    std::copy(corValVec.begin(),corValVec.end(),fftInputArray);
    
    fftr2c->SetPoints(fftInputArray);
    fftr2c->Transform();
    
    // In the second step we recover the power spectrum
    double realVals[fftDataSize];
    double imaginaryVals[fftDataSize];
    
    fftr2c->GetPointsComplex(realVals, imaginaryVals);
    
    size_t halfFFTDataSize(fftDataSize/2);
    
    double powerVec[halfFFTDataSize];
    
    std::transform(realVals, realVals + halfFFTDataSize, imaginaryVals, powerVec, [](const double& real, const double& imaginary){return std::sqrt(real*real + imaginary*imaginary);});

    // Third step is to zap those bins under threshold
    for(size_t idx = 0; idx < halfFFTDataSize; idx++)
    {
        if (powerVec[idx] < minPowerThreshold)
        {
            realVals[idx]                  = 0.;
            realVals[fftDataSize-idx]      = 0.;
            imaginaryVals[idx]             = 0.;
            imaginaryVals[fftDataSize-idx] = 0.;
        }
    }

    // Finally, we invert the resulting time domain values to recover the new waveform
    TVirtualFFT* fftc2r = TVirtualFFT::FFT(1, &fftDataSize, "C2R M K");
    
    fftc2r->SetPointsComplex(realVals,imaginaryVals);
    fftc2r->Transform();
    
    double* fftOutputArray = fftc2r->GetPointsReal();
    
    double normFctr = 1. / double(fftDataSize);
    
    std::transform(fftOutputArray, fftOutputArray + fftDataSize, corValVec.begin(), [normFctr](const double& real){return real * normFctr;});
    
    delete fftc2r;
    
    return;
}
    
template void RawDigitFFTAlg::getFFTCorrection<float>(std::vector<float>&, double) const;

template<class T> void RawDigitFFTAlg::getFFTCorrection(std::vector<T>& corValVec, size_t maxBin) const
{
    // This version will take FFT of input waveform and then remove bins in the time domain above the
    // cutoff frequency defined by maxBin passed in above
    int fftDataSize = corValVec.size();
    
    TVirtualFFT* fftr2c = TVirtualFFT::FFT(1, &fftDataSize, "R2C");
    
    // In the first step we copy the input array into a container of doubles to pass to the FFT
    double fftInputArray[fftDataSize];
    
    std::copy(corValVec.begin(),corValVec.end(),fftInputArray);
    
    fftr2c->SetPoints(fftInputArray);
    fftr2c->Transform();
    
    // In the second step we recover the power spectrum
    double realVals[fftDataSize];
    double imaginaryVals[fftDataSize];
    
    fftr2c->GetPointsComplex(realVals, imaginaryVals);
    
    size_t halfFFTDataSize(fftDataSize/2);
    
    double powerVec[halfFFTDataSize];
    
    std::transform(realVals, realVals + halfFFTDataSize, imaginaryVals, powerVec, [](const double& real, const double& imaginary){return std::sqrt(real*real + imaginary*imaginary);});
    
    // Zero all bins above selected frequency
    std::fill(realVals      + maxBin, realVals      + fftDataSize - maxBin, 0.);
    std::fill(imaginaryVals + maxBin, imaginaryVals + fftDataSize - maxBin, 0.);
    
    // Finally, we invert the resulting time domain values to recover the new waveform
    TVirtualFFT* fftc2r = TVirtualFFT::FFT(1, &fftDataSize, "C2R M K");
    
    fftc2r->SetPointsComplex(realVals,imaginaryVals);
    fftc2r->Transform();
    
    double* fftOutputArray = fftc2r->GetPointsReal();
    
    double normFctr = 1. / double(fftDataSize);
    
    std::transform(fftOutputArray, fftOutputArray + fftDataSize, corValVec.begin(), [normFctr](const double& real){return real * normFctr;});
    
    delete fftc2r;
    
    return;
}

template void RawDigitFFTAlg::getFFTCorrection<float>(std::vector<float>& corValVec, size_t maxBin) const;

    
}
