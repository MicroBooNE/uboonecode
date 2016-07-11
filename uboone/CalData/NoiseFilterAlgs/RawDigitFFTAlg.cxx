
#include "RawDigitFFTAlg.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

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
    fTransformViewVec   = pset.get<std::vector<bool>>("TransformViewVec",     std::vector<bool>() = {true,false,false});
    fZigZagCorrectVec   = pset.get<std::vector<bool>>("ZigZagCorrectVec",     std::vector<bool>() = {true,true,false} );
    fZigZagCorrectBin   = pset.get<size_t           >("ZigZagCorrectBin",                                         2800);
    fFillHistograms     = pset.get<bool             >("FillHistograms",                                          false);
    fHistDirName        = pset.get<std::string      >("HistDirName",                                       "FFT_hists");
    
}
    
//----------------------------------------------------------------------------
/// Begin job method.
void RawDigitFFTAlg::initializeHists(art::ServiceHandle<art::TFileService>& tfs)
{
    if (fFillHistograms)
    {
        // Define the histograms. Putting semi-colons around the title
        // causes it to be displayed as the x-axis label if the histogram
        // is drawn.
        
        // hijack hists here
        //    double sampleRate  = fDetectorProperties->SamplingRate();
        double readOutSize = fDetectorProperties->ReadOutWindowSize();
        //    double maxFreq     = 1000000. / (2. * sampleRate);
        //    double minFreq     = 1000000. / (2. * sampleRate * readOutSize);
        //    int    numSamples  = (readOutSize / 2 + 1) / 4;
        int numSamples     = readOutSize / 2;
        
        fCorValHistVec.resize(20);
        fFFTPowerVec.resize(20);
        fFFTCorValHistVec.resize(20);
        
        // Make a directory for these histograms
        art::TFileDirectory dir = tfs->mkdir(fHistDirName.c_str());
        
        
        for(size_t idx = 0; idx < 20; idx++)
        {
            std::string histName = "RawWaveform_" + std::to_string(idx);
            
            fCorValHistVec[idx] = dir.make<TProfile>(histName.c_str(), "Raw Waveform;Tick", readOutSize, 0., readOutSize, -100., 100.);
            
            histName = "FFTPower_" + std::to_string(idx);
            
            //fFFTPowerVec[idx] = tfs->make<TProfile>(histName.c_str(),  "Power Spectrum;kHz;Power", numSamples, minFreq, maxFreq, 0., 10000.);
            fFFTPowerVec[idx] = dir.make<TProfile>(histName.c_str(),  "Power Spectrum;kHz;Power", numSamples, 0, numSamples, 0., 10000.);
            
            histName = "FFTCorrected_" + std::to_string(idx);
            
            fFFTCorValHistVec[idx] = dir.make<TProfile>(histName.c_str(),  "Corrected Waveform;Tick", readOutSize, 0., readOutSize, -100., 100.);
        }
    }
    
    return;
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
    delete fftr2c;
    
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
    delete fftr2c;
    
    return;
}
    
void RawDigitFFTAlg::filterFFT(std::vector<short>& rawadc, size_t view, size_t wire, float pedestal) const
{
    if (fTransformViewVec[view] || fZigZagCorrectVec[view])
    {
        size_t lowWire(1500);
        size_t hiWire(1520);
    
        //                double sampleFreq  = 1000000. / fDetectorProperties->SamplingRate();
        //                double readOutSize = fDetectorProperties->ReadOutWindowSize();
        //                double binSize     = sampleFreq / readOutSize;
        int    fftDataSize = rawadc.size();
    
        TVirtualFFT* fftr2c = TVirtualFFT::FFT(1, &fftDataSize, "R2C");
    
        double fftInputArray[fftDataSize];
    
        for(size_t tick = 0; tick < size_t(fftDataSize); tick++)
        {
            fftInputArray[tick] = rawadc[tick] - pedestal;
        
            if (fFillHistograms && view == 0 && wire >= lowWire && wire < hiWire)
                fCorValHistVec[wire-lowWire]->Fill(tick, fftInputArray[tick], 1.);
        }
    
        fftr2c->SetPoints(fftInputArray);
        fftr2c->Transform();
    
        // Recover the power spectrum...
        double realVals[fftDataSize];
        double imaginaryVals[fftDataSize];
    
        fftr2c->GetPointsComplex(realVals, imaginaryVals);
    
        size_t halfFFTDataSize(fftDataSize/2);
    
        std::vector<double> powerVec(halfFFTDataSize);
    
        std::transform(realVals, realVals + halfFFTDataSize, imaginaryVals, powerVec.begin(), [](const double& real, const double& imaginary){return std::sqrt(real*real + imaginary*imaginary);});
    
        if (fFillHistograms && view == 0 && wire >= lowWire && wire < hiWire)
        {
            // Third step is to zap those bins under threshold
            for(size_t idx = 0; idx < halfFFTDataSize; idx++)
            {
                //                    double bin = (double(idx) + 0.5) * binSize;
            
                //fFFTPowerVec[channel-1500]->Fill(bin, powerVec[idx], 1.);
                fFFTPowerVec[wire-lowWire]->Fill(idx, std::min(powerVec[idx],9999.), 1.);
            }
        }
    
        if (fTransformViewVec[view])
        {
/*
            std::vector<double>::iterator powerItr = std::max_element(powerVec.begin(),powerVec.end());

            while(*powerItr > 4000.)
            {
                size_t idx = std::distance(powerVec.begin(),powerItr);
    
                if (idx > 0)
                {
                    realVals[idx]                  = 0.;
                    realVals[fftDataSize-idx]      = 0.;
                    imaginaryVals[idx]             = 0.;
                    imaginaryVals[fftDataSize-idx] = 0.;
                }
    
                *powerItr = 0.;
    
                powerItr = std::max_element(powerVec.begin(),powerVec.end());
            }
*/
            // These bins correspond to the 36 kHz line
            size_t magic_bins[] = {113, 114, 115, 116}; //, 342, 343, 344, 345, 346};
            
            for(const auto& idx : magic_bins)
            {
                realVals[idx]                  = 0.;
                realVals[fftDataSize-idx]      = 0.;
                imaginaryVals[idx]             = 0.;
                imaginaryVals[fftDataSize-idx] = 0.;
            }
        }

        if (fZigZagCorrectVec[view])
        {
            for(size_t idx = fZigZagCorrectBin; idx < halfFFTDataSize; idx++)
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
    
        std::transform(fftOutputArray, fftOutputArray + fftDataSize, rawadc.begin(), [normFctr,pedestal](const double& real){return std::round(real * normFctr + pedestal);});
        
        if (fFillHistograms && view == 0 && wire >= lowWire && wire < hiWire)
            for(int idx = 0; idx < fftDataSize; idx++) fFFTCorValHistVec[wire-lowWire]->Fill(idx, rawadc[idx] - pedestal, 1.);
        
        delete fftc2r;
        delete fftr2c;
    }
    
    return;
}

template void RawDigitFFTAlg::getFFTCorrection<float>(std::vector<float>& corValVec, size_t maxBin) const;

    
}
