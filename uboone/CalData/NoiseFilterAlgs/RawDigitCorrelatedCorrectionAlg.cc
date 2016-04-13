
#include <cmath>
#include <algorithm>
#include <vector>

#include "RawDigitCorrelatedCorrectionAlg.h"

#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

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
RawDigitCorrelatedCorrectionAlg::RawDigitCorrelatedCorrectionAlg(fhicl::ParameterSet const & pset) :
    fFFTAlg(pset)
{
    reconfigure(pset);

    // Report.
    mf::LogInfo("RawDigitCorrelatedCorrectionAlg") << "RawDigitCorrelatedCorrectionAlg configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
RawDigitCorrelatedCorrectionAlg::~RawDigitCorrelatedCorrectionAlg()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void RawDigitCorrelatedCorrectionAlg::reconfigure(fhicl::ParameterSet const & pset)
{
    fTruncMeanFraction     = pset.get<float>              ("TruncMeanFraction",                                        0.15);
    fApplyCorSmoothing     = pset.get<bool>               ("ApplyCorSmoothing",                                        true);
    fApplyFFTCorrection    = pset.get<bool>               ("ApplyFFTCorrection",                                       true);
    fFillFFTHistograms     = pset.get<bool>               ("FillFFTHistograms",                                       false);
    fFFTHistsWireGroup     = pset.get<std::vector<size_t>>("FFTHistsWireGroup",         std::vector<size_t>() = {1, 33, 34});
    fFFTNumHists           = pset.get<std::vector<size_t>>("FFTNumWaveHistograms",       std::vector<size_t>() = {10,48,48});
    fFFTHistsStartTick     = pset.get<std::vector<double>>("FFTWaveHistsStartTick", std::vector<double>() = {96.,96.,7670.});
    fFFTMinPowerThreshold  = pset.get<std::vector<double>>("FFTPowerThreshold",     std::vector<double>() = {100.,75.,500.});
    fNumWiresToGroup       = pset.get<std::vector<size_t>>("NumWiresToGroup",          std::vector<size_t>() = {48, 48, 96});
    fFillHistograms        = pset.get<bool>               ("FillHistograms",                                          false);
    fRunFFTCorrected       = pset.get<bool>               ("RunFFTCorrectedWires",                                    false);
    fNumRmsToSmoothVec     = pset.get<std::vector<float>> ("NumRmsToSmooth",          std::vector<float>() = {3.6, 3.6, 4.});
}

//----------------------------------------------------------------------------
/// Begin job method.
void RawDigitCorrelatedCorrectionAlg::initializeHists(art::ServiceHandle<art::TFileService>& tfs)
{
    // Following to determine min/max frequencies
    double sampleRate  = fDetectorProperties->SamplingRate();
    double readOutSize = fDetectorProperties->ReadOutWindowSize();
    double maxFreq     = 1000000. / (2. * sampleRate);
    double minFreq     = 1000000. / (2. * sampleRate * readOutSize);
    int    numSamples  = (readOutSize / 2 + 1) / 4;
    
    // quick aside
//    std::cout << "++> plane 0 offset: " << fDetectorProperties->GetXTicksOffset(0,0,0) << std::endl;
//    std::cout << "++> plane 1 offset: " << fDetectorProperties->GetXTicksOffset(1,0,0) << std::endl;
//    std::cout << "++> plane 2 offset: " << fDetectorProperties->GetXTicksOffset(2,0,0) << std::endl;
    
    fFFTHist[0]       = tfs->make<TProfile>("FFTPlaneU",    "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    fFFTHist[1]       = tfs->make<TProfile>("FFTPlaneV",    "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    fFFTHist[2]       = tfs->make<TProfile>("FFTPlaneW",    "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    
    fFFTHistLow[0]    = tfs->make<TProfile>("FFTLowPlaneU", "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    fFFTHistLow[1]    = tfs->make<TProfile>("FFTLowPlaneV", "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    fFFTHistLow[2]    = tfs->make<TProfile>("FFTLowPlaneW", "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    
    fFFTHistCor[0]    = tfs->make<TProfile>("FFTCorPlaneU", "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    fFFTHistCor[1]    = tfs->make<TProfile>("FFTCorPlaneV", "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    fFFTHistCor[2]    = tfs->make<TProfile>("FFTCorPlaneW", "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);

    fCorValHist[0]    = tfs->make<TProfile>("FFTCorValU",   "Raw Corrections;Tick",     readOutSize,     0., readOutSize, -100., 100.);
    fCorValHist[1]    = tfs->make<TProfile>("FFTCorValV",   "Raw Corrections;Tick",     readOutSize,     0., readOutSize, -100., 100.);
    fCorValHist[2]    = tfs->make<TProfile>("FFTCorValW",   "Raw Corrections;Tick",     readOutSize,     0., readOutSize, -100., 100.);
    
    fFFTCorValHist[0] = tfs->make<TProfile>("FFTHistValU",  "Power Spectrum;kHz;Power", numSamples, minFreq,     maxFreq,    0., 10000.);
    fFFTCorValHist[1] = tfs->make<TProfile>("FFTHistValV",  "Power Spectrum;kHz;Power", numSamples, minFreq,     maxFreq,    0., 10000.);
    fFFTCorValHist[2] = tfs->make<TProfile>("FFTHistValW",  "Power Spectrum;kHz;Power", numSamples, minFreq,     maxFreq,    0., 10000.);
    
    fFFTCorHist[0]    = tfs->make<TProfile>("FFTCorHistU",  "FFT Corrections;Tick",     readOutSize,     0., readOutSize, -100., 100.);
    fFFTCorHist[1]    = tfs->make<TProfile>("FFTCorHistV",  "FFT Corrections;Tick",     readOutSize,     0., readOutSize, -100., 100.);
    fFFTCorHist[2]    = tfs->make<TProfile>("FFTCorHistW",  "FFT Corrections;Tick",     readOutSize,     0., readOutSize, -100., 100.);
    
    if (fFillFFTHistograms)
    {
        fInputWaveHists.resize(3);
        fOutputWaveHists.resize(3);
        fStdCorWaveHists.resize(3);
        fWaveformProfHists.resize(3);
        fWaveCorHists.resize(3);
    
        for(size_t viewIdx = 0; viewIdx < 3; viewIdx++)
        {
            fInputWaveHists[viewIdx].resize(fNumWiresToGroup[viewIdx]);
            fOutputWaveHists[viewIdx].resize(fNumWiresToGroup[viewIdx]);
            fStdCorWaveHists[viewIdx].resize(fNumWiresToGroup[viewIdx]);
            fWaveformProfHists[viewIdx].resize(fNumWiresToGroup[viewIdx]);

            for(size_t histIdx = 0; histIdx < fNumWiresToGroup[viewIdx]; histIdx++)
            {
                std::string histName = "InputWaveform_" + std::to_string(viewIdx) + "_" + std::to_string(histIdx);
            
                fInputWaveHists[viewIdx][histIdx] = tfs->make<TProfile>(histName.c_str(), "Raw Waveform;Tick;dADC", readOutSize, 0., readOutSize, -500., 500.);
            
                histName = "OutputWaveform_" + std::to_string(viewIdx) + "_" + std::to_string(histIdx);
            
                fOutputWaveHists[viewIdx][histIdx] = tfs->make<TProfile>(histName.c_str(), "FFT Corrected Waveform;Tick;dADC", readOutSize, 0., readOutSize, -500., 500.);
            
                histName = "StdCorWaveform_" + std::to_string(viewIdx) + "_" + std::to_string(histIdx);
            
                fStdCorWaveHists[viewIdx][histIdx] = tfs->make<TProfile>(histName.c_str(), "Standard Corrected Waveform;Tick;dADC", readOutSize, 0., readOutSize, -500., 500.);
                
                histName = "InputWaveform_prof_" + std::to_string(viewIdx) + "_" + std::to_string(histIdx);
                
                fWaveformProfHists[viewIdx][histIdx] = tfs->make<TH1D>(histName.c_str(), "Waveform Profile;dADC;Count", 400, -200., 200.);
            }
            
            fWaveCorHists[viewIdx].resize(fFFTNumHists[viewIdx]);
            
            for(size_t histIdx = 0; histIdx < fFFTNumHists[viewIdx]; histIdx++)
            {
                std::string histName = "WaveCor_" + std::to_string(viewIdx) + "_tick_" + std::to_string(fFFTHistsStartTick[viewIdx]+histIdx);
                
                fWaveCorHists[viewIdx][histIdx] = tfs->make<TProfile>(histName.c_str(), "Tick;ADC", fNumWiresToGroup[viewIdx], 0., fNumWiresToGroup[viewIdx], -500., 500.);
            }
        }
    }
    
    fFFTvsMBProf[0]   = tfs->make<TProfile2D>("FFTvsMBPlaneU", "FFT;MB;kHz", numSamples, minFreq, maxFreq, 2400/16, 0., 2400/16);
    fFFTvsMBProf[1]   = tfs->make<TProfile2D>("FFTvsMBPlaneV", "FFT;MB;kHz", numSamples, minFreq, maxFreq, 2400/16, 0., 2400/16);
    fFFTvsMBProf[2]   = tfs->make<TProfile2D>("FFTvsMBPlaneW", "FFT;MB;kHz", numSamples, minFreq, maxFreq, 3456/16, 0., 3456/16);
}

void RawDigitCorrelatedCorrectionAlg::smoothCorrectionVec(std::vector<float>& corValVec, unsigned int& viewIdx) const
{
    // First get the truncated mean and rms for the input vector (noting that it is not in same format as raw data)
    // We need a local copy so we can sort it
    std::vector<float> localCorValVec = corValVec;
    
    std::sort(localCorValVec.begin(),localCorValVec.end());

    int   nTruncVal  = (1. - fTruncMeanFraction) * localCorValVec.size();
    float corValSum  = std::accumulate(localCorValVec.begin(),localCorValVec.begin() + nTruncVal,0.);
    float meanCorVal = corValSum / float(nTruncVal);
    
    std::vector<float> diffVec(nTruncVal);
    std::transform(localCorValVec.begin(),localCorValVec.begin() + nTruncVal, diffVec.begin(), std::bind2nd(std::minus<float>(),meanCorVal));
    
    float rmsValSq   = std::inner_product(diffVec.begin(),diffVec.end(),diffVec.begin(),0.);
    float rmsVal     = std::sqrt(rmsValSq / float(nTruncVal));
    
    // Now set up to run through and do a "simple" interpolation over outliers
    std::vector<float>::iterator lastGoodItr = corValVec.begin();
    
    bool wasOutlier(false);
    
    for(std::vector<float>::iterator corValItr = lastGoodItr+1; corValItr != corValVec.end(); corValItr++)
    {
        if (fabs(*corValItr - meanCorVal) < fNumRmsToSmoothVec.at(viewIdx)*rmsVal)
        {
            if (wasOutlier)
            {
                float lastVal  = *lastGoodItr;
                float curVal   = *corValItr;
                float numTicks = std::distance(lastGoodItr,corValItr);
                float slope    = (curVal - lastVal) / numTicks;
                
                while(lastGoodItr != corValItr)
                {
                    *lastGoodItr++ = (numTicks - std::distance(lastGoodItr,corValItr)) * slope + lastVal;
                }
            }
            
            wasOutlier  = false;
            lastGoodItr = corValItr;
        }
        else wasOutlier = true;
    }
    
    return;
}

void RawDigitCorrelatedCorrectionAlg::removeCorrelatedNoise(RawDigitAdcIdxPair& digitIdxPair,
                                                            unsigned int        viewIdx,
                                                            std::vector<float>& truncMeanWireVec,
                                                            std::vector<float>& truncRmsWireVec,
                                                            std::vector<short>& minMaxWireVec,
                                                            std::vector<short>& meanWireVec,
                                                            std::vector<float>& skewnessWireVec,
                                                            std::vector<float>& neighborRatioWireVec,
                                                            std::vector<float>& pedCorWireVec) const
{
    // This method represents and enhanced implementation of "Corey's Algorithm" for correcting the
    // correlated noise across a group of wires. The primary enhancement involves using a FFT to
    // "fit" for the underlying noise as a way to reduce the impact on the signal.
    WireToRawDigitVecMap& wireToRawDigitVecMap = digitIdxPair.first;
    WireToAdcIdxMap&      wireToAdcIdxMap      = digitIdxPair.second;
    
    size_t maxTimeSamples(wireToRawDigitVecMap.begin()->second.size());
    size_t baseWireIdx(wireToRawDigitVecMap.begin()->first - wireToRawDigitVecMap.begin()->first % fNumWiresToGroup[viewIdx]);
    
    // Should we turn on our diagnostics blocks?
    bool doFFTCorrection(fFillFFTHistograms && baseWireIdx / fNumWiresToGroup[viewIdx] == fFFTHistsWireGroup[viewIdx]);
    
    // Keep track of the "waveform" for the correction
    std::vector<float> origCorValVec;

    // Don't try to do correction if too few wires unless they have gaps
    if (wireToAdcIdxMap.size() > 2) // || largestGapSize > 2)
    {
        // What is the average min/max of this group?
//        short aveMinMax = std::round(float(std::accumulate(minMaxByWireVec.begin(),minMaxByWireVec.end(),0)) / float(minMaxByWireVec.size()));
        
        std::vector<float> corValVec;
        
        corValVec.resize(maxTimeSamples, 0.);
        
        // Build the vector of corrections for each time bin
        for(size_t sampleIdx = 0; sampleIdx < maxTimeSamples; sampleIdx++)
        {
            // Try to do the calculations for histogramming outside of the loop over channels
            bool   fillHists(false);
            size_t histIdx(0);

            // Diagnostics
            if (doFFTCorrection)
            {
                fillHists = sampleIdx >= fFFTHistsStartTick[viewIdx] && sampleIdx < fFFTHistsStartTick[viewIdx] + fFFTNumHists[viewIdx];
                histIdx   = sampleIdx - fFFTHistsStartTick[viewIdx];
            }
            
            // Define a vector for accumulating values...
            // Loop over the wires at this time bin and get their pedestal corrected ADC values
            // We'll use a simple stl vector for this
            std::vector<float> adcValuesVec;
            
            for(const auto& wireAdcItr : wireToAdcIdxMap)
            {
                // Check that we should be doing something in this range
                // Note that if the wire is not to be considered then the "start" bin will be after the last bin
                if (sampleIdx < wireAdcItr.second.first || sampleIdx >= wireAdcItr.second.second) continue;
                
                int wireIdx(wireAdcItr.first - baseWireIdx);

                // Accumulate
                adcValuesVec.push_back(float(wireToRawDigitVecMap.at(wireAdcItr.first)[sampleIdx]) - truncMeanWireVec[wireIdx]);
                
                // Make hists if requested
                if (fillHists)
                {
                    fWaveCorHists[viewIdx][histIdx]->Fill(wireIdx, adcValuesVec.back());
                }
            }
            
            float medianValue = getMedian(adcValuesVec, float(-10000.));
            
            corValVec[sampleIdx] = medianValue;
        }
        
        // Try to eliminate any real outliers
        if (fApplyCorSmoothing) smoothCorrectionVec(corValVec, viewIdx);
        
        // Diagnostics block
        if (doFFTCorrection)
        {
            origCorValVec = corValVec;
            
            for(size_t tick = 0; tick < maxTimeSamples; tick++)
            {
                double corVal = origCorValVec[tick];
                
                for(const auto& wireAdcItr : wireToRawDigitVecMap)
                {
                    short  inputWaveformVal = wireAdcItr.second[tick];
                    size_t wireIdx          = wireAdcItr.first % fNumWiresToGroup[viewIdx];
                    double inputVal         = inputWaveformVal - truncMeanWireVec[wireIdx];
                    double outputVal        = std::round(inputVal - corVal);
                    
                    fInputWaveHists[viewIdx][wireIdx]->Fill(tick, inputVal, 1.);
                    fStdCorWaveHists[viewIdx][wireIdx]->Fill(tick, outputVal, 1.);
                    fWaveformProfHists[viewIdx][wireIdx]->Fill(std::round(inputVal) + 0.5, 1.);
                }
            }
        }
        
        // Get the FFT correction
        if (fApplyFFTCorrection) fFFTAlg.getFFTCorrection(corValVec,fFFTMinPowerThreshold[viewIdx]);
        
        // Now go through and apply the correction
        for(size_t sampleIdx = 0; sampleIdx < maxTimeSamples; sampleIdx++)
        {
            // Now run through and apply correction
            for (const auto& wireAdcItr : wireToAdcIdxMap)
            {
                float corVal(corValVec[sampleIdx]);
                int   wireIdx(wireAdcItr.first);
                
                // If the "start" bin is after the "stop" bin then we are meant to skip this wire in the averaging process
                // Or if the sample index is in a chirping section then no correction is applied.
                // Both cases are handled by looking at the sampleIdx
                if (sampleIdx < wireAdcItr.second.first || sampleIdx >= wireAdcItr.second.second)
                    corVal = 0.;
                
                //RawDigitVector& rawDataTimeVec = wireToRawDigitVecMap.at(wireIdx);
                short& rawDataTimeVal = wireToRawDigitVecMap.at(wireIdx).at(sampleIdx);
                
                // Probably doesn't matter, but try to get slightly more accuracy by doing float math and rounding
                float newAdcValueFloat = float(rawDataTimeVal) - corVal - pedCorWireVec[wireIdx - baseWireIdx];
                
                rawDataTimeVal = std::round(newAdcValueFloat);
                
                if (doFFTCorrection)
                {
                    float stdCorWaveVal = std::round(newAdcValueFloat + pedCorWireVec[wireIdx - baseWireIdx] - truncMeanWireVec[wireIdx - baseWireIdx]);
                    
                    fOutputWaveHists[viewIdx][wireIdx - baseWireIdx]->Fill(sampleIdx, stdCorWaveVal, 1.);
                    
                    fFFTCorHist[viewIdx]->Fill(sampleIdx,corVal);
                }
            }
        }
    }
    
    // Final diagnostics block
    if (doFFTCorrection && !origCorValVec.empty())
    {
        double sampleFreq  = 1000000. / fDetectorProperties->SamplingRate();
        double readOutSize = fDetectorProperties->ReadOutWindowSize();
        int    fftDataSize = origCorValVec.size();
        
        TVirtualFFT* fftr2c = TVirtualFFT::FFT(1, &fftDataSize, "R2C");
        
        double fftInputArray[fftDataSize];
        
        for(size_t tick = 0; tick < size_t(fftDataSize); tick++)
        {
            fftInputArray[tick] = origCorValVec[tick];
            
            fCorValHist[viewIdx]->Fill(tick,origCorValVec[tick]);
        }
        
        fftr2c->SetPoints(fftInputArray);
        fftr2c->Transform();
        
        // Recover the power spectrum...
        double realVals[fftDataSize];
        double imaginaryVals[fftDataSize];
        double realPart(0.);
        double imaginaryPart(0.);
        
        fftr2c->GetPointsComplex(realVals, imaginaryVals);
            
        for(size_t idx = 0; idx < size_t(fftDataSize)/2; idx++)
        {
            realPart      = realVals[idx];
            imaginaryPart = imaginaryVals[idx];
                
            double bin   = (idx * sampleFreq) / readOutSize;
            double power = std::sqrt(realPart*realPart + imaginaryPart*imaginaryPart);
            
            fFFTCorValHist[viewIdx]->Fill(bin, power, 1.);
        }
    }
    
    // Run an FFT here to check our "corrected" wires
    if (fRunFFTCorrected)
    {
        double sampleFreq  = 1000000. / fDetectorProperties->SamplingRate();
        double readOutSize = fDetectorProperties->ReadOutWindowSize();
        
        for(const auto& wireAdcItr : wireToAdcIdxMap)
        {
            if (wireAdcItr.second.first > wireAdcItr.second.second) continue;

            RawDigitVector& rawDataTimeVec = wireToRawDigitVecMap.at(wireAdcItr.first);
            
            size_t wireIdx = wireAdcItr.first  - baseWireIdx;
            
            int fftDataSize = rawDataTimeVec.size();
        
            TVirtualFFT* fftr2c = TVirtualFFT::FFT(1, &fftDataSize, "R2C");
        
            double fftInputArray[fftDataSize];
        
            for(size_t idx = 0; idx < size_t(fftDataSize); idx++) fftInputArray[idx] = rawDataTimeVec[idx] - truncMeanWireVec[wireIdx];
        
            fftr2c->SetPoints(fftInputArray);
            fftr2c->Transform();
        
            // Recover the power spectrum...
            double realPart(0.);
            double imaginaryPart(0.);
        
            for(size_t idx = 0; idx < size_t(fftDataSize)/2; idx++)
            {
                fftr2c->GetPointComplex(idx+1, realPart, imaginaryPart);
            
                double bin   = (idx * sampleFreq) / readOutSize;
                double power = realPart*realPart + imaginaryPart*imaginaryPart;
            
                if (power > 0.) power = std::sqrt(power);
            
                fFFTHistCor[viewIdx]->Fill(bin, power, 1.);
            }
        }
    }
    
    return;
}
    
template<class T> T RawDigitCorrelatedCorrectionAlg::getMedian(std::vector<T>& valuesVec, T defaultValue) const
{
    T medianValue(defaultValue);
    
    if (!valuesVec.empty())
    {
        std::sort(valuesVec.begin(),valuesVec.end());
        
        size_t medianIdx = valuesVec.size() / 2;
        
        medianValue = valuesVec[medianIdx];
        
        if (valuesVec.size() > 1 && medianIdx % 2) medianValue = (medianValue + valuesVec[medianIdx+1]) / 2;
    }
    
    return std::max(medianValue,defaultValue);
}
    
}
