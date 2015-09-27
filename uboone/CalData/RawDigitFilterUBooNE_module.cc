////////////////////////////////////////////////////////////////////////
//
// Class:       RawDigitFilterUBooNE
// Module Type: producer
// File:        RawDigitFilterUBooNE_module.cc
//
//              The intent of this module is to filter out "bad" channels
//              in an input RawDigit data stream. In the current implementation,
//              "bad" is defined as the truncated rms for a channel crossing
//              a user controlled threshold
//
// Configuration parameters:
//
// DigitModuleLabel      - the source of the RawDigit collection
// TruncMeanFraction     - the fraction of waveform bins to discard when
//                         computing the means and rms
// RMSRejectionCutHi     - vector of maximum allowed rms values to keep channel
// RMSRejectionCutLow    - vector of lowest allowed rms values to keep channel
// RMSSelectionCut       - vector of rms values below which to not correct
// TheChoseWire          - Wire chosen for "example" hists
// MaxPedestalDiff       - Baseline difference to pedestal to flag
// SmoothCorrelatedNoise - Turns on the correlated noise suppression
// NumWiresToGroup       - When removing correlated noise, # wires to group
// RunFFTInputWires      - FFT analyze the input RawDigits if true
// RunFFTCorrectedWires  - FFT analyze the output RawDigits if true
//
//
// Created by Tracy Usher (usher@slac.stanford.edu) on August 17, 2015
//
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <algorithm>
#include <vector>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geometry/Geometry.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/TimeService.h"
#include "Utilities/SimpleTimeService.h"
#include "CalibrationDBI/Interface/IDetPedestalService.h"
#include "CalibrationDBI/Interface/IDetPedestalProvider.h"

#include "RawData/RawDigit.h"
#include "RawData/raw.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include "TVirtualFFT.h"

class Propagator;

class RawDigitFilterUBooNE : public art::EDProducer
{
public:

    // Copnstructors, destructor.
    explicit RawDigitFilterUBooNE(fhicl::ParameterSet const & pset);
    virtual ~RawDigitFilterUBooNE();

    // Overrides.
    virtual void reconfigure(fhicl::ParameterSet const & pset);
    virtual void produce(art::Event & e);
    virtual void beginJob();
    virtual void endJob();

private:

    // Fcl parameters.
    std::string           fDigitModuleLabel;      ///< The full collection of hits
    float                 fTruncMeanFraction;     ///< Fraction for truncated mean
    std::vector<double>   fRmsRejectionCutHi;     ///< Maximum rms for input channels, reject if larger
    std::vector<double>   fRmsRejectionCutLow;    ///< Minimum rms to consider channel "alive"
    std::vector<double>   fRmsSelectionCut;       ///< Don't use/apply correction to wires below this
    unsigned int          fTheChosenWire;         ///< For example hist
    double                fMaxPedestalDiff;       ///< Max pedestal diff to db to warn
    bool                  fSmoothCorrelatedNoise; ///< Should we smooth the noise?
    std::vector<size_t>   fNumWiresToGroup;       ///< If smoothing, the number of wires to look at
    bool                  fRunFFTInput;           ///< Should we run FFT's on input wires?
    bool                  fRunFFTCorrected;       ///< Should we run FFT's on corrected wires?

    // Statistics.
    int fNumEvent;        ///< Number of events seen.
    
    // Pointers to the histograms we'll create for monitoring what is happening
    TH1D*                 fAdcCntHist[3];
    TH1D*                 fAveValHist[3];
    TH1D*                 fRmsValHist[3];
    TH1D*                 fPedValHist[3];
    TH1D*                 fAverageHist[3];
    TProfile*             fRmsValProf[3];
    TProfile*             fPedValProf[3];

    TProfile*             fFFTHist[3];
    TProfile*             fFFTHistLow[3];
    TProfile*             fFFTHistCor[3];
    
    TProfile2D*           fFFTvsMBProf[3];
    
    bool                  fFirstEvent;
    
    // Useful services, keep copies for now (we can update during begin run periods)
    art::ServiceHandle<geo::Geometry>            fGeometry;             ///< pointer to Geometry service
    art::ServiceHandle<util::DetectorProperties> fDetectorProperties;   ///< Detector properties service
    const lariov::IDetPedestalProvider&          fPedestalRetrievalAlg; ///< Keep track of an instance to the pedestal retrieval alg
};

DEFINE_ART_MODULE(RawDigitFilterUBooNE)

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
RawDigitFilterUBooNE::RawDigitFilterUBooNE(fhicl::ParameterSet const & pset) :
                      fNumEvent(0),
                      fFirstEvent(true),
                      fPedestalRetrievalAlg(art::ServiceHandle<lariov::IDetPedestalService>()->GetPedestalProvider())

{
    reconfigure(pset);
    produces<std::vector<raw::RawDigit> >();

    // Report.
    mf::LogInfo("RawDigitFilterUBooNE") << "RawDigitFilterUBooNE configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
RawDigitFilterUBooNE::~RawDigitFilterUBooNE()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void RawDigitFilterUBooNE::reconfigure(fhicl::ParameterSet const & pset)
{
    fDigitModuleLabel      = pset.get<std::string>        ("DigitModuleLabel",                                       "daq");
    fTruncMeanFraction     = pset.get<float>              ("TruncMeanFraction",                                        0.1);
    fRmsRejectionCutHi     = pset.get<std::vector<double>>("RMSRejectionCutHi",   std::vector<double>() = {25.0,25.0,25.0});
    fRmsRejectionCutLow    = pset.get<std::vector<double>>("RMSRejectionCutLow",  std::vector<double>() = {0.70,0.70,0.70});
    fRmsSelectionCut       = pset.get<std::vector<double>>("RMSSelectionCut",     std::vector<double>() = {1.40,1.40,1.00});
    fTheChosenWire         = pset.get<unsigned int>       ("TheChosenWire",                                           1200);
    fMaxPedestalDiff       = pset.get<double>             ("MaxPedestalDiff",                                          10.);
    fSmoothCorrelatedNoise = pset.get<bool>               ("SmoothCorrelatedNoise",                                   true);
    fNumWiresToGroup       = pset.get<std::vector<size_t>>("NumWiresToGroup",         std::vector<size_t>() = {48, 48, 96});
    fRunFFTInput           = pset.get<bool>               ("RunFFTInputWires",                                       false);
    fRunFFTCorrected       = pset.get<bool>               ("RunFFTCorrectedWires",                                   false);
}

//----------------------------------------------------------------------------
/// Begin job method.
void RawDigitFilterUBooNE::beginJob()
{
    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us.
    art::ServiceHandle<art::TFileService> tfs;
    
    // Define the histograms. Putting semi-colons around the title
    // causes it to be displayed as the x-axis label if the histogram
    // is drawn.
    fAdcCntHist[0]  = tfs->make<TH1D>("CntUPlane", ";#adc",  200, 9000., 10000.);
    fAdcCntHist[1]  = tfs->make<TH1D>("CntVPlane", ";#adc",  200, 9000., 10000.);
    fAdcCntHist[2]  = tfs->make<TH1D>("CntWPlane", ";#adc",  200, 9000., 10000.);
    fAveValHist[0]  = tfs->make<TH1D>("AveUPlane", ";Ave",   120,  -30.,    30.);
    fAveValHist[1]  = tfs->make<TH1D>("AveVPlane", ";Ave",   120,  -30.,    30.);
    fAveValHist[2]  = tfs->make<TH1D>("AveWPlane", ";Ave",   120,  -30.,    30.);
    fRmsValHist[0]  = tfs->make<TH1D>("RmsUPlane", ";RMS",   200,    0.,    50.);
    fRmsValHist[1]  = tfs->make<TH1D>("RmsVPlane", ";RMS",   200,    0.,    50.);
    fRmsValHist[2]  = tfs->make<TH1D>("RmsWPlane", ";RMS",   200,    0.,    50.);
    fPedValHist[0]  = tfs->make<TH1D>("PedUPlane", ";Ped",   200,  1950,  2150.);
    fPedValHist[1]  = tfs->make<TH1D>("PedVPlane", ";Ped",   200,  1950,  2150.);
    fPedValHist[2]  = tfs->make<TH1D>("PedWPlane", ";Ped",   200,   350,   550.);
    
    fRmsValProf[0]  = tfs->make<TProfile>("RmsUPlaneProf",  ";Wire #",  2400, 0., 2400., 0., 100.);
    fRmsValProf[1]  = tfs->make<TProfile>("RmsVPlaneProf",  ";Wire #",  2400, 0., 2400., 0., 100.);
    fRmsValProf[2]  = tfs->make<TProfile>("RmsWPlaneProf",  ";Wire #",  3456, 0., 3456., 0., 100.);

    fPedValProf[0]  = tfs->make<TProfile>("PedUPlaneProf",  ";Wire #",  2400, 0., 2400., 1500., 2500.);
    fPedValProf[1]  = tfs->make<TProfile>("PedVPlaneProf",  ";Wire #",  2400, 0., 2400., 1500., 2500.);
    fPedValProf[2]  = tfs->make<TProfile>("PedWPlaneProf",  ";Wire #",  3456, 0., 3456.,    0., 1000.);
    
    fAverageHist[0] = tfs->make<TH1D>("AverageU", ";Bin", 1000, 1500., 2500.);
    fAverageHist[1] = tfs->make<TH1D>("AverageV", ";Bin", 1000, 1500., 2500.);
    fAverageHist[2] = tfs->make<TH1D>("AverageW", ";Bin", 1000,    0., 1000.);
    
    // Following to determine min/max frequencies
    double sampleRate  = fDetectorProperties->SamplingRate();
    double readOutSize = fDetectorProperties->ReadOutWindowSize();
    double maxFreq     = 1000000. / (2. * sampleRate);
    double minFreq     = 1000000. / (2. * sampleRate * readOutSize);
    int    numSamples  = (readOutSize / 2 + 1) / 4;
    
    fFFTHist[0]     = tfs->make<TProfile>("FFTPlaneU",    "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    fFFTHist[1]     = tfs->make<TProfile>("FFTPlaneV",    "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    fFFTHist[2]     = tfs->make<TProfile>("FFTPlaneW",    "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    
    fFFTHistLow[0]  = tfs->make<TProfile>("FFTLowPlaneU", "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    fFFTHistLow[1]  = tfs->make<TProfile>("FFTLowPlaneV", "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    fFFTHistLow[2]  = tfs->make<TProfile>("FFTLowPlaneW", "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    
    fFFTHistCor[0]  = tfs->make<TProfile>("FFTCorPlaneU", "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    fFFTHistCor[1]  = tfs->make<TProfile>("FFTCorPlaneV", "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    fFFTHistCor[2]  = tfs->make<TProfile>("FFTCorPlaneW", "FFT;kHz",  numSamples, minFreq, maxFreq, 0., 10000.);
    
    fFFTvsMBProf[0] = tfs->make<TProfile2D>("FFTvsMBPlaneU", "FFT;MB;kHz", numSamples, minFreq, maxFreq, 2400/16, 0., 2400/16);
    fFFTvsMBProf[1] = tfs->make<TProfile2D>("FFTvsMBPlaneV", "FFT;MB;kHz", numSamples, minFreq, maxFreq, 2400/16, 0., 2400/16);
    fFFTvsMBProf[2] = tfs->make<TProfile2D>("FFTvsMBPlaneW", "FFT;MB;kHz", numSamples, minFreq, maxFreq, 3456/16, 0., 3456/16);
}

//----------------------------------------------------------------------------
/// Produce method.
///
/// Arguments:
///
/// evt - Art event.
///
/// This is the primary method.
///
void RawDigitFilterUBooNE::produce(art::Event & event)
{
    ++fNumEvent;
    
    // Agreed convention is to ALWAYS output to the event store so get a pointer to our collection
    std::unique_ptr<std::vector<raw::RawDigit> > filteredRawDigit(new std::vector<raw::RawDigit>);
    
    // Read in the digit List object(s).
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
    event.getByLabel(fDigitModuleLabel, digitVecHandle);
    
    // Require a valid handle
    if (digitVecHandle.isValid())
    {
        unsigned int maxChannels    = fGeometry->Nchannels();
        unsigned int maxTimeSamples = fDetectorProperties->NumberTimeSamples();
        double       sampleFreq     = 1000000. / fDetectorProperties->SamplingRate();
        double       readOutSize    = fDetectorProperties->ReadOutWindowSize();
    
        // Ok, to do the correlated noise removal we are going to need a rather impressive data structure...
        // Because we need to unpack each wire's data, we will need to "explode" it out into a data structure
        // here... with the good news that we'll release the memory at the end of the module so should not
        // impact downstream processing (I hope!).
        // What we are going to do is make a vector over views of vectors over wires of vectors over time samples
        std::vector<std::vector<raw::RawDigit::ADCvector_t>> rawDataViewWireTimeVec;
        std::vector<std::vector<float>>                      rawDataViewWireNoiseVec;
        std::vector<std::vector<float>>                      pedestalViewWireVec;
        std::vector<std::vector<float>>                      pedCorViewWireVec;
        std::vector<std::vector<raw::ChannelID_t>>           channelViewWireVec;
    
        // Initialize outer range to number of views
        rawDataViewWireTimeVec.resize(3);
        rawDataViewWireNoiseVec.resize(3);
        pedestalViewWireVec.resize(3);
        pedCorViewWireVec.resize(3);
        channelViewWireVec.resize(3);
    
        // Basic initialization goes here:
        for(size_t viewIdx = 0; viewIdx < 3; viewIdx++)
        {
            // For each view we need to presize the vector to the number of wires
            rawDataViewWireTimeVec[viewIdx].resize(fGeometry->Nwires(viewIdx));
            rawDataViewWireNoiseVec[viewIdx].resize(fGeometry->Nwires(viewIdx));
            pedestalViewWireVec[viewIdx].resize(fGeometry->Nwires(viewIdx));
            pedCorViewWireVec[viewIdx].resize(fGeometry->Nwires(viewIdx));
            channelViewWireVec[viewIdx].resize(fGeometry->Nwires(viewIdx));
        }
    
        // Commence looping over raw digits
        for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter)
        {
            // get the reference to the current raw::RawDigit
            art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
        
            raw::ChannelID_t channel = digitVec->Channel();
        
            bool goodChan(true);
        
            // The below try-catch block may no longer be necessary
            // Decode the channel and make sure we have a valid one
            std::vector<geo::WireID> wids;
            try {
                wids = fGeometry->ChannelToWire(channel);
            }
            catch(...)
            {
                //std::cout << "===>> Found illegal channel with id: " << channel << std::endl;
                goodChan = false;
            }
        
            if (channel >= maxChannels || !goodChan) continue;
        
            // Recover plane and wire in the plane
            unsigned int view = wids[0].Plane;
            unsigned int wire = wids[0].Wire;
        
            unsigned int dataSize = digitVec->Samples();
            
            maxTimeSamples = std::min(maxTimeSamples, dataSize);
        
            // vector holding uncompressed adc values
            std::vector<short>& rawadc = rawDataViewWireTimeVec[view][wire];
        
            channelViewWireVec[view][wire] = channel;
        
            rawadc.resize(maxTimeSamples);
        
            // And now uncompress
            raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());

            // The strategy for finding the average for a given wire will be to
            // find the most populated bin and the average using the neighboring bins
            // To do this we'll use a map with key the bin number and data the count in that bin
            // Define the map first
            std::map<short,short> binAdcMap;
        
            // Populate the mape
            for(const auto& adcVal : rawadc)
            {
                binAdcMap[adcVal]++;
            }
        
            // Find the max bin
            short binMax(-1);
            short binMaxCnt(0);
        
            for(const auto& binAdcItr : binAdcMap)
            {
                if (binAdcItr.second > binMaxCnt)
                {
                    binMax    = binAdcItr.first;
                    binMaxCnt = binAdcItr.second;
                }
            }
        
            // fill example hists - throw away code
            if (fFirstEvent && wire == fTheChosenWire)
            {
                for(const auto& binAdcItr : binAdcMap)
                {
                    fAverageHist[view]->Fill(binAdcItr.first, binAdcItr.second);
                }
            }
        
            // Armed with the max bin and its count, now set up to get an average
            // about this bin. We'll want to cut off at some user defined fraction
            // of the total bins on the wire
            int minNumBins = (1. - fTruncMeanFraction) * dataSize - 1;
            int curBinCnt(binMaxCnt);
        
            double peakValue(curBinCnt * binMax);
            double truncMean(peakValue);
        
            short binOffset(1);
        
            // This loop to develop the average
            // In theory, we could also keep the sum of the squares for the rms but I had problems doing
            // it that way so will loop twice... (potential time savings goes here!)
            while(curBinCnt < minNumBins)
            {
                if (binAdcMap[binMax-binOffset])
                {
                    curBinCnt += binAdcMap[binMax-binOffset];
                    truncMean += double(binAdcMap[binMax-binOffset] * (binMax - binOffset));
                }
            
                if (binAdcMap[binMax+binOffset])
                {
                    curBinCnt += binAdcMap[binMax+binOffset];
                    truncMean += double(binAdcMap[binMax+binOffset] * (binMax + binOffset));
                }
            
                binOffset++;
            }
        
            truncMean /= double(curBinCnt);

            // do rms calculation - the old fashioned way and over all adc values
            double rmsVal = 0.;
            
            for(const auto& adcVal : rawadc)
            {
                double adcLessPed = adcVal - truncMean;
                rmsVal += adcLessPed * adcLessPed;
            }
            
            rmsVal = std::sqrt(std::max(0.,rmsVal / double(rawadc.size())));
        
            rawDataViewWireNoiseVec[view][wire] = rmsVal;

            // Recover the database version of the pedestal
            float pedestal = fPedestalRetrievalAlg.PedMean(channel);
        
            pedestalViewWireVec[view][wire] = truncMean; //pedestal;
            pedCorViewWireVec[view][wire]   = truncMean - pedestal;
        
            // Fill some histograms here
            fAdcCntHist[view]->Fill(curBinCnt, 1.);
            fAveValHist[view]->Fill(std::max(-29.9, std::min(29.9,truncMean - pedestal)), 1.);
            fRmsValHist[view]->Fill(std::min(49.9, rmsVal), 1.);
            fRmsValProf[view]->Fill(wire, rmsVal, 1.);
            fPedValProf[view]->Fill(wire, truncMean, 1.);
            fPedValHist[view]->Fill(truncMean, 1.);

            // Output a message is there is significant different to the pedestal
            if (abs(truncMean - pedestal) > fMaxPedestalDiff)
            {
                mf::LogInfo("RawDigitFilterUBooNE") << ">>> Pedestal mismatch, channel: " << channel << ", new value: " << truncMean << ", original: " << pedestal << ", rms: " << rmsVal << std::endl;
            }
            
            // Keep the RawDigit if below our rejection cut
            if (rmsVal < fRmsRejectionCutHi[view])
            {
                if (!fSmoothCorrelatedNoise) filteredRawDigit->emplace_back(*digitVec);
            }
            else
            {
                rawadc.clear();
                
                mf::LogInfo("RawDigitFilterUBooNE") <<  "--> Rejecting channel for large rms, channel: " << channel << ", rmsVal: " << rmsVal << ", truncMean: " << truncMean << ", pedestal: " << pedestal << std::endl;
            }

            // Plug in FFT here
            if (fRunFFTInput && rmsVal < fRmsRejectionCutHi[view])
            {
                int fftDataSize = rawadc.size();
            
                TVirtualFFT* fftr2c = TVirtualFFT::FFT(1, &fftDataSize, "R2C");
            
                double fftInputArray[fftDataSize];
            
                for(size_t idx = 0; idx < size_t(fftDataSize); idx++) fftInputArray[idx] = rawadc[idx] - truncMean;
            
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

                    if (rmsVal > fRmsRejectionCutLow[view])
                    {
                        int mbIdx = wire/ 16;
                        fFFTHist[view]->Fill(bin, power, 1.);
                        fFFTvsMBProf[view]->Fill(bin, mbIdx, power);
                    }
                    else fFFTHistLow[view]->Fill(bin, power, 1.);
                }
            }
        }
    
        // Try to implement Corey's algorithm here
        // The basic idea is to try to take groups of wires and find a metric within a given time bin
        // to use to correct the adc values on the wire
        // Make sure we want to do this...
        if (fSmoothCorrelatedNoise)
        {
            // Perform the outer loop over views
            for(size_t viewIdx = 0; viewIdx < 3; viewIdx++)
            {
                size_t nWiresPerMotherBoard(fNumWiresToGroup[viewIdx]);
                
                // How many groups of wires this view?
                size_t nMotherBoards = fGeometry->Nwires(viewIdx) / nWiresPerMotherBoard;
                
                // Get vectors for this view
                std::vector<raw::RawDigit::ADCvector_t>& rawDataWireTimeVec  = rawDataViewWireTimeVec[viewIdx];
                std::vector<float>&                      pedestalWireVec     = pedestalViewWireVec[viewIdx];
                std::vector<float>&                      pedCorWireVec       = pedCorViewWireVec[viewIdx];
                std::vector<float>&                      rawDataWireNoiseVec = rawDataViewWireNoiseVec[viewIdx];
        
                // Loop over wires in group (probably a motherboard's worth)
                for(size_t mbIdx = 0; mbIdx < nMotherBoards; mbIdx++)
                {
                    // Get wire offset this section
                    size_t wireBaseOffset = mbIdx * nWiresPerMotherBoard;
                    
                    // Try to optimize the loop a bit by pre-checking if a wire is good or bad
                    std::map<size_t, raw::RawDigit::ADCvector_t&> wireToAdcMap;
                    
                    // Keep track of gap sizes
                    size_t lastWireIdx(wireBaseOffset);
                    size_t largestGapSize(0);
                    
                    // Finally, inside of here we are looping over wires on a motherboard
                    for(size_t wireIdx = 0; wireIdx < nWiresPerMotherBoard; wireIdx++)
                    {
                        // Recover the physical wire
                        size_t physWireIdx = wireBaseOffset + wireIdx;
                        
                        // If the channel has been marked bad we simply ignore
                        if (rawDataWireTimeVec[physWireIdx].empty()) continue;
                        
                        // If this wire is too noisy, or not enough noisy, reject
                        double rmsNoise(rawDataWireNoiseVec[physWireIdx]);
                        
                        // We don't want to adjust low noise channels but do want them in the output
                        if (rmsNoise < fRmsSelectionCut[viewIdx])
                        {
                            filteredRawDigit->emplace_back(raw::RawDigit(channelViewWireVec[viewIdx][physWireIdx], maxTimeSamples, rawDataWireTimeVec[physWireIdx], raw::kNone));
                            filteredRawDigit->back().SetPedestal(pedestalWireVec[physWireIdx],rmsNoise);
                            continue;
                        }
                        
                        // Don't select "bad" wires, they are lost anyway
                        // Also, it is pointless to try to include the "ultra low noise" channels (or correct them)
                        if (rmsNoise < fRmsRejectionCutHi[viewIdx])
                        {
                            wireToAdcMap.insert(std::pair<size_t,raw::RawDigit::ADCvector_t&>(physWireIdx,rawDataWireTimeVec[physWireIdx]));
                            
                            if (physWireIdx - lastWireIdx > largestGapSize)
                            {
                                largestGapSize = physWireIdx - lastWireIdx;
                            }
                            
                            lastWireIdx = physWireIdx;
                        }
                    }
                    
                    // Don't try to do correction if too few wires unless they have gaps
                    if (wireToAdcMap.size() > 5 || largestGapSize > 2)
                    {
                        // Now we loop over the number of time bins (samples)
                        for(size_t sampleIdx = 0; sampleIdx < maxTimeSamples; sampleIdx++)
                        {
                            // Define a vector for accumulating values...
                            // Loop over the wires at this time bin and get their pedestal corrected ADC values
                            // We'll use a simple stl vector for this
                            std::vector<float> adcValuesVec;
                            
                            for(const auto& wireAdcItr : wireToAdcMap)
                                adcValuesVec.push_back(float(wireAdcItr.second[sampleIdx]) - pedestalWireVec[wireAdcItr.first]);
                            
                            // Sort to get the right order...
                            std::sort(adcValuesVec.begin(),adcValuesVec.end());
                            
                            // Extract the median
                            size_t medianIdx   = adcValuesVec.size() / 2;
                            float  medianValue = adcValuesVec[medianIdx];
                            
                            // What if we have an odd number of wires we're working with?
                            if (adcValuesVec.size() % 2)
                                medianValue = 0.5 * (medianValue + adcValuesVec[medianIdx-1]);
                            
                            //float mostProbableValue = 0.1 * float(maxAdcValue);
                            float mostProbableValue = medianValue; //float(maxAdcValue);

                            // Now run through and apply correction
                            for (const auto& wireAdcItr : wireToAdcMap)
                            {
                                // Probably doesn't matter, but try to get slightly more accuracy by doing float math and rounding
                                float     newAdcValueFloat = float(wireAdcItr.second[sampleIdx]) - mostProbableValue - pedCorWireVec[wireAdcItr.first];
                                short int newAdcValue      = std::round(newAdcValueFloat);
                            
                                wireAdcItr.second[sampleIdx] = newAdcValue;
                            }
                        }
                    }

                    // One more pass through to store the good channels
                    for (const auto& wireAdcItr : wireToAdcMap)
                    {
                        // recalculate rms
                        double rmsVal   = 0.;
                        double pedestal = pedestalWireVec[wireAdcItr.first];
                        
                        for(const auto& adcVal : wireAdcItr.second)
                        {
                            double adcLessPed = adcVal - pedestal;
                            rmsVal += adcLessPed * adcLessPed;
                        }
                        
                        rmsVal = std::sqrt(std::max(0.,rmsVal / double(wireAdcItr.second.size())));
                        
                        filteredRawDigit->emplace_back(raw::RawDigit(channelViewWireVec[viewIdx][wireAdcItr.first], maxTimeSamples, wireAdcItr.second, raw::kNone));
                        filteredRawDigit->back().SetPedestal(pedestal,rmsVal);
                        
                        // Run an FFT here to check our "corrected" wires
                        if (fRunFFTCorrected)
                        {
                            int fftDataSize = wireAdcItr.second.size();
                            
                            TVirtualFFT* fftr2c = TVirtualFFT::FFT(1, &fftDataSize, "R2C");
                            
                            double fftInputArray[fftDataSize];
                            
                            for(size_t idx = 0; idx < size_t(fftDataSize); idx++) fftInputArray[idx] = wireAdcItr.second[idx] - pedestal;
                            
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
                }
            }
        }
    }
    
    // Reset this silly flag so we only fill our example hists once...
    fFirstEvent = false;
    
    // Add tracks and associations to event.
    event.put(std::move(filteredRawDigit));
}

//----------------------------------------------------------------------------
/// End job method.
void RawDigitFilterUBooNE::endJob()
{
    mf::LogInfo("RawDigitFilterUBooNE") << "Looked at " << fNumEvent << " events" << std::endl;
}
