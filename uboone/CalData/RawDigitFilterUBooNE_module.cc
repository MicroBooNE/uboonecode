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
// FillHistograms        - Turn on histogram filling for diagnostics
// RunFFTInputWires      - FFT analyze the input RawDigits if true - diagnostics
// RunFFTCorrectedWires  - FFT analyze the output RawDigits if true - diagnostics
// TruncateTicks:        - Provide mechanism to truncate a readout window to a smaller size
// WindowSize:           - The desired size of the output window
// NumTicksToDropFront:  - The number ticks dropped off the front of the original window
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
#include "art/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"

#include "lardata/RawData/RawDigit.h"
#include "lardata/RawData/raw.h"

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
    
    // Basic waveform mean, rms and pedestal offset
    void getMeanAndRms(const raw::RawDigit::ADCvector_t& rawWaveform,
                       unsigned int                      channel,
                       unsigned int                      view,
                       unsigned int                      wire,
                       float&                            aveVal,
                       float&                            rmsVal,
                       float&                            pedCorVal);
    
    void removeCorrelatedNoise(std::vector<raw::RawDigit::ADCvector_t>&    rawDataWireTimeVec,
                               unsigned int                                viewIdx,
                               unsigned int                                wireBaseOffset,
                               std::vector<float>&                         pedestalWireVec,
                               std::vector<float>&                         pedCorWireVec,
                               std::vector<float>&                         rawDataWireNoiseVec,
                               std::map<int, raw::RawDigit::ADCvector_t&>& wireToAdcMap);
    
    void saveRawDigits(std::unique_ptr<std::vector<raw::RawDigit> >&, raw::ChannelID_t&, raw::RawDigit::ADCvector_t&, float, float);

    template<class T> T getMedian(std::vector<T>&, T) const;
    
    // Fcl parameters.
    std::string          fDigitModuleLabel;      ///< The full collection of hits
    float                fTruncMeanFraction;     ///< Fraction for truncated mean
    std::vector<float>   fRmsRejectionCutHi;     ///< Maximum rms for input channels, reject if larger
    std::vector<float>   fRmsRejectionCutLow;    ///< Minimum rms to consider channel "alive"
    std::vector<float>   fRmsSelectionCut;       ///< Don't use/apply correction to wires below this
    unsigned int         fTheChosenWire;         ///< For example hist
    double               fMaxPedestalDiff;       ///< Max pedestal diff to db to warn
    bool                 fSmoothCorrelatedNoise; ///< Should we smooth the noise?
    std::vector<size_t>  fNumWiresToGroup;       ///< If smoothing, the number of wires to look at
    bool                 fFillHistograms;        ///< if true then will fill diagnostic hists
    bool                 fRunFFTInput;           ///< Should we run FFT's on input wires?
    bool                 fRunFFTCorrected;       ///< Should we run FFT's on corrected wires?
    bool                 fTruncateTicks;         ///< If true then drop channels off ends of wires
    unsigned int         fWindowSize;            ///< # ticks to keep in window
    unsigned int         fNumTicksToDropFront;   ///< # ticks to drop from front of waveform

    // Statistics.
    int fNumEvent;        ///< Number of events seen.
    
    // Pointers to the histograms we'll create for monitoring what is happening
    TH1D*                fAdcCntHist[3];
    TH1D*                fAveValHist[3];
    TH1D*                fRmsValHist[3];
    TH1D*                fPedValHist[3];
    TH1D*                fAverageHist[3];
    TProfile*            fRmsValProf[3];
    TProfile*            fPedValProf[3];

    TProfile*            fFFTHist[3];
    TProfile*            fFFTHistLow[3];
    TProfile*            fFFTHistCor[3];
    
    TProfile2D*          fFFTvsMBProf[3];
    
    bool                 fFirstEvent;
    
    // Useful services, keep copies for now (we can update during begin run periods)
    geo::GeometryCore const* fGeometry;                         ///< pointer to Geometry service
    detinfo::DetectorProperties const* fDetectorProperties;   ///< Detector properties service
    const lariov::DetPedestalProvider&          fPedestalRetrievalAlg; ///< Keep track of an instance to the pedestal retrieval alg
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
                      fPedestalRetrievalAlg(*lar::providerFrom<lariov::DetPedestalService>())
{
    
    fGeometry = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
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
    fDigitModuleLabel      = pset.get<std::string>        ("DigitModuleLabel",                                      "daq");
    fTruncMeanFraction     = pset.get<float>              ("TruncMeanFraction",                                       0.1);
    fRmsRejectionCutHi     = pset.get<std::vector<float>> ("RMSRejectionCutHi",   std::vector<float>() = {25.0,25.0,25.0});
    fRmsRejectionCutLow    = pset.get<std::vector<float>> ("RMSRejectionCutLow",  std::vector<float>() = {0.70,0.70,0.70});
    fRmsSelectionCut       = pset.get<std::vector<float>> ("RMSSelectionCut",     std::vector<float>() = {1.40,1.40,1.00});
    fTheChosenWire         = pset.get<unsigned int>       ("TheChosenWire",                                          1200);
    fMaxPedestalDiff       = pset.get<double>             ("MaxPedestalDiff",                                         10.);
    fSmoothCorrelatedNoise = pset.get<bool>               ("SmoothCorrelatedNoise",                                  true);
    fNumWiresToGroup       = pset.get<std::vector<size_t>>("NumWiresToGroup",        std::vector<size_t>() = {48, 48, 96});
    fFillHistograms        = pset.get<bool>               ("FillHistograms",                                        false);
    fRunFFTInput           = pset.get<bool>               ("RunFFTInputWires",                                      false);
    fRunFFTCorrected       = pset.get<bool>               ("RunFFTCorrectedWires",                                  false);
    fTruncateTicks         = pset.get<bool>               ("TruncateTicks",                                         false);
    fWindowSize            = pset.get<size_t>             ("WindowSize",                                             6400);
    fNumTicksToDropFront   = pset.get<size_t>             ("NumTicksToDropFront",                                    2250);
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
        
        // Sadly, the RawDigits come to us in an unsorted condition which is not optimal for
        // what we want to do here. So we make a vector of pointers to the input raw digits and sort them
        std::vector<const raw::RawDigit*> rawDigitVec;
        
        // Ugliness to fill the pointer vector...
        for(size_t idx = 0; idx < digitVecHandle->size(); idx++) rawDigitVec.push_back(&digitVecHandle->at(idx)); //art::Ptr<raw::RawDigit>(digitVecHandle, idx).get());
        
        // Sort (use a lambda to sort by channel id)
        std::sort(rawDigitVec.begin(),rawDigitVec.end(),[](const raw::RawDigit* left, const raw::RawDigit* right) {return left->Channel() < right->Channel();});
    
        // Ok, to do the correlated noise removal we are going to need a rather impressive data structure...
        // Because we need to unpack each wire's data, we will need to "explode" it out into a data structure
        // here... with the good news that we'll release the memory at the end of the module so should not
        // impact downstream processing (I hope!).
        // What we are going to do is make a vector over views of vectors over wires of vectors over time samples
        std::vector<raw::RawDigit::ADCvector_t> rawDataWireTimeVec;
        std::vector<float>                      rawDataWireNoiseVec;
        std::vector<float>                      pedestalWireVec;
        std::vector<float>                      pedCorWireVec;
        std::vector<raw::ChannelID_t>           channelWireVec;
        
        // Declare a temporary digit holder and resize it if downsizing the waveform
        raw::RawDigit::ADCvector_t tempVec;
        if (fTruncateTicks) tempVec.resize(9595);
    
        // Commence looping over raw digits
        for(const auto& rawDigit : rawDigitVec)
        {
            raw::ChannelID_t channel = rawDigit->Channel();
        
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
        
            unsigned int dataSize = rawDigit->Samples();
            unsigned int wireIdx  = wire % fNumWiresToGroup[view];
            
            // Cross check that our storage arrays are the correct size
            // (note there is a possible boundary issue here that we are going to ignore...)
            if (rawDataWireTimeVec.size() != fNumWiresToGroup[view])
            {
                // For each view we need to presize the vector to the number of wires
                rawDataWireTimeVec.resize(fNumWiresToGroup[view]);
                rawDataWireNoiseVec.resize(fNumWiresToGroup[view]);
                pedestalWireVec.resize(fNumWiresToGroup[view]);
                pedCorWireVec.resize(fNumWiresToGroup[view]);
                channelWireVec.resize(fNumWiresToGroup[view]);
            }
        
            // vector holding uncompressed adc values
            std::vector<short>& rawadc = rawDataWireTimeVec[wireIdx];
        
            channelWireVec[wireIdx] = channel;
            
            // If we are trying to truncate the incoming RawDigit collection then we need to do so when we extract from the input raw digits
            // This causes a small division here...
            if (fTruncateTicks)
            {
                maxTimeSamples = fWindowSize;
                
                if (rawadc.size() != maxTimeSamples) rawadc.resize(maxTimeSamples);
                
                // And now uncompress
                raw::Uncompress(rawDigit->ADCs(), tempVec, rawDigit->Compression());
                
                std::copy(tempVec.begin() + fNumTicksToDropFront, tempVec.begin() + fNumTicksToDropFront + fWindowSize, rawadc.begin());
            }
            else
            {
                maxTimeSamples = dataSize;
                
                if (rawadc.size() != dataSize) rawadc.resize(maxTimeSamples);
                
                // And now uncompress
                raw::Uncompress(rawDigit->ADCs(), rawadc, rawDigit->Compression());
            }
            
            // do the mean and rms calculation
            getMeanAndRms(rawadc, channel, view, wireIdx, pedestalWireVec[wireIdx], rawDataWireNoiseVec[wireIdx], pedCorWireVec[wireIdx]);
            
            // If we are not performing noise corrections then we are done with this wire
            // Store it and move on
            if (!fSmoothCorrelatedNoise)
            {
                // Filter out the very high noise wires
                if (rawDataWireNoiseVec[wireIdx] < fRmsRejectionCutHi[view])
                    saveRawDigits(filteredRawDigit, channel, rawadc, pedestalWireVec[wireIdx], rawDataWireNoiseVec[wireIdx]);
                else
                {
                    // Eventually we'll interface to some sort of channel status communication mechanism.
                    // For now use the log file
                    mf::LogInfo("RawDigitFilterUBooNE") <<  "--> Rejecting channel for large rms, channel: " << channel << ", rmsVal: " << pedCorWireVec[wireIdx] << ", truncMean: " << pedestalWireVec[wireIdx] << ", pedestal: " << pedCorWireVec[wireIdx] << std::endl;
                }
                
                continue;
            }

            // Are we at the correct boundary for dealing with the noise?
            if (!((wireIdx + 1) % fNumWiresToGroup[view]))
            {
                std::map<int, raw::RawDigit::ADCvector_t&> wireToAdcMap;
                size_t baseWireIdx = wire - wire % fNumWiresToGroup[view];
                    
                removeCorrelatedNoise(rawDataWireTimeVec,
                                      view,
                                      baseWireIdx,
                                      pedestalWireVec,
                                      pedCorWireVec,
                                      rawDataWireNoiseVec,
                                      wireToAdcMap);
                    
                // One more pass through to store the good channels
                for (const auto& wireAdcItr : wireToAdcMap)
                {
                    size_t wireIdx = abs(wireAdcItr.first) % fNumWiresToGroup[view];
                        
                    // recalculate rms for the output
                    double rmsVal   = 0.;
                    double pedestal = pedestalWireVec[wireIdx];
                    double pedCor   = pedCorWireVec[wireIdx];
                        
                    for(const auto& adcVal : wireAdcItr.second)
                    {
                        double adcLessPed = adcVal - pedestal + pedCor;
                        rmsVal += adcLessPed * adcLessPed;
                    }
                        
                    rmsVal = std::sqrt(std::max(0.,rmsVal / double(wireAdcItr.second.size())));
                    
                    // The ultra high noise channels are simply zapped
                    if (rmsVal < fRmsRejectionCutHi[view])
                    {
                        saveRawDigits(filteredRawDigit, channelWireVec[wireIdx], wireAdcItr.second, pedestal, rmsVal);
                    }
                    else
                    {
                        mf::LogInfo("RawDigitFilterUBooNE") <<  "--> Rejecting channel for large rms, channel: " << channel << ", rmsVal: " << rawDataWireNoiseVec[wireIdx] << ", truncMean: " << pedestalWireVec[wireIdx] << ", pedestal: " << pedCorWireVec[wireIdx] << std::endl;
                    }
                }
            }

            // Plug in FFT here
            if (fRunFFTInput && pedCorWireVec[wireIdx] < fRmsRejectionCutHi[view])
            {
                int fftDataSize = rawadc.size();
            
                TVirtualFFT* fftr2c = TVirtualFFT::FFT(1, &fftDataSize, "R2C");
            
                double fftInputArray[fftDataSize];
            
                for(size_t idx = 0; idx < size_t(fftDataSize); idx++) fftInputArray[idx] = rawadc[idx] - pedestalWireVec[wireIdx];
            
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

                    if (pedCorWireVec[wire] > fRmsRejectionCutLow[view])
                    {
                        int mbIdx = wire/ 16;
                        fFFTHist[view]->Fill(bin, power, 1.);
                        fFFTvsMBProf[view]->Fill(bin, mbIdx, power);
                    }
                    else fFFTHistLow[view]->Fill(bin, power, 1.);
                }
            }
        }
    }
    
    // Reset this silly flag so we only fill our example hists once...
    fFirstEvent = false;
    
    // Add tracks and associations to event.
    event.put(std::move(filteredRawDigit));
}

void RawDigitFilterUBooNE::getMeanAndRms(const raw::RawDigit::ADCvector_t& rawWaveform,
                                         unsigned int                      channel,
                                         unsigned int                      view,
                                         unsigned int                      wire,
                                         float&                            aveVal,
                                         float&                            rmsVal,
                                         float&                            pedCorVal)
{
    // The strategy for finding the average for a given wire will be to
    // find the most populated bin and the average using the neighboring bins
    // To do this we'll use a map with key the bin number and data the count in that bin
    // Define the map first
    std::map<short,short> binAdcMap;
    
    // Populate the map
    for(const auto& adcVal : rawWaveform)
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
    if (fFillHistograms && fFirstEvent && wire == fTheChosenWire)
    {
        for(const auto& binAdcItr : binAdcMap)
        {
            fAverageHist[view]->Fill(binAdcItr.first, binAdcItr.second);
        }
    }
    
    // Armed with the max bin and its count, now set up to get an average
    // about this bin. We'll want to cut off at some user defined fraction
    // of the total bins on the wire
    size_t maxTimeSamples(rawWaveform.size());
    int    minNumBins = (1. - fTruncMeanFraction) * maxTimeSamples - 1;
    int    curBinCnt(binMaxCnt);
    
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
    rmsVal = 0.;
    
    for(const auto& adcVal : rawWaveform)
    {
        double adcLessPed = adcVal - truncMean;
        rmsVal += adcLessPed * adcLessPed;
    }
    
    rmsVal = std::sqrt(std::max(0.,rmsVal / double(maxTimeSamples)));
    
    // Recover the database version of the pedestal
    float pedestal = fPedestalRetrievalAlg.PedMean(channel);
    
    aveVal    = truncMean;
    pedCorVal = truncMean - pedestal;
    
    // Fill some histograms here
    if (fFillHistograms)
    {
        fAdcCntHist[view]->Fill(curBinCnt, 1.);
        fAveValHist[view]->Fill(std::max(-29.9, std::min(29.9,truncMean - pedestal)), 1.);
        fRmsValHist[view]->Fill(std::min(49.9, double(rmsVal)), 1.);
        fRmsValProf[view]->Fill(wire, double(rmsVal), 1.);
        fPedValProf[view]->Fill(wire, truncMean, 1.);
        fPedValHist[view]->Fill(truncMean, 1.);
    }
    
    // Output a message is there is significant different to the pedestal
    if (abs(truncMean - pedestal) > fMaxPedestalDiff)
    {
        mf::LogInfo("RawDigitFilterUBooNE") << ">>> Pedestal mismatch, channel: " << channel << ", new value: " << truncMean << ", original: " << pedestal << ", rms: " << rmsVal << std::endl;
    }
    
    return;
}

void RawDigitFilterUBooNE::removeCorrelatedNoise(std::vector<raw::RawDigit::ADCvector_t>&    rawDataWireTimeVec,
                                                 unsigned int                                viewIdx,
                                                 unsigned int                                wireBaseOffset,
                                                 std::vector<float>&                         pedestalWireVec,
                                                 std::vector<float>&                         pedCorWireVec,
                                                 std::vector<float>&                         rawDataWireNoiseVec,
                                                 std::map<int, raw::RawDigit::ADCvector_t&>& wireToAdcMap)
{
    // This is a straightforward implementation of "Corey's Algorithm" to reduce the noise that is
    // correlated across channels.
    size_t lastWireIdx(wireBaseOffset);
    size_t largestGapSize(0);
    size_t maxTimeSamples(rawDataWireTimeVec[0].size());
    
    // What is median rms of wires input?
    std::vector<float> rmsValsVec;
    
    for(size_t wireIdx = 0; wireIdx < fNumWiresToGroup[viewIdx]; wireIdx++)
    {
        float rmsNoise(rawDataWireNoiseVec[wireIdx]);
        
        if (rmsNoise > fRmsSelectionCut[viewIdx] && rmsNoise < fRmsRejectionCutHi[viewIdx])
            rmsValsVec.push_back(rmsNoise);
    }
    
    float meanRms = getMedian(rmsValsVec, fRmsSelectionCut[viewIdx]);
    
    // Finally, inside of here we are looping over wires on a motherboard
    for(size_t wireIdx = 0; wireIdx < fNumWiresToGroup[viewIdx]; wireIdx++)
    {
        // Recover the physical wire
        size_t physWireIdx = wireBaseOffset + wireIdx;
        
        // If the channel has been marked bad we simply ignore
        if (rawDataWireTimeVec[wireIdx].empty()) continue;
        
        // If this wire is too noisy, or not enough noisy, reject
        double rmsNoise(rawDataWireNoiseVec[wireIdx]);
        
        // We don't want to adjust low noise channels or high noise channels but do want them in the output
        if (rmsNoise < fRmsSelectionCut[viewIdx] || rmsNoise < 0.5 * meanRms || rmsNoise > fRmsRejectionCutHi[viewIdx])
        {
            wireToAdcMap.insert(std::pair<int,raw::RawDigit::ADCvector_t&>(-physWireIdx,rawDataWireTimeVec[wireIdx]));

            continue;
        }
        
        // Don't select "bad" wires, they are lost anyway
        // Also, it is pointless to try to include the "ultra low noise" channels (or correct them)
        wireToAdcMap.insert(std::pair<int,raw::RawDigit::ADCvector_t&>(physWireIdx,rawDataWireTimeVec[wireIdx]));
        
        if (physWireIdx - lastWireIdx > largestGapSize)
        {
            largestGapSize = physWireIdx - lastWireIdx;
        }
        
        lastWireIdx = physWireIdx;
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
            {
                // A negative physical wire idx means we are meant to skip this wire in the averaging processing
                if (wireAdcItr.first < 0) continue;
                
                // Accumulate
                adcValuesVec.push_back(float(wireAdcItr.second[sampleIdx]) - pedestalWireVec[wireAdcItr.first%fNumWiresToGroup[viewIdx]]);
            }
            
            float medianValue = getMedian(adcValuesVec, float(0.));
            
            // Now run through and apply correction
            for (const auto& wireAdcItr : wireToAdcMap)
            {
                float corVal(medianValue);
                int   wireIdx(wireAdcItr.first);
                
                // A negative wire means we are meant to skip this wire in the averaging processing
                if (wireIdx < 0)
                {
                    corVal  = 0.;
                    wireIdx = -wireIdx;
                }
                
                // Probably doesn't matter, but try to get slightly more accuracy by doing float math and rounding
                float     newAdcValueFloat = float(wireAdcItr.second[sampleIdx]) - corVal - pedCorWireVec[wireIdx%fNumWiresToGroup[viewIdx]];
                short int newAdcValue      = std::round(newAdcValueFloat);
                
                wireAdcItr.second[sampleIdx] = newAdcValue;
            }
        }
    }
    
    // Run an FFT here to check our "corrected" wires
    if (fRunFFTCorrected)
    {
        double       sampleFreq     = 1000000. / fDetectorProperties->SamplingRate();
        double       readOutSize    = fDetectorProperties->ReadOutWindowSize();
        
        for(const auto& wireAdcItr : wireToAdcMap)
        {
            if (wireAdcItr.first < 0) continue;
            
            size_t wireIdx = wireAdcItr.first % fNumWiresToGroup[viewIdx];
            
            int fftDataSize = wireAdcItr.second.size();
        
            TVirtualFFT* fftr2c = TVirtualFFT::FFT(1, &fftDataSize, "R2C");
        
            double fftInputArray[fftDataSize];
        
            for(size_t idx = 0; idx < size_t(fftDataSize); idx++) fftInputArray[idx] = wireAdcItr.second[idx] - pedestalWireVec[wireIdx];
        
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

template<class T> T RawDigitFilterUBooNE::getMedian(std::vector<T>& valuesVec, T defaultValue) const
{
    T medianValue(defaultValue);
    
    if (!valuesVec.empty())
    {
        std::sort(valuesVec.begin(),valuesVec.end());
        
        size_t medianIdx = valuesVec.size() / 2;
        
        medianValue = valuesVec[medianIdx];
        
        if (valuesVec.size() > 1 && medianIdx % 2) medianValue = (medianValue + valuesVec[medianIdx+1]) / 2;
    }
    
    return medianValue;
}

void RawDigitFilterUBooNE::saveRawDigits(std::unique_ptr<std::vector<raw::RawDigit> >& filteredRawDigit,
                                         raw::ChannelID_t&                             channel,
                                         raw::RawDigit::ADCvector_t&                   rawDigitVec,
                                         float                                         pedestal,
                                         float                                         rms)
{
    filteredRawDigit->emplace_back(raw::RawDigit(channel, rawDigitVec.size(), rawDigitVec, raw::kNone));
    filteredRawDigit->back().SetPedestal(pedestal,rms);
    
    return;
}

//----------------------------------------------------------------------------
/// End job method.
void RawDigitFilterUBooNE::endJob()
{
    mf::LogInfo("RawDigitFilterUBooNE") << "Looked at " << fNumEvent << " events" << std::endl;
}
