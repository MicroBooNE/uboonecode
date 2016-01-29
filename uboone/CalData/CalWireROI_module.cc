////////////////////////////////////////////////////////////////////////
//
// CalWireROI class - variant of CalWire that deconvolves in 
// Regions Of Interest
//
// baller@fnal.gov
//
//
////////////////////////////////////////////////////////////////////////

// C/C++ standard libraries
#include <string>
#include <vector>
#include <utility> // std::pair<>
#include <memory> // std::unique_ptr<>
#include <iomanip>

// ROOT libraries
#include "TComplex.h"
#include "TH1D.h"

// framework libraries
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Utilities/Exception.h"

// LArSoft libraries
#include "SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "Geometry/Geometry.h"
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "RecoBase/Wire.h"
#include "RecoBaseArt/WireCreator.h"
#include "Utilities/LArFFT.h"
#include "Utilities/AssociationUtil.h"
#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"
#include "CalibrationDBI/Interface/IDetPedestalService.h"
#include "CalibrationDBI/Interface/IDetPedestalProvider.h"
#include "CalibrationDBI/Interface/IChannelStatusService.h"
#include "CalibrationDBI/Interface/IChannelStatusProvider.h"

#include "WaveformPropertiesAlg.h"

///creation of calibrated signals on wires
namespace caldata {

class CalWireROI : public art::EDProducer
{
  public:
    // create calibrated signals on wires. this class runs 
    // an fft to remove the electronics shaping.     
    explicit CalWireROI(fhicl::ParameterSet const& pset); 
    virtual ~CalWireROI();
    
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob();                 
    void reconfigure(fhicl::ParameterSet const& p);
    
    void reconfFFT(int temp_fftsize);
    
  private:
    
    std::string                 fDigitModuleLabel;     ///< module that made digits
    std::string                 fSpillName;            ///< nominal spill is an empty string
                                                       ///< it is set by the DigitModuleLabel
                                                       ///< ex.:  "daq:preSpill" for prespill data
    unsigned short              fNoiseSource;          ///< Used to determine ROI threshold
    unsigned short              fNumBinsHalf;          ///< Determines # bins in ROI running sum
    std::vector<unsigned short> fThreshold;            ///< abs(threshold) ADC counts for ROI
    std::vector<int>            fNumSigma;             ///< "# sigma" rms noise for ROI threshold
    int                         fFFTSize;              ///< FFT size for ROI deconvolution
    std::vector<unsigned short> fPreROIPad;            ///< ROI padding
    std::vector<unsigned short> fPostROIPad;           ///< ROI padding
    bool                        fDoBaselineSub;        ///< Do baseline subtraction after deconvolution?
    bool                        fuPlaneRamp;           ///< set true for correct U plane wire response
    int                         fSaveWireWF;           ///< Save recob::wire object waveforms
    size_t                      fEventCount;           ///< count of event processed
    int                         fMinAllowedChanStatus; ///< Don't consider channels with lower status
    
    void doDecon(std::vector<float>&                                       holder,
                 raw::ChannelID_t                                          channel,
                 unsigned int                                              thePlane,
                 const std::vector<std::pair<size_t, size_t>>&             rois,
                 const std::vector<std::pair<size_t, size_t>>&             holderInfo,
                 recob::Wire::RegionsOfInterest_t&                         ROIVec,
                 art::ServiceHandle<util::SignalShapingServiceMicroBooNE>& sss);
    
    float SubtractBaseline(std::vector<float>& holder,
                           float               basePre,
                           float               basePost,
                           size_t              roiStart,
                           size_t              roiLen,
                           size_t              dataSize);

    bool                               fDoBaselineSub_WaveformPropertiesAlg;
    util::WaveformPropertiesAlg<float> fROIPropertiesAlg;
    float SubtractBaseline(const std::vector<float>& holder);
    
  protected: 
    
  }; // class CalWireROI

  DEFINE_ART_MODULE(CalWireROI)
  
  //-------------------------------------------------
  CalWireROI::CalWireROI(fhicl::ParameterSet const& pset):
    fROIPropertiesAlg(pset.get<fhicl::ParameterSet>("ROIPropertiesAlg"))
  {
    this->reconfigure(pset);

    produces< std::vector<recob::Wire> >(fSpillName);
    produces<art::Assns<raw::RawDigit, recob::Wire>>(fSpillName);
  }
  
  //-------------------------------------------------
  CalWireROI::~CalWireROI()
  {
  }

//////////////////////////////////////////////////////
void CalWireROI::reconfigure(fhicl::ParameterSet const& p)
{
    // Get signal shaping service.
    art::ServiceHandle<util::SignalShapingServiceMicroBooNE> sss;
    bool doInducedChargeDeconv = false;
    std::vector<std::vector<size_t> > respNums = sss->GetNResponses();
    for (size_t i = 0; i < respNums.at(1).size(); i++)
    {
        if (respNums.at(1).at(i) > 1) {
            doInducedChargeDeconv = true;
        }
    }

    // Throw exception if deconvolution should include dynamic induced charge effects (not yet implemented in CalROI) - M. Mooney
    if (doInducedChargeDeconv == true)
    {
        throw art::Exception(art::errors::Configuration)
            << "CalWireROI can not yet handle deconvolution with dynamic induced charge effects turned on.  Please use CalWireMicroBooNE instead.";
    }

    std::vector<unsigned short> uin;    std::vector<unsigned short> vin;
    std::vector<unsigned short> zin;
    
    fDigitModuleLabel     = p.get< std::string >                   ("DigitModuleLabel", "daq");
    fNoiseSource          = p.get< unsigned short >                ("NoiseSource",          3);
    fNumBinsHalf          = p.get< unsigned short >                ("NumBinsHalf",          3);
    fThreshold            = p.get< std::vector<unsigned short> >   ("Threshold"              );
    fNumSigma             = p.get< std::vector<int> >              ("NumSigma"               );
    uin                   = p.get< std::vector<unsigned short> >   ("uPlaneROIPad"           );
    vin                   = p.get< std::vector<unsigned short> >   ("vPlaneROIPad"           );
    zin                   = p.get< std::vector<unsigned short> >   ("zPlaneROIPad"           );
    fDoBaselineSub        = p.get< bool >                          ("DoBaselineSub"          );
    fuPlaneRamp           = p.get< bool >                          ("uPlaneRamp"             );
    fFFTSize              = p.get< int  >                          ("FFTSize"                );
    fSaveWireWF           = p.get< int >                           ("SaveWireWF"             );
    fMinAllowedChanStatus = p.get< int >                           ("MinAllowedChannelStatus");

    fDoBaselineSub_WaveformPropertiesAlg = p.get< bool >("DoBaselineSub_WaveformPropertiesAlg");
        
    if(uin.size() != 2 || vin.size() != 2 || zin.size() != 2) {
      throw art::Exception(art::errors::Configuration)
        << "u/v/z plane ROI pad size != 2";
    }

    fPreROIPad.resize(3);
    fPostROIPad.resize(3);
    
    // put the ROI pad sizes into more convenient vectors
    fPreROIPad[0]  = uin[0];
    fPostROIPad[0] = uin[1];
    fPreROIPad[1]  = vin[0];
    fPostROIPad[1] = vin[1];
    fPreROIPad[2]  = zin[0];
    fPostROIPad[2] = zin[1];
    
    fSpillName.clear();
    
    size_t pos = fDigitModuleLabel.find(":");
    if( pos!=std::string::npos ) {
      fSpillName = fDigitModuleLabel.substr( pos+1 );
      fDigitModuleLabel = fDigitModuleLabel.substr( 0, pos );
    }
    
    // re-initialize the FFT service for the request size
    // art::ServiceHandle<util::LArFFT> fFFT;
    // std::string options = fFFT->FFTOptions();
    // int fitbins = fFFT->FFTFitBins();
    // fFFT->ReinitializeFFT(fFFTSize, options, fitbins);
    //reconfFFT(fFFTSize);
    
}


void CalWireROI::reconfFFT(int temp_fftsize)
{
    // re-initialize the FFT service for the request size
    art::ServiceHandle<util::LArFFT> fFFT;

    if(fFFT->FFTSize() >= temp_fftsize) return;

    std::string options = fFFT->FFTOptions();
    int fitbins = fFFT->FFTFitBins();
    fFFT->ReinitializeFFT(temp_fftsize, options, fitbins);
}

//-------------------------------------------------
void CalWireROI::beginJob()
{
    fEventCount = 0;
} // beginJob

//////////////////////////////////////////////////////
void CalWireROI::endJob()
{
}
  
//////////////////////////////////////////////////////
void CalWireROI::produce(art::Event& evt)
{
    //get pedestal conditions
    const lariov::IDetPedestalProvider& pedestalRetrievalAlg = art::ServiceHandle<lariov::IDetPedestalService>()->GetPedestalProvider();
  
    // get the geometry
    art::ServiceHandle<geo::Geometry> geom;

    // get the FFT service to have access to the FFT size
    art::ServiceHandle<util::LArFFT> fFFT;
    reconfFFT(fFFTSize);
    
    // make a collection of Wires
    std::unique_ptr<std::vector<recob::Wire> > wirecol(new std::vector<recob::Wire>);
    // ... and an association set
    std::unique_ptr<art::Assns<raw::RawDigit,recob::Wire> > WireDigitAssn(new art::Assns<raw::RawDigit,recob::Wire>);

    // Read in the digit List object(s). 
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
    if(fSpillName.size()>0) evt.getByLabel(fDigitModuleLabel, fSpillName, digitVecHandle);
    else evt.getByLabel(fDigitModuleLabel, digitVecHandle);

    if (!digitVecHandle->size())  return;
    
    raw::ChannelID_t channel = raw::InvalidChannelID; // channel number
    
    const lariov::IChannelStatusProvider& chanFilt = art::ServiceHandle<lariov::IChannelStatusService>()->GetProvider();

    art::ServiceHandle<util::SignalShapingServiceMicroBooNE> sss;
    double deconNorm = sss->GetDeconNorm();

    // We'll need to set the transform size once we get the waveform and know its size
    size_t transformSize = 0;
    
    // loop over all wires
    wirecol->reserve(digitVecHandle->size());
    for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter)
    {
        // vector that will be moved into the Wire object
        recob::Wire::RegionsOfInterest_t ROIVec;
      
        // vector of ROI begin and end bins
        std::vector<std::pair<size_t, size_t>> rois;
      
        // get the reference to the current raw::RawDigit
        art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
        channel = digitVec->Channel();

        // The following test is meant to be temporary until the "correct" solution is implemented
        if (!chanFilt.IsPresent(channel)) continue;

        // Testing an idea about rejecting channels
        if (digitVec->GetPedestal() < 0.) continue;
        
        // skip bad channels
        if( chanFilt.Status(channel) >= fMinAllowedChanStatus)
        {
            size_t dataSize = digitVec->Samples();
            
            // vector holding uncompressed adc values
            std::vector<short> rawadc(dataSize);
            
            std::vector<geo::WireID> wids     = geom->ChannelToWire(channel);
            size_t                   thePlane = wids[0].Plane;
            
            // Set up the deconvolution and the vector to deconvolve
          
            // Set up the deconvolution and the vector to deconvolve
            // This is called only once per event, but under the hood nothing happens
            //   unless the FFT vector length changes (which it shouldn't for a run)
            if (!transformSize)
            {
                sss->SetDecon(dataSize, channel);
                transformSize = fFFT->FFTSize();
            }
            //sss->SetDecon(dataSize, channel);
            //size_t transformSize = fFFT->FFTSize();
            
            std::vector<float> rawAdcLessPedVec;
            
            rawAdcLessPedVec.resize(transformSize,0.);
            
            // uncompress the data
            raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
            
            // loop over all adc values and subtract the pedestal
            // When we have a pedestal database, can provide the digit timestamp as the third argument of GetPedestalMean
            size_t numBins   = 2 * fNumBinsHalf + 1;
            float  pdstl     = pedestalRetrievalAlg.PedMean(channel);
            float  rms_noise = digitVec->GetSigma();
            float  raw_noise = sss->GetRawNoise(channel);
            
            if      (fNoiseSource == 1) raw_noise = rms_noise;
            else if (fNoiseSource != 2) raw_noise = std::max(raw_noise,rms_noise);
            
            size_t binOffset(0);
            
            if (transformSize > dataSize) binOffset = (transformSize - dataSize) / 2;
            
            size_t startBin(binOffset);
            size_t stopBin(binOffset+numBins);
            
            float startThreshold = sqrt(float(numBins)) * (fNumSigma[thePlane] * raw_noise + fThreshold[thePlane]);
            float stopThreshold  = startThreshold;
            
            // Get the pedestal subtracted data, centered in the deconvolution vector
            std::transform(rawadc.begin(),rawadc.end(),rawAdcLessPedVec.begin()+startBin,[pdstl](const short& adc){return std::round(float(adc) - pdstl);});
            std::fill(rawAdcLessPedVec.begin(),rawAdcLessPedVec.begin()+startBin,0.);        //rawAdcLessPedVec.at(startBin));
            std::fill(rawAdcLessPedVec.begin()+startBin+dataSize,rawAdcLessPedVec.end(),0.); //rawAdcLessPedVec.at(startBin+dataSize-1));
            
            float runningSum = std::accumulate(rawAdcLessPedVec.begin()+startBin,rawAdcLessPedVec.begin()+stopBin, 0.);
            
            size_t roiStartBin(0);
            bool   roiCandStart(false);

            // search for ROIs - follow prescription from Bruce B using a running sum to make faster
            // Note that we start in the middle of the running sum... if we find an ROI padding will extend
            // past this to take care of ends of the waveform
            for(size_t bin = fNumBinsHalf + 1; bin < dataSize - fNumBinsHalf; bin++)
            {
                // handle the running sum
                // Case, we are at start of waveform
                runningSum -= rawAdcLessPedVec[startBin++];
                
                // Case, we are at end of waveform
                runningSum += rawAdcLessPedVec[stopBin++];
                
                // We have already started a candidate ROI
                if (roiCandStart)
                {
                    if (fabs(runningSum) < stopThreshold)
                    {
                        if (bin - roiStartBin > 2) rois.push_back(std::make_pair(roiStartBin, bin));
                        
                        roiCandStart = false;
                    }
                }
                // Not yet started a candidate ROI
                else
                {
                    if (fabs(runningSum) > startThreshold)
                    {
                        roiStartBin  = bin;
                        roiCandStart = true;
                    }
                }
            } // bin
            // add the last ROI if existed
            if (roiCandStart)
            {
                //unsigned int roiLen = dataSize -1 - roiStart;
                // if(roiLen > fMinWid)
                rois.push_back(std::make_pair(roiStartBin, dataSize-1));
            }

            // skip deconvolution if there are no ROIs
            if(rois.size() == 0) continue;

            // pad the ROIs
            for(size_t ii = 0; ii < rois.size(); ++ii)
            {
                // low ROI end
                rois[ii].first  = std::max(int(rois[ii].first - fPreROIPad[thePlane]),0);
                // high ROI end
                rois[ii].second = std::min(rois[ii].second + fPostROIPad[thePlane],dataSize - 1);
            }

            // merge the ROIs?
            if(rois.size() > 1)
            {
                // temporary vector for merged ROIs
                std::vector<std::pair<size_t,size_t>> trois;
          
                // Loop through candidate roi's
                size_t startRoi = rois.front().first;
                size_t stopRoi  = rois.front().second;
                
                for(size_t idx = 1; idx < rois.size(); idx++)
                {
                    if (rois[idx].first < stopRoi)
                    {
                        stopRoi = rois[idx].second;
                    }
                    else
                    {
                        trois.push_back(std::pair<size_t,size_t>(startRoi,stopRoi));
                        
                        startRoi = rois[idx].first;
                        stopRoi  = rois[idx].second;
                    }
                }
                
                // Make sure to get the last one
                trois.push_back(std::pair<size_t,size_t>(startRoi,stopRoi));
	  
                rois = trois;
            }
            
            // Strategy is to run deconvolution on the entire channel and then pick out the ROI's we found above
            sss->Deconvolute(channel,rawAdcLessPedVec);
            
            std::vector<float> holder;
            
            for(const auto roi : rois)
            {
                // First up: copy out the relevent ADC bins into the ROI holder
                size_t roiLen = roi.second - roi.first;
                
                holder.resize(roiLen);
                
                std::copy(rawAdcLessPedVec.begin()+binOffset+roi.first, rawAdcLessPedVec.begin()+binOffset+roi.second, holder.begin());
                std::transform(holder.begin(),holder.end(),holder.begin(),[deconNorm](float& deconVal){return deconVal/deconNorm;});

                // Now we do the baseline determination (and I'm left wondering if there is a better way using the entire waveform?)
                bool  baseSet(false);
                float base(0.);
                if(fDoBaselineSub && fPreROIPad[thePlane] > 0 )
                {
                    //1. Check Baseline match?
                    // If not, include next ROI(if none, go to the end of signal)
                    // If yes, proceed
                    size_t binsToAve(20);
                    float  basePre  = std::accumulate(holder.begin(),holder.begin()+binsToAve,0.) / float(binsToAve);
                    float  basePost = std::accumulate(holder.end()-binsToAve,holder.end(),0.) / float(binsToAve);
                   
                    baseSet = true;
                    base    = SubtractBaseline(holder, basePre,basePost,roi.first,roiLen,dataSize);
                } // fDoBaselineSub ...
                else if(fDoBaselineSub_WaveformPropertiesAlg)
                {
                    baseSet = true;
                    base    = fROIPropertiesAlg.GetWaveformPedestal(holder);
                }
                
                if (baseSet) std::transform(holder.begin(),holder.end(),holder.begin(),[base](float& adcVal){return adcVal - base;});

                // add the range into ROIVec
                ROIVec.add_range(roi.first, std::move(holder));
            }
        } // end if not a bad channel

        // create the new wire directly in wirecol
        wirecol->push_back(recob::WireCreator(std::move(ROIVec),*digitVec).move());
        
        // add an association between the last object in wirecol
        // (that we just inserted) and digitVec
        if (!util::CreateAssn(*this, evt, *wirecol, digitVec, *WireDigitAssn, fSpillName)) {
            throw art::Exception(art::errors::InsertFailure)
                << "Can't associate wire #" << (wirecol->size() - 1)
                << " with raw digit #" << digitVec.key();
        } // if failed to add association
        //  DumpWire(wirecol->back()); // for debugging
    }

    if(wirecol->size() == 0)
      mf::LogWarning("CalWireROI") << "No wires made for this event.";

    //Make Histogram of recob::wire objects from Signal() vector
    // get access to the TFile service
    if ( fSaveWireWF ){
        art::ServiceHandle<art::TFileService> tfs;
        for (size_t wireN = 0; wireN < wirecol->size(); wireN++){
            std::vector<float> sigTMP = wirecol->at(wireN).Signal();
            TH1D* fWire = tfs->make<TH1D>(Form("Noise_Evt%04zu_N%04zu",fEventCount,wireN), ";Noise (ADC);",
				      sigTMP.size(),-0.5,sigTMP.size()-0.5);
            for (size_t tick = 0; tick < sigTMP.size(); tick++){
                fWire->SetBinContent(tick+1, sigTMP.at(tick) );
            }
        }
    }
    
    evt.put(std::move(wirecol), fSpillName);
    evt.put(std::move(WireDigitAssn), fSpillName);

    fEventCount++;

    return;
} // produce


float CalWireROI::SubtractBaseline(std::vector<float>& holder,
                                   float               basePre,
                                   float               basePost,
                                   size_t              roiStart,
                                   size_t              roiLen,
                                   size_t              dataSize)
{
    float base=0;

    //can not trust the early part
    if (roiStart < 20 && roiStart + roiLen < dataSize - 20){
        base = basePost;
    // can not trust the later part
    }else if (roiStart >= 20 && roiStart + roiLen >= dataSize - 20){
        base = basePre;
    // can trust both
    }else if (roiStart >= 20 && roiStart + roiLen < dataSize - 20){
        if (fabs(basePre-basePost)<3){
            base = (basePre+basePost)/2.;
        }else{
            if (basePre < basePost){
                base = basePre;
            }else{
                base = basePost;
            }
        }
    // can not use both
    }else{
        float min = 0,max=0;
        for (unsigned int bin = 0; bin < roiLen; bin++){
            if (holder[bin] > max) max = holder[bin];
            if (holder[bin] < min) min = holder[bin];
        }
        int nbin = max - min;
        if (nbin!=0){
            TH1F *h1 = new TH1F("h1","h1",nbin,min,max);
            for (unsigned int bin = 0; bin < roiLen; bin++){
                h1->Fill(holder[bin]);
            }
            float ped = h1->GetMaximum();
            float ave=0,ncount = 0;
            for (unsigned int bin = 0; bin < roiLen; bin++){
                if (fabs(holder[bin]-ped)<2){
                    ave +=holder[bin];
                    ncount ++;
                }
            }
            if (ncount==0) ncount=1;
            ave = ave/ncount;
            h1->Delete();
            base = ave;
        }
    }
    
    return base;
}


void CalWireROI::doDecon(std::vector<float>&                                       holder,
                         raw::ChannelID_t                                          channel,
                         unsigned int                                              thePlane,
                         const std::vector<std::pair<size_t,size_t>>&              rois,
                         const std::vector<std::pair<size_t,size_t>>&              holderInfo,
                         recob::Wire::RegionsOfInterest_t&                         ROIVec,
                         art::ServiceHandle<util::SignalShapingServiceMicroBooNE>& sss)
{
    sss->Deconvolute(channel,holder);

    // transfer the ROIs and start bins into the vector that will be
    // put into the event
    for(size_t jr = 0; jr < holderInfo.size(); ++jr)
    {
        std::vector<float> sigTemp;
        size_t bBegin = holderInfo[jr].first;
        size_t theROI = holderInfo[jr].second;
        size_t roiLen = rois[theROI].second - rois[theROI].first;
        size_t bEnd = bBegin + roiLen;
        float basePre = 0., basePost = 0.;
        // Baseline subtraction if requested and the ROIs are padded.
        // Can't baseline subtract signals when the ROI start is close to 0 either
        if(fDoBaselineSub && fPreROIPad[thePlane] > 0 &&
           rois[theROI].first > fPreROIPad[thePlane])
        {
            // find the baseline from the first few bins in the leading Pad region
            unsigned short bbins = fPreROIPad[thePlane];
            size_t bin;
            if(bbins > 5) bbins = 5;
            for(bin = 0; bin < bbins; ++bin) {
                basePre  += holder[bBegin + bin];
                basePost += holder[bEnd - bin];
            }
            basePre /= (float)bbins;
            basePost /= (float)bbins;
            float slp = (basePost - basePre) / (float)(roiLen - bbins);
            float base;
            for(size_t jj = bBegin; jj < bEnd; ++jj) {
                base = basePre + slp * (jj - bBegin);
                sigTemp.push_back(holder[jj] - base);
            } // jj
        } // fDoBaselineSub ...
        else {
            for(size_t jj = bBegin; jj < bEnd; ++jj) sigTemp.push_back(holder[jj]);
        } // !fDoBaselineSub ...
      
        // add the range into ROIVec
        ROIVec.add_range(rois[theROI].first, std::move(sigTemp));
    } // jr
} // doDecon


} // end namespace caldata
