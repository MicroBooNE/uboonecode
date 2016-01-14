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

/* unused function
namespace {

	void DumpWire(const recob::Wire& wire) {
		
		const size_t pagesize = 10;
		mf::LogDebug log("DumpWire");
		
		const recob::Wire::RegionsOfInterest_t& wireSignalRoI = wire.SignalROI();
		
		log << "\nDumpWire: wire on view " << ((int) wire.View())
			<< " channel " << wire.Channel()
			<< " with " << wireSignalRoI.n_ranges() << " regions of interest:";
		size_t iRoI = 0;
		auto RoI = wireSignalRoI.begin_range(), rend = wireSignalRoI.end_range();
		while (RoI != rend) {
			++iRoI;
			log << "\nDumpWire: [RoI " << iRoI << "] starts at " << RoI->begin_index()
				<< ", " << RoI->size() << " samples:";
			size_t iSample = 0;
			for (auto sample: *RoI) {
				if (iSample % pagesize == 0)
					log << "\nDumpWire: [RoI " << iRoI << "/" << iSample << "]";
				log << '\t' << sample;
				++iSample;
			} // for sample
			++RoI;
		} // for RoI
		
		
		const recob::Wire::RegionsOfInterest_t& wireSignal = wireSignalRoI;
		std::vector<float> buffer(pagesize), prev_buffer;
		log << "\nDumpWire: wire on view " << ((int)wire.View())
			<< " channel " << wire.Channel()
			<< " with " << wireSignal.size() << " samples:";
		size_t i = 0, nSame = 0;
		auto iSample = wireSignal.begin(), send = wireSignal.end();
		while (iSample < send) {
			i += prev_buffer.size();
			buffer.assign(iSample, std::min(iSample + pagesize, send));
			iSample += buffer.size();
			if (buffer == prev_buffer) {
				++nSame;
			}
			else {
				if (nSame > 0) {
					log << "\nDumpWire: [" << i << "]  ... and " << nSame << " more";
					nSame = 0;
				}
				
				if (i % pagesize == 0) log << "\nDumpWire: [" << i << "]";
				for (auto value: buffer) log << '\t' << value;
				buffer.swap(prev_buffer);
			}
		} // while
		if (nSame > 0) {
			log << "\nDumpWire: [" << i << "]  ... and " << nSame << " more to the end";
			nSame = 0;
		}
	} // DumpWire()
}
*/

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
    unsigned short              fMinWid;               ///< min width for a ROI
    unsigned short              fMinSep;               ///< min separation between ROIs
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
    fMinWid               = p.get< unsigned short >                ("MinWid"                 );
    fMinSep               = p.get< unsigned short >                ("MinSep"                 );
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
    unsigned int bin(0);     // time bin loop variable
    
    const lariov::IChannelStatusProvider& chanFilt = art::ServiceHandle<lariov::IChannelStatusService>()->GetProvider();

    art::ServiceHandle<util::SignalShapingServiceMicroBooNE> sss;
    double DeconNorm = sss->GetDeconNorm();
    
    // loop over all wires
    wirecol->reserve(digitVecHandle->size());
    for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter)
    {
        // vector that will be moved into the Wire object
        recob::Wire::RegionsOfInterest_t ROIVec;
      
        // the starting position and length of each ROI in the packed holder vector
        std::vector<std::pair<size_t, size_t>> holderInfo;
        // vector of ROI begin and end bins
        std::vector<std::pair<size_t, size_t>> rois;
      
        // get the reference to the current raw::RawDigit
        art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
        channel = digitVec->Channel();

        // The following test is meant to be temporary until the "correct" solution is implemented
        if (!chanFilt.IsPresent(channel)) continue;

        // Testing an idea about rejecting channels
        if (digitVec->GetPedestal() < 0.) continue;

        size_t dataSize = digitVec->Samples();
        // vector holding uncompressed adc values
        std::vector<short> rawadc(dataSize);
      
        std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
        size_t thePlane = wids[0].Plane;
        //unsigned int theWire = wids[0].Wire;
      
        // use short pedestal for testing the 1st induction plane
        // unsigned short sPed = abs(fThreshold[thePlane]);
    
        // Find minimum separation between ROIs including the padding.
        // Use the longer pad
        // unsigned int minSepPad = fMinSep + fPreROIPad[thePlane];
        //if(fPostROIPad[thePlane] > fPreROIPad[thePlane])
        //  minSepPad = fMinSep + fPostROIPad[thePlane];
      
        // skip bad channels
        if( chanFilt.Status(channel) >= fMinAllowedChanStatus)
        {
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
            
            size_t startBin(0);
            size_t stopBin(numBins);
            
            float startThreshold = sqrt(float(numBins))*6.*raw_noise;
            float stopThreshold  = 0.5 * startThreshold;
            
            std::vector<float> rawAdcLessPedVec(rawadc.size());
            
            std::transform(rawadc.begin(),rawadc.end(),rawAdcLessPedVec.begin(),[pdstl](const short& adc){return float(adc) - pdstl;});
            
            float runningSum = std::accumulate(rawAdcLessPedVec.begin(),rawAdcLessPedVec.begin() + stopBin, 0.); //, [](float sum, float add){return sum + fabs(add);});
            
            size_t roiStartBin(0);
            bool   roiCandStart(false);

            // search for ROIs - follow prescription from Bruce B using a running sum to make faster
            // Note that we start in the middle of the running sum... if we find an ROI padding will extend
            // past this to take care of ends of the waveform
            for(bin = fNumBinsHalf; bin < dataSize - fNumBinsHalf; ++bin)
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
                        rois.push_back(std::make_pair(roiStartBin, bin));
                        
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
            
            holderInfo.clear();

            for (unsigned int ir = 0; ir < rois.size(); ++ir)
            {
                size_t roiLen = rois[ir].second - rois[ir].first + 1;
                size_t roiStart = rois[ir].first;
	  
                int flag =1;
                float tempPre=0,tempPost=0;
                std::vector<float> holder;
                while(flag)
                {
                    unsigned int transformSize = fFFTSize; //current transformsize
                    //if ROI length is longer, take ROI length
                    if (roiLen > transformSize) transformSize = roiLen;
	    
                    // Get signal shaping service.
                    sss->SetDecon(transformSize);
                    transformSize = fFFT->FFTSize();
                    
                    // temporary vector of signals
                    holder.resize(transformSize,0.);
	    
                    // Apparently the transformSize gets reset in the call to SetDecon so make sure we are not running off the
                    // end of the waveform here...
                    size_t copySize(holder.size());
                    size_t leftoverSize(0);
                    
                    if (roiStart + copySize > dataSize)
                    {
                        leftoverSize  = roiStart + copySize - dataSize;
                        copySize     -= leftoverSize;
                        flag          = 0;
                    }
                    
                    std::copy(rawAdcLessPedVec.begin() + roiStart, rawAdcLessPedVec.begin() + roiStart + copySize, holder.begin());
                    
                    if (leftoverSize > 0)
                        std::copy(rawAdcLessPedVec.begin(), rawAdcLessPedVec.begin() + leftoverSize, holder.begin() + copySize);

                    sss->Deconvolute(channel,holder);
                    for(bin = 0; bin < holder.size(); ++bin) holder[bin]=holder[bin]/DeconNorm;

                    //1. Check Baseline match?
                    // If not, include next ROI(if none, go to the end of signal)
                    // If yes, proceed
                    tempPre=0,tempPost=0;
                    for(unsigned int bin = 0; bin < 20; ++bin)
                    {
                        tempPre  += holder[bin];
                        tempPost += holder[roiLen -1 - bin];
                    }
                    tempPre = tempPre/20.;
                    tempPost = tempPost/20.;
	    
                    double deconNoise = sss->GetDeconNoise(channel)/sqrt(10.)*4;
	    
                    if (fabs(tempPost-tempPre)<deconNoise){
                        flag = 0;
                    }else{
                        if (tempPre > tempPost && roiStart <= 2){
                            //if (tempPre > tempPost){
                            flag = 0;
                        }else{
                            if (ir+1<rois.size()){
                                roiLen += 100;
                                if (roiLen >= rois[ir+1].first - roiStart + 1)
                                    roiLen = rois[++ir].second - roiStart + 1;
                            }else{
                                roiLen += 100;
                                if (roiLen>dataSize-roiStart)
                                    roiLen = dataSize - roiStart;
                            }
                        }
                    }
                }
	  
                // transfer the ROI and start bins into the vector that will be
                // put into the event
                std::vector<float> sigTemp;
                unsigned int bBegin = 0;
                //unsigned int theROI =ir;
                unsigned int bEnd = bBegin + roiLen;
                float basePre = 0., basePost = 0.;

                float base=0;
                if(fDoBaselineSub && fPreROIPad[thePlane] > 0 )
                {
                    basePre =tempPre;
                    basePost=tempPost;
                    base = SubtractBaseline(holder, basePre,basePost,roiStart,roiLen,dataSize);
                } // fDoBaselineSub ...
                else if(fDoBaselineSub_WaveformPropertiesAlg)
                {
                    holder.resize(roiLen);
                    base = fROIPropertiesAlg.GetWaveformPedestal(holder);
                }

                for(unsigned int jj = bBegin; jj < bEnd; ++jj)
                    sigTemp.push_back(holder[jj]-base);
	        
                // add the range into ROIVec
                ROIVec.add_range(roiStart, std::move(sigTemp));
	  
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
