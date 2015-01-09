////////////////////////////////////////////////////////////////////////
//
// CalWireROI class - variant of CalWire that deconvolves in 
// Regions Of Interest
//
// baller@fnal.gov
//
//
////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <iomanip>

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
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "art/Utilities/Exception.h"

#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"
#include "Geometry/Geometry.h"
#include "Filters/ChannelFilter.h"
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "RecoBase/Wire.h"
#include "Utilities/LArFFT.h"


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

  class CalWireROI : public art::EDProducer {

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
    
    std::string  fDigitModuleLabel;  ///< module that made digits

    std::string  fSpillName;  ///< nominal spill is an empty string
                              ///< it is set by the DigitModuleLabel
                              ///< ex.:  "daq:preSpill" for prespill data

    std::vector<unsigned short> fThreshold;   ///< abs(threshold) ADC counts for ROI
    unsigned short fMinWid;         ///< min width for a ROI
    unsigned short fMinSep;          ///< min separation between ROIs
    int fFFTSize;        ///< FFT size for ROI deconvolution
    std::vector<unsigned short> fPreROIPad; ///< ROI padding
    std::vector<unsigned short> fPostROIPad; ///< ROI padding
    bool fDoBaselineSub;  ///< Do baseline subtraction after deconvolution?
    bool fuPlaneRamp;     ///< set true for correct U plane wire response
    int  fSaveWireWF;     ///< Save recob::wire object waveforms
    size_t fEventCount;  ///< count of event processed

    void doDecon(std::vector<float>& holder, 
      uint32_t channel, unsigned int thePlane,
      std::vector<std::pair<unsigned int, unsigned int>> rois,
      std::vector<std::pair<unsigned int, unsigned int>> holderInfo,
      recob::Wire::RegionsOfInterest_t& ROIVec,
      art::ServiceHandle<util::SignalShapingServiceMicroBooNE>& sss);
    
  protected: 
    
  }; // class CalWireROI

  DEFINE_ART_MODULE(CalWireROI)
  
  //-------------------------------------------------
  CalWireROI::CalWireROI(fhicl::ParameterSet const& pset)
  {
    fSpillName="";
    this->reconfigure(pset);

    if(fSpillName.size()<1) produces< std::vector<recob::Wire> >();
    else produces< std::vector<recob::Wire> >(fSpillName);
  }
  
  //-------------------------------------------------
  CalWireROI::~CalWireROI()
  {
  }

  //////////////////////////////////////////////////////
  void CalWireROI::reconfigure(fhicl::ParameterSet const& p)
  {
  
    std::vector<unsigned short> uin;    std::vector<unsigned short> vin;
    std::vector<unsigned short> zin;
    
    fDigitModuleLabel = p.get< std::string >          ("DigitModuleLabel", "daq");
    fThreshold        = p.get< std::vector<unsigned short> >   ("Threshold");
    fMinWid           = p.get< unsigned short >       ("MinWid");
    fMinSep           = p.get< unsigned short >       ("MinSep");
    uin               = p.get< std::vector<unsigned short> >   ("uPlaneROIPad");
    vin               = p.get< std::vector<unsigned short> >   ("vPlaneROIPad");
    zin               = p.get< std::vector<unsigned short> >   ("zPlaneROIPad");
    fDoBaselineSub    = p.get< bool >                 ("DoBaselineSub");
    fuPlaneRamp       = p.get< bool >                 ("uPlaneRamp");
    fFFTSize          = p.get< int  >                 ("FFTSize");
    fSaveWireWF       = p.get< int >                  ("SaveWireWF");
    
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
    
    fSpillName="";
    
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
    reconfFFT(fFFTSize);
    
  }


  void CalWireROI::reconfFFT(int temp_fftsize){
    // re-initialize the FFT service for the request size
    art::ServiceHandle<util::LArFFT> fFFT;
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
    // get the geometry
    art::ServiceHandle<geo::Geometry> geom;

    // get the FFT service to have access to the FFT size
    art::ServiceHandle<util::LArFFT> fFFT;
  
    
    // make a collection of Wires
    std::unique_ptr<std::vector<recob::Wire> > wirecol(new std::vector<recob::Wire>);

   
    
    // Read in the digit List object(s). 
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
    if(fSpillName.size()>0) evt.getByLabel(fDigitModuleLabel, fSpillName, digitVecHandle);
    else evt.getByLabel(fDigitModuleLabel, digitVecHandle);

    if (!digitVecHandle->size())  return;
    
    uint32_t     channel(0); // channel number
    unsigned int bin(0);     // time bin loop variable
    
    filter::ChannelFilter *chanFilt = new filter::ChannelFilter();  
    
    // loop over all wires
    wirecol->reserve(digitVecHandle->size());
    for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter){
      
      // vector that will be moved into the Wire object
      recob::Wire::RegionsOfInterest_t ROIVec;
      
      // the starting position and length of each ROI in the packed holder vector
      std::vector<std::pair<unsigned int, unsigned int>> holderInfo;
      // vector of ROI begin and end bins
      std::vector<std::pair<unsigned int, unsigned int>> rois;
      
      // get the reference to the current raw::RawDigit
      art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
      channel = digitVec->Channel();
      unsigned int dataSize = digitVec->Samples();
      // vector holding uncompressed adc values
      std::vector<short> rawadc(dataSize);
      
      std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
      unsigned int thePlane = wids[0].Plane;
      //unsigned int theWire = wids[0].Wire;
      
      // use short pedestal for testing the 1st induction plane
      unsigned short sPed = abs(fThreshold[thePlane]);
    
      // Find minimum separation between ROIs including the padding.
      // Use the longer pad
      // unsigned int minSepPad = fMinSep + fPreROIPad[thePlane];
      //if(fPostROIPad[thePlane] > fPreROIPad[thePlane]) 
      //  minSepPad = fMinSep + fPostROIPad[thePlane];
      
      // skip bad channels
      if(!chanFilt->BadChannel(channel)) {
        
        // uncompress the data
        raw::Uncompress(digitVec->fADC, rawadc, digitVec->Compression());
        // loop over all adc values and subtract the pedestal
        float pdstl = digitVec->GetPedestal();
	//subtract time-offset added in SimWireMicroBooNE_module
	// Xin, remove the time offset
	//int time_offset = 0.;//sss->FieldResponseTOffset(channel);
        unsigned int roiStart = 0;
       
        // search for ROIs
        for(bin = 1; bin < dataSize; ++bin) {
          float SigVal = fabs(rawadc[bin] - pdstl);
          if(roiStart == 0) {
            // not in a ROI
            // Handle the onset of a ROI differently for the 1st induction plane
            // if it has been modeled correctly
            if(fuPlaneRamp && thePlane == 0) {
              if(rawadc[bin] - rawadc[bin - 1] < sPed) roiStart = bin;
            } else {
              if(SigVal > fThreshold[thePlane]) roiStart = bin;
            }
	  } else {
            // leaving a ROI?
            if(SigVal < fThreshold[thePlane]) {
              // is the ROI wide enough?
              unsigned int roiLen = bin - roiStart;
              // if(roiLen > transformSize) {
              //   mf::LogError("CalWireROI")<<"ROI too long "
              //     <<roiLen<<" on plane:wire "<<thePlane<<":"<<theWire;
              //   break;
              // }
              //if(roiLen > fMinWid && roiLen < transformSize) 
	      if(roiLen > fMinWid) 
                rois.push_back(std::make_pair(roiStart, bin));
              roiStart = 0;
            }
          } // roiStart test
        } // bin
	// add the last ROI if existed
	if (roiStart!=0){
	   unsigned int roiLen = dataSize -1 - roiStart;
	   if(roiLen > fMinWid) 
	     rois.push_back(std::make_pair(roiStart, dataSize-1));
	   roiStart = 0;
	}

        // skip deconvolution if there are no ROIs
        if(rois.size() == 0) continue;

        holderInfo.clear();
	//        for(unsigned int bin = 0; bin < holder.size(); ++bin) holder[bin] = 0;
        
	 // pad the ROIs
        for(unsigned int ii = 0; ii < rois.size(); ++ii) {
          // low ROI end
          int low = rois[ii].first - fPreROIPad[thePlane];
          if(low < 0) low = 0;
          rois[ii].first = low;
          // high ROI end
          unsigned int high = rois[ii].second + fPostROIPad[thePlane];
          if(high >= dataSize) high = dataSize-1;
          rois[ii].second = high;
        }

	// if (channel==3218){
	//   std::cout << "Xin " << " " << channel << " " << rois.size() << std::endl;
	//   for(unsigned int ii = 0; ii < rois.size(); ++ii) {
	//     std::cout << rois[ii].first << " " << rois[ii].second << std::endl;
	//   }
	// }

        // merge the ROIs?
        if(rois.size() > 1) {
          // temporary vector for merged ROIs
          std::vector<std::pair<unsigned int, unsigned int>> trois;
          
	  for (unsigned int ii = 0; ii<rois.size();ii++){
	    unsigned int roiStart = rois[ii].first;
	    unsigned int roiEnd = rois[ii].second;
          
	    int flag1 = 1;
	    unsigned int jj=ii+1;
	    while(flag1){	
	      if (jj<rois.size()){
		if(rois[jj].first - roiEnd >=0 ) {
		  roiEnd = rois[jj].second;
		  ii = jj;
		  jj = ii+1;
		}else{
		  flag1 = 0;
		}
	      }else{
		flag1 = 0;
	      }
	    }
	    

	    trois.push_back(std::make_pair(roiStart,roiEnd));	    
	  }
	  
	  rois = trois;
	}
	  
	

	for (unsigned int ir = 0; ir < rois.size(); ++ir) {
	  unsigned int roiLen = rois[ir].second - rois[ir].first;
	  unsigned int roiStart = rois[ir].first;
	  //treat FFT Size

	  int flag =1;
	  float tempPre=0,tempPost=0;
	  std::vector<float> holder;
	  while(flag){
	    
	    unsigned int transformSize = fFFTSize; //current transformsize
	    //if ROI length is longer, take ROI length
	    if (roiLen > transformSize) transformSize = roiLen;
	    
	    // if (transformSize < 8192) transformSize=8192;

	    // Get signal shaping service.
	    art::ServiceHandle<util::SignalShapingServiceMicroBooNE> sss;
	    sss->SetDecon(transformSize);
	    transformSize = fFFT->FFTSize();
	    // temporary vector of signals
	    holder.resize(transformSize,0);
	    
	    unsigned int hBin = 0;
	    for(unsigned int bin = roiStart; bin < roiStart + holder.size(); ++bin) {
	      if (bin < dataSize){
		holder[hBin] = rawadc[bin]-pdstl;
	      }else{
		holder[hBin] = rawadc[bin-dataSize]-pdstl;
		flag = 0;
	      }
	      ++hBin;
	    } // bin

	    

	    sss->Deconvolute(channel,holder);
	    if (channel==3218){
	      for(unsigned int bin = 0; bin <holder.size(); ++bin) {
		std::cout << bin << " " <<  holder[bin] << std::endl;
	      }
	    }
	    //1. Check Baseline match?
	    // If not, include next ROI(if none, go to the end of signal)
	    // If yes, proceed
	    tempPre=0,tempPost=0;
	    for(unsigned int bin = 0; bin < 5; ++bin) {
	      tempPre  += holder[bin];
	      tempPost += holder[roiLen - bin];
	    }
	    tempPre = tempPre/5.;
	    tempPost = tempPost/5.;

	    // if (channel==3218){
	    //   std::cout << flag << " " << tempPre << " " << tempPost << std::endl;
	    // }

	    if (fabs(tempPost-tempPre)<2){
	      flag = 0;
	    }else{
	      ir++;
	      if (ir<rois.size()){
		roiLen = rois[ir].second - roiStart;
	      }else{
		roiLen = dataSize -1 - roiStart;
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

	  if(fDoBaselineSub && fPreROIPad[thePlane] > 0 ) {
	    basePre =tempPre;
	    basePost=tempPost;
	    float base;
	    if (fabs(basePre-basePost)<3){
	      base = (basePre+basePost)/2.;
	    }else{
	      base = basePre;
	    }
	    for(unsigned int jj = bBegin; jj < bEnd; ++jj) {
	      sigTemp.push_back(holder[jj] - base);
	    } // jj
	  } // fDoBaselineSub ...
	  else {
	    for(unsigned int jj = bBegin; jj < bEnd; ++jj) sigTemp.push_back(holder[jj]);
	  } // !fDoBaselineSub ...
	        
	  // add the range into ROIVec 
	  ROIVec.add_range(roiStart, std::move(sigTemp));
	  
	}
	

       // unsigned int hBin = 0;
       //  // max number of holder bins that should be packed. This is the
       //  // FFT size - width of the response function. Do this crudely for
       //  // now since the response functions are < 20 time bins in length
       //  unsigned int packLen = transformSize - 20;
       //  for(unsigned int ir = 0; ir < rois.size(); ++ir) {
       //    // enough room to continue packing?
       //    unsigned int roiLen = rois[ir].second - rois[ir].first;
       //    if(roiLen > transformSize) {
       //      mf::LogWarning("CalWireROI")<<"ROI length "<<roiLen
       //        <<" > FFT size "<<transformSize<<" on plane:wire "
       //        <<thePlane<<":"<<theWire;
       //      continue;
       //    }
       //    // not enough room to continue packing. Deconvolute this batch, stuff it
       //    // into ROIVec and continue
       //    if(hBin + roiLen > packLen) {
       //      // need to de-convolute holder and continue packing
       //      doDecon(holder, channel, thePlane, rois, holderInfo, ROIVec, sss);
       //      holderInfo.clear();
       //      hBin = 0;
       //      for(unsigned int bin = 0; bin < holder.size(); ++bin) holder[bin] = 0;
       //    } // hBin + roiLen > packLen
       //    // save the position of this ROI in the holder vector and the ROI index
       //    holderInfo.push_back(std::make_pair(hBin, ir));
       //    for(unsigned int bin = rois[ir].first; bin < rois[ir].second; ++bin) {
       // 	    if ( (hBin-time_offset >= 0) and (hBin-time_offset < holder.size()) ){
       // 	      holder[hBin-time_offset] = rawadc[bin]-pdstl;
       // 	      ++hBin;
       // 	    }//if time-offset does not force to go out of holder bounds
       //    } // bin
       //  } // ir < rois.size
       //  // do the last deconvolution if needed
       //  if(holderInfo.size() > 0)
       //    doDecon(holder, channel, thePlane, rois, holderInfo, ROIVec, sss);
      



      } // end if not a bad channel 

      // create the new wire directly in wirecol
      wirecol->emplace_back(std::move(ROIVec), digitVec);
      
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
    

    if(fSpillName.size()>0)
      evt.put(std::move(wirecol), fSpillName);
    else evt.put(std::move(wirecol));

    delete chanFilt;

    fEventCount++;

    return;
  } // produce

  void CalWireROI::doDecon(std::vector<float>& holder, 
    uint32_t channel, unsigned int thePlane,
    std::vector<std::pair<unsigned int, unsigned int> > rois,
    std::vector<std::pair<unsigned int, unsigned int> > holderInfo,
    recob::Wire::RegionsOfInterest_t& ROIVec,
    art::ServiceHandle<util::SignalShapingServiceMicroBooNE>& sss)
  {
    

    sss->Deconvolute(channel,holder);

   
    // transfer the ROIs and start bins into the vector that will be
    // put into the event
    for(unsigned int jr = 0; jr < holderInfo.size(); ++jr) {
      std::vector<float> sigTemp;
      unsigned int bBegin = holderInfo[jr].first;
      unsigned int theROI = holderInfo[jr].second;
      unsigned int roiLen = rois[theROI].second - rois[theROI].first;
      unsigned int bEnd = bBegin + roiLen;
      float basePre = 0., basePost = 0.;
      // Baseline subtraction if requested and the ROIs are padded.
      // Can't baseline subtract signals when the ROI start is close to 0 either
      if(fDoBaselineSub && fPreROIPad[thePlane] > 0 && 
          rois[theROI].first > fPreROIPad[thePlane]) {
        // find the baseline from the first few bins in the leading Pad region
        unsigned short bbins = fPreROIPad[thePlane];
        unsigned int bin;
        if(bbins > 5) bbins = 5;
        for(bin = 0; bin < bbins; ++bin) {
          basePre  += holder[bBegin + bin];
          basePost += holder[bEnd - bin];
        }
        basePre /= (float)bbins;
        basePost /= (float)bbins;
        float slp = (basePost - basePre) / (float)(roiLen - bbins);
        float base;
        for(unsigned int jj = bBegin; jj < bEnd; ++jj) {
          base = basePre + slp * (jj - bBegin);
          sigTemp.push_back(holder[jj] - base);
        } // jj
      } // fDoBaselineSub ...
      else {
        for(unsigned int jj = bBegin; jj < bEnd; ++jj) sigTemp.push_back(holder[jj]);
      } // !fDoBaselineSub ...
      
      // add the range into ROIVec 
      ROIVec.add_range(rois[theROI].first, std::move(sigTemp));
    } // jr
  } // doDecon


} // end namespace caldata
