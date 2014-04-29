////////////////////////////////////////////////////////////////////////
//
// CalWireROI class - variant of CalWire that deconvolues in 
// Regions Of Interest
//
// baller@fnal.gov
//
//
////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>

#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Utilities/Exception.h"

#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"
#include "Geometry/Geometry.h"
#include "Filters/ChannelFilter.h"
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "RecoBase/Wire.h"
#include "Utilities/LArFFT.h"


namespace {

/* // old recob::Wire version
	void DumpWire(const recob::Wire& wire) {
		
		const size_t pagesize = 10;
		mf::LogDebug log("DumpWire");
		
		const recob::Wire::RegionsOfInterest_t& wireSignalRoI = wire.SignalROI();
		
		log << "\nDumpWire: wire on view " << ((int) wire.View())
			<< " channel " << wire.Channel()
			<< " with " << wireSignalRoI.size() << " regions of interest:";
		size_t iRoI = 0;
		for (auto RoI: wireSignalRoI) {
			log << "\nDumpWire: [RoI " << iRoI << "] starts at " << RoI.first
				<< ", " << RoI.second.size() << " samples:";
			++iRoI;
			for (size_t i = 0; i < RoI.second.size(); ++i) {
				if (i % pagesize == 0)
					log << "\nDumpWire: [RoI " << iRoI << "/" << i << "]";
				log << '\t' << RoI.second[i];
			} // for sample
		} // for RoI
		
		
		const std::vector<float> wireSignal{wire.Signal()};
		std::vector<float> buffer(pagesize), prev_buffer;
		log << "\nDumpWire: wire on view " << ((int)wire.View())
			<< " channel " << wire.Channel()
			<< " with " << wireSignal.size() << " samples:";
		size_t i = 0, nSame = 0;
		auto iSample = wireSignal.begin(), send = wireSignal.end();
		while (iSample != send) {
			i += prev_buffer.size();
			buffer.assign(iSample, iSample + pagesize);
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
*/
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
 
  private:
    
    std::string  fDigitModuleLabel;  ///< module that made digits

    std::string  fSpillName;  ///< nominal spill is an empty string
                              ///< it is set by the DigitModuleLabel
                              ///< ex.:  "daq:preSpill" for prespill data

    std::vector<unsigned short> fThreshold;   ///< abs(threshold) ADC counts for ROI
    unsigned short fMinWid;         ///< min width for a ROI
    unsigned short fMinSep;          ///< min separation between ROIs
    std::vector<unsigned short> fPreROIPad; ///< ROI padding
    std::vector<unsigned short> fPostROIPad; ///< ROI padding
    bool fDoBaselineSub;  ///< Do baseline subtraction after deconvolution?
    bool fuPlaneRamp;     ///< set true for correct U plane wire response

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
  
    std::vector<unsigned short> uin;
    std::vector<unsigned short> vin;
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
    
  }

  //-------------------------------------------------
  void CalWireROI::beginJob()
  {
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
    unsigned int transformSize = fFFT->FFTSize();
    
    // make a collection of Wires
    std::unique_ptr<std::vector<recob::Wire> > wirecol(new std::vector<recob::Wire>);

    // Get signal shaping service.
    art::ServiceHandle<util::SignalShapingServiceMicroBooNE> sss;
    
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

      // temporary vector of signals
      std::vector<float> holder(transformSize, 0.);
      
      // get the reference to the current raw::RawDigit
      art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
      channel = digitVec->Channel();
      unsigned int dataSize = digitVec->Samples();
      // vector holding uncompressed adc values
      std::vector<short> rawadc(dataSize);
      
      std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
      unsigned int thePlane = wids[0].Plane;
      unsigned int theWire = wids[0].Wire;
      
      // use short pedestal for testing the 1st induction plane
      unsigned short sPed = abs(fThreshold[thePlane]);
    
      // Find minimum separation between ROIs including the padding.
      // Use the longer pad
      unsigned int minSepPad = fMinSep + fPreROIPad[thePlane];
      if(fPostROIPad[thePlane] > fPreROIPad[thePlane]) 
        minSepPad = fMinSep + fPostROIPad[thePlane];
      
      // skip bad channels
      if(!chanFilt->BadChannel(channel)) {
        
        // uncompress the data
        raw::Uncompress(digitVec->fADC, rawadc, digitVec->Compression());
        // loop over all adc values and subtract the pedestal
        float pdstl = digitVec->GetPedestal();
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
              if(roiLen > transformSize) {
                mf::LogError("CalWireROI")<<"ROI too long "
                  <<roiLen<<" on plane:wire "<<thePlane<<":"<<theWire;
                break;
              }
              if(roiLen > fMinWid && roiLen < transformSize) 
                rois.push_back(std::make_pair(roiStart, bin));
              roiStart = 0;
            }
          } // roiStart test
        } // bin

        // skip deconvolution if there are no ROIs
        if(rois.size() == 0) continue;
        
        // merge the ROIs?
        if(rois.size() > 1) {
          // temporary vector for merged ROIs
          std::vector<std::pair<unsigned int, unsigned int>> trois;
          unsigned int ii = 0;
          unsigned int roiStart = rois[0].first;
          unsigned int roiEnd = rois[0].second;
          bool inMerge = false;
          while(ii < rois.size() - 1) {
            for(unsigned int jj = ii + 1; jj < rois.size(); ++jj) {
              if(rois[jj].first - rois[jj-1].second < minSepPad) {
                // merge them
                ii = jj;
                roiEnd = rois[jj].second;
                inMerge = true;
              } else {
                // don't merge them
                trois.push_back(std::make_pair(roiStart,roiEnd));
                ii = jj;
                roiStart = rois[ii].first;
                roiEnd = rois[ii].second;
                inMerge = false;
                break;
              }
            } // jj
          } // ii
          // finish an unfinished merge
          if(inMerge) trois.push_back(std::make_pair(roiStart,roiEnd));
          // make the last one if necessary
          if(trois.size() > 0) {
            if(trois[trois.size()-1].second != rois[rois.size()-1].second) {
              roiStart = rois[rois.size() - 1].first;
              trois.push_back(std::make_pair(roiStart,roiEnd));
            } // fix up
          } // trois.size() > 0
          rois = trois;
        } // merge ROIs

        // pad the ROIs
        for(unsigned int ii = 0; ii < rois.size(); ++ii) {
          // low ROI end
          int low = rois[ii].first - fPreROIPad[thePlane];
          if(low < 0) low = 0;
          rois[ii].first = low;
          // high ROI end
          unsigned int high = rois[ii].second + fPostROIPad[thePlane];
          if(high > dataSize) high = dataSize;
          rois[ii].second = high;
        }

       unsigned int hBin = 0;
        // max number of holder bins that should be packed. This is the
        // FFT size - width of the response function. Do this crudely for
        // now since the response functions are < 20 time bins in length
        unsigned int packLen = transformSize - 20;
        for(unsigned int ir = 0; ir < rois.size(); ++ir) {
          // enough room to continue packing?
          unsigned int roiLen = rois[ir].second - rois[ir].first;
          if(roiLen > transformSize) {
            mf::LogWarning("CalWireROI")<<"ROI length "<<roiLen
              <<" > FFT size "<<transformSize<<" on plane:wire "
              <<thePlane<<":"<<theWire;
            continue;
          }
          // not enough room to continue packing. Deconvolute this batch, stuff it
          // into ROIVec and continue
          if(hBin + roiLen > packLen) {
            // need to de-convolute holder and continue packing
            doDecon(holder, channel, thePlane, rois, holderInfo, ROIVec, sss);
            holderInfo.clear();
            hBin = 0;
          } // hBin + roiLen > packLen
          // save the position of this ROI in the holder vector and the ROI index
          holderInfo.push_back(std::make_pair(hBin, ir));
          for(unsigned int bin = rois[ir].first; bin < rois[ir].second; ++bin) {
            holder[hBin] = rawadc[bin]-pdstl;
            ++hBin;
          } // bin
        } // ir < rois.size
        // do the last deconvolution if needed
        if(holderInfo.size() > 0)
          doDecon(holder, channel, thePlane, rois, holderInfo, ROIVec, sss);
      } // end if not a bad channel 
      
      // create the new wire directly in wirecol
      wirecol->emplace_back(std::move(ROIVec), digitVec);
      
      DumpWire(wirecol->back());
    }

    if(wirecol->size() == 0)
      mf::LogWarning("CalWireROI") << "No wires made for this event.";

    if(fSpillName.size()>0)
      evt.put(std::move(wirecol), fSpillName);
    else evt.put(std::move(wirecol));

    delete chanFilt;

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
      float baseLine = 0.;
      // Baseline subtraction if requested and the ROIs are padded.
      // Can't baseline subtract signals when the ROI start is close to 0 either
      if(fDoBaselineSub && fPreROIPad[thePlane] > 0 && 
          rois[theROI].first > fPreROIPad[thePlane]) {
        // find the baseline from the first few bins in the leading Pad region
        unsigned short bbins = fPreROIPad[thePlane];
        if(bbins > 5) bbins = 5;
        for(unsigned int bin = 0; bin < bbins; ++bin) 
          baseLine += holder[bBegin + bin];
        baseLine /= (float)bbins;
      }
      for(unsigned int jj = bBegin; jj < bEnd; ++jj) {
        sigTemp.push_back(holder[jj] - baseLine);
      } // jj
      
      // add the range into ROIVec 
      ROIVec.add_range(rois[theROI].first, std::move(sigTemp));
    } // jr
  } // doDecon


} // end namespace caldata
