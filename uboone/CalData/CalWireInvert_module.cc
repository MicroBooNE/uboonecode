////////////////////////////////////////////////////////////////////////
//
// CalWireInvert class - variant of CalWire that inverts RawDigit signals
//                       on the W plane. The primary purpose of this
//                       module is to provide inverted waveform ROI's in
//                       which to search for PMT/Laser pickup pulses
// usher@slac.stanford.edu
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
#include "Filters/ChannelFilter.h"
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "RecoBase/Wire.h"
#include "RecoBaseArt/WireCreator.h"
#include "Utilities/AssociationUtil.h"
#include "CalibrationDBI/Interface/IDetPedestalService.h"
#include "CalibrationDBI/Interface/IDetPedestalProvider.h"

///creation of calibrated signals on wires
namespace caldata {

class CalWireInvert : public art::EDProducer
{
public:
    
    // create calibrated signals on wires. this class runs 
    // an fft to remove the electronics shaping.     
    explicit CalWireInvert(fhicl::ParameterSet const& pset); 
    virtual ~CalWireInvert();
    
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
    std::vector<double>         fMinRmsNoise;          ///< minimum noise level for ROI
    std::vector<unsigned short> fThreshold;            ///< abs(threshold) ADC counts for ROI
    std::vector<unsigned short> fPreROIPad;            ///< ROI padding
    std::vector<unsigned short> fPostROIPad;           ///< ROI padding
    size_t                      fEventCount;           ///< count of event processed
    int                         fMaxAllowedChanStatus; ///< Examine channels with status below this
    
    const lariov::IDetPedestalProvider& fPedestalRetrievalAlg; ///< Keep track of an instance to the pedestal retrieval alg
    
    float SubtractBaseline(std::vector<float>& holder, float basePre,
                           float basePost,unsigned int roiStart,unsigned int roiLen,
                           unsigned int dataSize);
    
}; // class CalWireInvert

DEFINE_ART_MODULE(CalWireInvert)
  
//-------------------------------------------------
CalWireInvert::CalWireInvert(fhicl::ParameterSet const& pset) :
    fPedestalRetrievalAlg(art::ServiceHandle<lariov::IDetPedestalService>()->GetPedestalProvider())
{
    this->reconfigure(pset);

    produces< std::vector<recob::Wire> >(fSpillName);
    produces<art::Assns<raw::RawDigit, recob::Wire>>(fSpillName);
}
  
//-------------------------------------------------
CalWireInvert::~CalWireInvert()
{
}

//////////////////////////////////////////////////////
void CalWireInvert::reconfigure(fhicl::ParameterSet const& p)
{
    std::vector<unsigned short> uin;    std::vector<unsigned short> vin;
    std::vector<unsigned short> zin;
    
    fDigitModuleLabel     = p.get< std::string >                   ("DigitModuleLabel", "daq");
    fMinRmsNoise          = p.get< std::vector<double>>            ("MinRMSNoise",      {1., 1., 1.});
    fThreshold            = p.get< std::vector<unsigned short> >   ("Threshold");
    uin                   = p.get< std::vector<unsigned short> >   ("uPlaneROIPad");
    vin                   = p.get< std::vector<unsigned short> >   ("vPlaneROIPad");
    zin                   = p.get< std::vector<unsigned short> >   ("zPlaneROIPad");
    fMaxAllowedChanStatus = p.get< int >                           ("MaxAllowedChannelStatus");
    
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

    return;
}

//-------------------------------------------------
void CalWireInvert::beginJob()
{
    fEventCount = 0;
} // beginJob

//////////////////////////////////////////////////////
void CalWireInvert::endJob()
{  
}

//////////////////////////////////////////////////////
void CalWireInvert::produce(art::Event& evt)
{
    // The following code almost exactly parallels the code in CalWireROI for finding
    // Regions of Interest to snip out. It removes the deconvolution as we're only looking
    // at the collection plane here and should not need it (for looking at PMT pickup).
    // get the geometry
    art::ServiceHandle<geo::Geometry> geom;
  
    // make a collection of Wires
    std::unique_ptr<std::vector<recob::Wire> > wirecol(new std::vector<recob::Wire>);
    // ... and an association set
    std::unique_ptr<art::Assns<raw::RawDigit,recob::Wire> > WireDigitAssn(new art::Assns<raw::RawDigit,recob::Wire>);
  
    // Read in the digit List object(s).
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
    if(fSpillName.size()>0) evt.getByLabel(fDigitModuleLabel, fSpillName, digitVecHandle);
    else                    evt.getByLabel(fDigitModuleLabel, digitVecHandle);

    if (digitVecHandle->size() > 0)
    {
        raw::ChannelID_t channel = raw::InvalidChannelID; // channel number
        unsigned int bin(0);     // time bin loop variable
  
        std::unique_ptr<filter::ChannelFilter> chanFilt(new filter::ChannelFilter());
  
        // loop over all wires
        wirecol->reserve(digitVecHandle->size());
  
        for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter)
        {
            // vector that will be moved into the Wire object
            recob::Wire::RegionsOfInterest_t ROIVec;
    
            // the starting position and length of each ROI in the packed holder vector
            std::vector<std::pair<unsigned int, unsigned int>> holderInfo;
            // vector of ROI begin and end bins
            std::vector<std::pair<unsigned int, unsigned int>> rois;
    
            // get the reference to the current raw::RawDigit
            art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
            channel = digitVec->Channel();

            // The following test is meant to be temporary until the "correct" solution is implemented
            if (chanFilt->GetChannelStatus(channel) == filter::ChannelFilter::NOTPHYSICAL) continue;
      
            // Testing an idea about rejecting channels
            if (digitVec->GetPedestal() < 0.) continue;

            unsigned int dataSize = digitVec->Samples();
            // vector holding uncompressed adc values
            std::vector<short> rawadc(dataSize);
    
            std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
            unsigned int thePlane = wids[0].Plane;
            //unsigned int theWire = wids[0].Wire;
      
            // Skip induction planes
            if (thePlane != 2) continue;
    
            // skip bad channels
            if(!(chanFilt->GetChannelStatus(channel) > fMaxAllowedChanStatus))
            {
                // uncompress the data
                raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
                // loop over all adc values and subtract the pedestal
                // When we have a pedestal database, can provide the digit timestamp as the third argument of GetPedestalMean
                float        pdstl    = fPedestalRetrievalAlg.PedMean(channel);
                unsigned int roiStart = 0;
        
                double rms_noise = digitVec->GetSigma();
                double raw_noise = std::max(fMinRmsNoise[thePlane],rms_noise);

                // search for ROIs
                for(bin = 1; bin < dataSize; ++bin)
                {
                    float SigVal = fabs(rawadc[bin] - pdstl);
                
                    if(roiStart == 0)
                    {
                        // not in a ROI
                        // Handle the onset of a ROI differently for the 1st induction plane
                        // if it has been modeled correctly
                        unsigned int sbin[7];
                  
                        if (bin>=3)
                        {
                            sbin[0] = bin -3;
                            sbin[1] = bin -2;
                            sbin[2] = bin -1;
                        }
                        else if (bin>=2)
                        {
                            sbin[0] = 0;
                            sbin[1] = bin-2;
                            sbin[2] = bin-1;
                        }
                        else if (bin>=1)
                        {
                            sbin[0] =0;
                            sbin[1] =0;
                            sbin[2] = bin-1;
                        }
                        else if (bin==0)
                        {
                            sbin[0] =0;
                            sbin[1] =0;
                            sbin[2] =0;
                        }
                  
                        sbin[3] = bin ;
                        sbin[4] = bin + 1; if (sbin[4]>dataSize-1) sbin[4] =dataSize-1;
                        sbin[5] = bin + 2; if (sbin[5]>dataSize-1) sbin[5] =dataSize-1;
                        sbin[6] = bin + 3; if (sbin[6]>dataSize-1) sbin[6] =dataSize-1;
                        float sum = 0;
                        for (int qx = 0; qx!=7;qx++)
                            sum += rawadc[sbin[qx]]-pdstl;
                    
                        sum = fabs(sum);

                        if (sum > raw_noise*sqrt(7.)*6.) roiStart = bin;
                    }
                    else
                    {
                        // leaving a ROI?
                        if(SigVal < fThreshold[thePlane])
                        {
                            // is the ROI wide enough?
                            rois.push_back(std::make_pair(roiStart, bin));
                            roiStart = 0;
                        }
                    } // roiStart test
                } // bin
          
                // add the last ROI if existed
                if (roiStart!=0)
                {
                    //unsigned int roiLen = dataSize -1 - roiStart;
                    rois.push_back(std::make_pair(roiStart, dataSize-1));
                    roiStart = 0;
                }

                // skip deconvolution if there are no ROIs
                if(rois.size() == 0) continue;

                holderInfo.clear();
        
                // pad the ROIs
                for(unsigned int ii = 0; ii < rois.size(); ++ii)
                {
                    // low ROI end
                    int low = rois[ii].first - fPreROIPad[thePlane];
                    if(low < 0) low = 0;
                    rois[ii].first = low;
                    // high ROI end
                    unsigned int high = rois[ii].second + fPostROIPad[thePlane];
                    if(high >= dataSize) high = dataSize-1;
                    rois[ii].second = high;
                }

                // merge the ROIs?
                if(rois.size() > 1)
                {
                    // temporary vector for merged ROIs
                    std::vector<std::pair<unsigned int, unsigned int>> trois;
          
                    for (unsigned int ii = 0; ii<rois.size();ii++)
                    {
                        unsigned int roiStart = rois[ii].first;
                        unsigned int roiEnd = rois[ii].second;

                        int flag1 = 1;
                        unsigned int jj=ii+1;
                        while(flag1)
                        {
                            if (jj<rois.size())
                            {
                                if(rois[jj].first <= roiEnd  )
                                {
                                    roiEnd = rois[jj].second;
                                    ii = jj;
                                    jj = ii+1;
                                }
                                else flag1 = 0;
                            }
                            else flag1 = 0;
                        }
	 
                        trois.push_back(std::make_pair(roiStart,roiEnd));
                    }
	  
                    rois = trois;
                }

                for (unsigned int ir = 0; ir < rois.size(); ++ir)
                {
                    unsigned int roiLen = rois[ir].second - rois[ir].first + 1;
                    unsigned int roiStart = rois[ir].first;
                    float        tempPre=0,tempPost=0;
                    std::vector<float> holder;
                
                    holder.resize(roiLen);
                
                    size_t hBin(0);
                
                    for(size_t roiIdx = roiStart; roiIdx < roiStart + roiLen; roiIdx++)
                        holder[hBin++] = -(rawadc[roiIdx] - pdstl);
	  
                    // transfer the ROI and start bins into the vector that will be
                    // put into the event
                    std::vector<float> sigTemp;
                    unsigned int bBegin = 0;
                    //unsigned int theROI =ir;
                    unsigned int bEnd = bBegin + roiLen;
                    float basePre = 0., basePost = 0.;

                    float base=0;
                    if(fPreROIPad[thePlane] > 0 )
                    {
                        basePre =tempPre;
                        basePost=tempPost;

                        base = SubtractBaseline(holder, basePre,basePost,roiStart,roiLen,dataSize);
                    } // fDoBaselineSub ...

                    for(unsigned int jj = bBegin; jj < bEnd; ++jj)
                    {
                        sigTemp.push_back(holder[jj]-base);
                    } // jj
	        
                    // add the range into ROIVec
                    ROIVec.add_range(roiStart, std::move(sigTemp));
                }
            } // end if not a bad channel

            // create the new wire directly in wirecol
            wirecol->push_back(recob::WireCreator(std::move(ROIVec),*digitVec).move());
            // add an association between the last object in wirecol
            // (that we just inserted) and digitVec
            if (!util::CreateAssn(*this, evt, *wirecol, digitVec, *WireDigitAssn, fSpillName))
            {
                throw art::Exception(art::errors::InsertFailure)
                    << "Can't associate wire #" << (wirecol->size() - 1)
                    << " with raw digit #" << digitVec.key();
            } // if failed to add association
        //  DumpWire(wirecol->back()); // for debugging
        }
    }

    if(wirecol->size() == 0)
        mf::LogWarning("CalWireInvert") << "No wires made for this event.";
    
    evt.put(std::move(wirecol), fSpillName);
    evt.put(std::move(WireDigitAssn), fSpillName);

    fEventCount++;

    return;
} // produce


float CalWireInvert::SubtractBaseline(std::vector<float>& holder, float basePre,
                                      float basePost,unsigned int roiStart,
                                      unsigned int roiLen,unsigned int dataSize)
{
    float base=0;

    //can not trust the early part
    if (roiStart < 20 && roiStart + roiLen < dataSize - 20)
    {
        base = basePost;
    }
    // can not trust the later part
    else if (roiStart >= 20 && roiStart + roiLen >= dataSize - 20)
    {
        base = basePre;
    }
    // can trust both
    else if (roiStart >= 20 && roiStart + roiLen < dataSize - 20)
    {
        if (fabs(basePre-basePost)<3)
        {
            base = (basePre+basePost)/2.;
        }
        else
        {
            if (basePre < basePost)
            {
                base = basePre;
            }
            else
            {
                base = basePost;
            }
        }
    }
    // can not use both
    else
    {
        float min = 0,max=0;
        for (unsigned int bin = 0; bin < roiLen; bin++)
        {
            if (holder[bin] > max) max = holder[bin];
            if (holder[bin] < min) min = holder[bin];
        }
        int nbin = max - min;
        if (nbin!=0)
        {
            TH1F *h1 = new TH1F("h1","h1",nbin,min,max);
            for (unsigned int bin = 0; bin < roiLen; bin++)
            {
                h1->Fill(holder[bin]);
            }
            float ped = h1->GetMaximum();
            float ave=0,ncount = 0;
            for (unsigned int bin = 0; bin < roiLen; bin++)
            {
                if (fabs(holder[bin]-ped)<2)
                {
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


} // end namespace caldata
