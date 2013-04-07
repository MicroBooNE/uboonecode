////////////////////////////////////////////////////////////////////////
//
// CalWireMicroBooNE class
//
// brebel@fnal.gov
//
// 11-3-09 Pulled all FFT code out and put into Utilitiess/LArFFT
//
////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <stdint.h>

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "cetlib/exception.h"
#include "cetlib/search_path.h"

#include "Utilities/SignalShapingServiceMicroBooNE.h"
#include "Geometry/Geometry.h"
#include "Filters/ChannelFilter.h"
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "RecoBase/Wire.h"
#include "Utilities/LArFFT.h"

#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"

///creation of calibrated signals on wires
namespace caldata {

  class CalWireMicroBooNE : public art::EDProducer {

  public:
    
    // create calibrated signals on wires. this class runs 
    // an fft to remove the electronics shaping.     
    explicit CalWireMicroBooNE(fhicl::ParameterSet const& pset); 
    virtual ~CalWireMicroBooNE();
    
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob();                 
    void reconfigure(fhicl::ParameterSet const& p);
 
  private:
    
    int          fDataSize;          ///< size of raw data on one wire
    int          fPostsample;        ///< number of postsample bins
    std::string  fDigitModuleLabel;  ///< module that made digits
                                                       ///< constants
    std::string  fSpillName;  ///< nominal spill is an empty string
                              ///< it is set by the DigitModuleLabel
                              ///< ex.:  "daq:preSpill" for prespill data

  protected: 
    
  }; // class CalWireMicroBooNE

  DEFINE_ART_MODULE(CalWireMicroBooNE);
  
  //-------------------------------------------------
  CalWireMicroBooNE::CalWireMicroBooNE(fhicl::ParameterSet const& pset)
  {
    fSpillName="";
    this->reconfigure(pset);

    if(fSpillName.size()<1) produces< std::vector<recob::Wire> >();
    else produces< std::vector<recob::Wire> >(fSpillName);
  }
  
  //-------------------------------------------------
  CalWireMicroBooNE::~CalWireMicroBooNE()
  {
  }

  //////////////////////////////////////////////////////
  void CalWireMicroBooNE::reconfigure(fhicl::ParameterSet const& p)
  {
    fDigitModuleLabel = p.get< std::string >("DigitModuleLabel", "daq");
    fPostsample       = p.get< int >        ("PostsampleBins");
    
    fSpillName="";
    
    size_t pos = fDigitModuleLabel.find(":");
    if( pos!=std::string::npos ) {
      fSpillName = fDigitModuleLabel.substr( pos+1 );
      fDigitModuleLabel = fDigitModuleLabel.substr( 0, pos );
    }
    
  }

  //-------------------------------------------------
  void CalWireMicroBooNE::beginJob()
  {  
  }

  //////////////////////////////////////////////////////
  void CalWireMicroBooNE::endJob()
  {  
  }
  
  //////////////////////////////////////////////////////
  void CalWireMicroBooNE::produce(art::Event& evt)
  {      
    // get the geometry
    art::ServiceHandle<geo::Geometry> geom;

    // get the FFT service to have access to the FFT size
    art::ServiceHandle<util::LArFFT> fFFT;
    int transformSize = fFFT->FFTSize();

    // Get signal shaping service.
    art::ServiceHandle<util::SignalShapingServiceMicroBooNE> sss;

    // make a collection of Wires
    std::unique_ptr<std::vector<recob::Wire> > wirecol(new std::vector<recob::Wire>);
    
    // Read in the digit List object(s). 
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
    if(fSpillName.size()>0) evt.getByLabel(fDigitModuleLabel, fSpillName, digitVecHandle);
    else evt.getByLabel(fDigitModuleLabel, digitVecHandle);

    if (!digitVecHandle->size())  return;
    mf::LogInfo("CalWireMicroBooNE") << "CalWireMicroBooNE:: digitVecHandle size is " << digitVecHandle->size();

    // Use the handle to get a particular (0th) element of collection.
    art::Ptr<raw::RawDigit> digitVec0(digitVecHandle, 0);
        
    unsigned int dataSize = digitVec0->Samples(); //size of raw data vectors
    
    uint32_t     channel(0); // channel number
    unsigned int bin(0);     // time bin loop variable
    
    filter::ChannelFilter *chanFilt = new filter::ChannelFilter();  

    std::vector<float> holder;                // holds signal data
    std::vector<short> rawadc(transformSize);  // vector holding uncompressed adc values
    std::vector<TComplex> freqHolder(transformSize+1); // temporary frequency data
    
    // loop over all wires    
    for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter){ // ++ move
      holder.clear();
      
      // get the reference to the current raw::RawDigit
      art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
      channel = digitVec->Channel();

      // skip bad channels
      if(!chanFilt->BadChannel(channel)) {
	holder.resize(transformSize);
	
	// uncompress the data
	raw::Uncompress(digitVec->fADC, rawadc, digitVec->Compression());
	
	// loop over all adc values and subtract the pedestal
	for(bin = 0; bin < dataSize; ++bin) 
	  holder[bin]=(rawadc[bin]-digitVec->GetPedestal());

	// Do deconvolution.
	sss->Deconvolute(channel, holder);

      } // end if not a bad channel 
      
      holder.resize(dataSize,1e-5);

      //This restores the DC component to signal removed by the deconvolution.
      if(fPostsample) {
        double average=0.0;
	for(bin=0; bin < (unsigned int)fPostsample; ++bin) 
	  average+=holder[holder.size()-1-bin]/(double)fPostsample;
        for(bin = 0; bin < holder.size(); ++bin) holder[bin]-=average;
      }  
      wirecol->push_back(recob::Wire(holder,digitVec));
    }
    
    if(wirecol->size() == 0)
      mf::LogWarning("CalWireMicroBooNE") << "No wires made for this event.";

    if(fSpillName.size()>0)
      evt.put(std::move(wirecol), fSpillName);
    else evt.put(std::move(wirecol));


    delete chanFilt;
    return;
  }
  
} // end namespace caldata
