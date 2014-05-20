/**
 * \file OpticalDRAMReadout_module.cc
 *
 * \ingroup OpticalDetectorSim
 * 
 * \brief Class def header for a class OpticalDRAMReadout_module.cc
 *
 * @author kazuhiro
 */

/** \addtogroup OpticalDetectorSim

@{*/

#ifndef OPTICALDRAMREADOUT_MODULE_CC
#define OPTICALDRAMREADOUT_MODULE_CC

// LArSoft includes
#include "OpticalDetectorData/FIFOChannel.h"
#include "OpticalDetectorData/OpticalRawDigit.h"
#include "RawData/TriggerData.h"
#include "Utilities/TimeService.h"
#include "uboonecode/uboone/OpticalDetectorSim/UBOpticalException.h"

// ART includes.
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Utilities/Exception.h"

// C++ language includes
#include <iostream>
#include <sstream>
#include <cstring>
#include <vector>
#include <map>
#include <memory>
#include <cmath>

namespace opdet {

  class OpticalDRAMReadout : public art::EDProducer{
  public:
    
    OpticalDRAMReadout(const fhicl::ParameterSet&);
    virtual ~OpticalDRAMReadout() {};
    
    // This method reads in any parameters from the .fcl files. 
    void reconfigure(fhicl::ParameterSet const& pset);
    
    void produce(art::Event&);
      
  private:  

    /// Producer module label for raw::FIFOChannel
    std::string fFIFOModuleName;

    /// Producer module label for raw::Trigger
    std::string fTrigModuleName;

    /// Readout frame offset (w.r.t. trigger frame)
    std::vector<optdata::Frame_t> fReadoutFrameOffset;

  };
} // namespace opdet


// Required for any LArSoft module.
namespace opdet {
  DEFINE_ART_MODULE(OpticalDRAMReadout)
} 


// Implementation
namespace opdet {
  
  OpticalDRAMReadout::OpticalDRAMReadout(fhicl::ParameterSet const& parameterSet)
  {
    // Describe the data products we can write.
    produces< std::vector< optdata::OpticalRawDigit > >();

    // Read in the parameters from the .fcl file.
    this->reconfigure(parameterSet);

  }
  
  void OpticalDRAMReadout::reconfigure(fhicl::ParameterSet const& p)
  {
    fReadoutFrameOffset = p.get<std::vector<optdata::Frame_t> >("ReadoutFrameOffset");

    if(fReadoutFrameOffset.size()!=2)

      throw UBOpticalException("ReadoutFrameOffset must be a vector of length 2!");

    return;
  }

  //-------------------------------------------------

  void OpticalDRAMReadout::produce(art::Event& event)
  {
    // Obtain optical clock to be used for sample/frame number generation
    art::ServiceHandle<util::TimeService> ts;
    ::util::ElecClock clock = ts->OpticalClock();

    std::unique_ptr< std::vector<optdata::OpticalRawDigit> > wf_array( new std::vector<optdata::OpticalRawDigit> );

    // Read in raw::Trigger
    art::Handle< std::vector<raw::Trigger> > trig_array;
    event.getByLabel(fTrigModuleName, trig_array);
    if(!trig_array.isValid()) 

      throw UBOpticalException(Form("Did not find raw::Trigger data product from module %s",fTrigModuleName.c_str()));

    // Store readout frame boundaries
    std::vector<std::pair<optdata::Frame_t,optdata::Frame_t> > readout_frames;
    readout_frames.reserve(trig_array->size());
    for(size_t i=0; i<trig_array->size(); ++i) {

      const art::Ptr<raw::Trigger> trig_ptr(trig_array,i);

      optdata::Frame_t trig_frame = clock.Frame(trig_ptr->TriggerTime());

      readout_frames.push_back(std::pair<optdata::Frame_t,optdata::Frame_t>(trig_frame - fReadoutFrameOffset.at(0),
									    trig_frame + fReadoutFrameOffset.at(1) )
			       );
    }

    // Read in optdata::FIFOChannel & store relevant portion
    art::Handle< std::vector<optdata::FIFOChannel> > fifo_array;
    event.getByLabel(fFIFOModuleName, fifo_array);

    if(fifo_array.isValid()) {

      wf_array->reserve(fifo_array->size());

      for(size_t i=0; i<fifo_array->size(); ++i) {

	const art::Ptr<optdata::FIFOChannel> fifo_ptr(fifo_array,i);
	
	bool store=false;
	for(auto const &period : readout_frames)
	  
	  if( period.first <= fifo_ptr->Frame() && fifo_ptr->Frame() < period.second ) {
	    
	    store=true;
	    break;
	    
	  }
	
	if(store){
	  
	  wf_array->push_back(optdata::OpticalRawDigit(fifo_ptr->Category(),
						       fifo_ptr->TimeSlice(),
						       fifo_ptr->Frame(),
						       fifo_ptr->ChannelNumber(),
						       fifo_ptr->size()
						       )
			      );

	  auto wf = wf_array->rbegin();

	  for(auto const& v : *fifo_ptr)

	    wf->push_back(v);

	}
      }

    }
    
    event.put( std::move( wf_array ) ); 

  }
}

#endif
/** @} */ // end of doxygen group 
