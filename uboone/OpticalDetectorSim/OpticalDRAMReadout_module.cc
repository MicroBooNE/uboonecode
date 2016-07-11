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
#include "lardata/OpticalDetectorData/FIFOChannel.h"
#include "lardata/OpticalDetectorData/OpticalRawDigit.h"
#include "lardata/RawData/TriggerData.h"
#include "lardata/RawData/OpDetWaveform.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h" // FIXME: this code is non-portable
#include "uboone/OpticalDetectorSim/UBOpticalException.h"
#include "larcore/Geometry/Geometry.h" // larcore
#include "uboone/Geometry/UBOpChannelTypes.h"  // uboonecode
#include "uboone/Geometry/UBOpReadoutMap.h"  // uboonecode

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
    std::string fFIFOModuleName; // (deprecated)

    /// Producer module label for raw::Trigger
    std::string fTrigModuleName;

    // Stem name for data products: defaults to opdrammcreadout
    std::string fDataProductsStemName;

    /// Readout frame offset (w.r.t. trigger frame)
    std::vector<optdata::Frame_t> fReadoutFrameOffset;

    /// Names of data products: one for each Optical Category
    std::map< opdet::UBOpticalChannelCategory_t, std::string > fPMTdataProductNames;
    void registerOpticalData();
    void putPMTDigitsIntoEvent( std::map< opdet::UBOpticalChannelCategory_t, std::unique_ptr< std::vector<raw::OpDetWaveform> > >& pmtdigitlist, art::Event& event);

  };
} // namespace opdet


// Required for any LArSoft module.
namespace opdet {
  DEFINE_ART_MODULE(OpticalDRAMReadout)
} 


// Implementation
namespace opdet {
  
  /// ------------------------------------------------------------------------------------
  /// Constructor
  OpticalDRAMReadout::OpticalDRAMReadout(fhicl::ParameterSet const& parameterSet)
  {

    // Read in the parameters from the .fcl file.
    this->reconfigure(parameterSet);

    // Describe the data products we can write.
    registerOpticalData();

  }

  /// ------------------------------------------------------------------------------------

  void OpticalDRAMReadout::registerOpticalData()
  {
    // we make a data product for each category of channels
    fPMTdataProductNames.clear();
    for ( unsigned int cat=0; cat<(unsigned int)opdet::NumUBOpticalChannelCategories; cat++ ) {
      produces< std::vector<raw::OpDetWaveform> >( opdet::UBOpChannelEnumName( (opdet::UBOpticalChannelCategory_t)cat ) );
      fPMTdataProductNames.insert( std::make_pair( (opdet::UBOpticalChannelCategory_t)cat, opdet::UBOpChannelEnumName( (opdet::UBOpticalChannelCategory_t)cat ) ) );
    }
  }


  /// ------------------------------------------------------------------------------------

  void OpticalDRAMReadout::putPMTDigitsIntoEvent( std::map< opdet::UBOpticalChannelCategory_t, 
						  std::unique_ptr< std::vector<raw::OpDetWaveform> > >& pmtdigitlist, 
						  art::Event& event ) {
    for ( unsigned int cat=0; cat<(unsigned int)opdet::NumUBOpticalChannelCategories; cat++ )
      event.put( std::move( pmtdigitlist[(opdet::UBOpticalChannelCategory_t)cat] ), fPMTdataProductNames[ (opdet::UBOpticalChannelCategory_t)cat ] );    
  }

  /// ------------------------------------------------------------------------------------

  void OpticalDRAMReadout::reconfigure(fhicl::ParameterSet const& p)
  {
    
    fFIFOModuleName     = p.get<std::string> ("FIFOModuleName");

    fTrigModuleName     = p.get<std::string> ("TrigModuleName");
    
    fReadoutFrameOffset = p.get<std::vector<optdata::Frame_t> >("ReadoutFrameOffset");

    fDataProductsStemName = p.get<std::string>( "OpDataProductStemName", "opdrammcreadout");

    // do we want to give user option to name data?
    // do this here

    if(fReadoutFrameOffset.size()!=2)

      throw UBOpticalException("ReadoutFrameOffset must be a vector of length 2!");

    return;
  }

  //-------------------------------------------------

  void OpticalDRAMReadout::produce(art::Event& event)
  {

    // -----------------------------------
    // Get needed services

    // Obtain optical clock to be used for sample/frame number generation
    // FIXME: this code is non-portable
    art::ServiceHandle<detinfo::DetectorClocksServiceStandard> tss;
    tss->preProcessEvent(event); // sets trigger time
    auto const* ts = lar::providerFrom<detinfo::DetectorClocksServiceStandard>();
    ::detinfo::ElecClock clock = ts->OpticalClock();
    //std::cout << "OpticalDRAM: Trigger time=" << ts->TriggerTime() << " Beam gate time=" << ts->BeamGateTime() << std::endl;

    // geometry and channel map services
    ::art::ServiceHandle<geo::Geometry> geom;
    ::art::ServiceHandle<geo::UBOpReadoutMap> ub_pmt_channel_map;

    // -----------------------------------
    // Create out container of waveforms: one container per readout channel type
    std::map< opdet::UBOpticalChannelCategory_t, std::unique_ptr< std::vector<raw::OpDetWaveform> > > pmt_raw_digits;
    for ( unsigned int opdetcat=0; opdetcat<(unsigned int)opdet::NumUBOpticalChannelCategories; opdetcat++ ) {
      pmt_raw_digits.insert( std::make_pair( (opdet::UBOpticalChannelCategory_t)opdetcat, std::unique_ptr< std::vector<raw::OpDetWaveform> >(  new std::vector<raw::OpDetWaveform> ) ) );
    }
    //std::unique_ptr< std::vector<optdata::OpticalRawDigit> > wf_array( new std::vector<optdata::OpticalRawDigit> ); // deprecated

    // -----------------------------------
    // PMT Trigger Sim

    // Check if trigger data product exists or not. If not, throw a warning
    art::Handle< std::vector<raw::Trigger> > trig_array;
    event.getByLabel(fTrigModuleName, trig_array);
    if(!trig_array.isValid()) 

      std::cout << std::endl << "  "
		<< "\033[95m" << "<<" << __PRETTY_FUNCTION__ << ">>" << "\033[00m"
		<< std::endl << "  "
		<< "\033[93m"
		<< " No trigger data exists => will use the default trigger time set in DetectorClocksService..."
		<< "\033[00m"
		<< std::endl;

    // Figure out the range of frame numbers to be readout
    std::vector<std::pair<optdata::Frame_t,optdata::Frame_t> > readout_frames;

    optdata::Frame_t trig_frame = clock.Frame(ts->TriggerTime());
    
    optdata::Frame_t start_frame = 0;
    optdata::Frame_t end_frame = trig_frame + fReadoutFrameOffset.at(1);

    if(trig_frame < fReadoutFrameOffset.at(0)) start_frame = 0;
    else start_frame = trig_frame - fReadoutFrameOffset.at(0);

    readout_frames.push_back(std::make_pair(start_frame,end_frame));

    /*
    //
    // In case of multipl trigger data products ... commented out for now
    //
    if(!trig_array.isValid())

      throw UBOpticalException(Form("Did not find raw::Trigger data product from module %s",fTrigModuleName.c_str()));

    readout_frames.reserve(trig_array->size());
    for(size_t i=0; i<trig_array->size(); ++i) {

      const art::Ptr<raw::Trigger> trig_ptr(trig_array,i);

      optdata::Frame_t trig_frame = clock.Frame(trig_ptr->TriggerTime());

      optdata::Frame_t start_frame = 0;
      optdata::Frame_t end_frame = trig_frame + fReadoutFrameOffset.at(1);

      if(trig_frame < fReadoutFrameOffset.at(0)) start_frame = 0;
      else start_frame = trig_frame - fReadoutFrameOffset.at(0);

      readout_frames.push_back(std::pair<optdata::Frame_t,optdata::Frame_t>(start_frame,end_frame) );
    }
    */   
    
    // -----------------------------------    
    // Read in optdata::FIFOChannel & store relevant portion

    art::Handle< std::vector<optdata::FIFOChannel> > fifo_array;
    event.getByLabel(fFIFOModuleName, fifo_array);

    if(!fifo_array.isValid()) {
      throw UBOpticalException("OpticalDRAMReadout_module::produce. Input data invalid!");
    }
    else {

      // we want to reserve enough space, once. trading speed at the cost of wasting memory
      for ( auto &it : pmt_raw_digits )
	it.second->reserve( fifo_array->size() );


      for(size_t i=0; i<fifo_array->size(); ++i) {

	const art::Ptr<optdata::FIFOChannel> fifo_ptr(fifo_array,i);
	
	bool store=false;
	for(auto const &period : readout_frames)
	  
	  if( period.first <= fifo_ptr->Frame() && fifo_ptr->Frame() < period.second ) {
	    
	    store=true;
	    break;
	    
	  }
	
	if(store){

	  // if FIFO has correct categories, then fine.
	  // but for now, get it using the channel number
	  unsigned int data_product_ch_num = fifo_ptr->ChannelNumber();
	  opdet::UBOpticalChannelCategory_t category = ub_pmt_channel_map->GetChannelCategory( data_product_ch_num );
	  auto it_wfarray = pmt_raw_digits.find( category );
	  double window_timestamp = clock.Time( fifo_ptr->TimeSlice(), fifo_ptr->Frame() );

	  raw::OpDetWaveform rd( window_timestamp, data_product_ch_num, fifo_ptr->size() );
	  for ( size_t iadc=0; iadc< fifo_ptr->size(); iadc++ )
	    rd.push_back( fifo_ptr->at(iadc) );
	  (*it_wfarray).second->emplace_back( rd );
	  // std::cout << "DRAM: Store FIFO from CH=" << data_product_ch_num 
	  // 	    << " put into " << opdet::UBOpChannelEnumName( (opdet::UBOpticalChannelCategory_t)category ) 
	  // 	    << " timestamp=" << window_timestamp 
	  // 	    << " frame=" << fifo_ptr->Frame()
	  // 	    << " timeslice/sample=" << fifo_ptr->TimeSlice()
	  // 	    << " size=" << fifo_ptr->size()
	  // 	    << std::endl;
	}// if store
      }// loop over FIFO data

    }// if valid

    // hand over data to event
    //event.put( std::move( wf_array ) );
    putPMTDigitsIntoEvent( pmt_raw_digits, event );

  }
  
}

#endif
/** @} */ // end of doxygen group 
