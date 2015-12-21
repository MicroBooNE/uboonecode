// \file OpticalFEM_module.cc 
// \author William Seligman <seligman@nevis.columbia.edu>
//
// This module models the behavior of the MicroBooNE PMT Front End
// Modules (FEMs). It reads in the output of the PMTs as ADC counts,
// and outputs the readout as produced by the on-line system.

// The logic of this routine is described in
// <http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=2465>

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/OpticalDetectorData/OpticalTypes.h"
#include "lardata/OpticalDetectorData/ChannelData.h"
#include "lardata/OpticalDetectorData/ChannelDataGroup.h"
#include "lardata/OpticalDetectorData/FIFOChannel.h"
#include "lardata/OpticalDetectorData/PMTTrigger.h"
#include "larsim/Simulation/BeamGateInfo.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "UBOpticalChConfig.h"
#include "UBOpticalConstants.h"
#include "uboone/Geometry/UBOpChannelTypes.h"
#include "uboone/Geometry/UBOpReadoutMap.h"

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

// ROOT includes (for diagnostic histograms)
#include <TH1S.h>
#include <TH1D.h>

// C++ language includes
#include <iostream>
#include <sstream>
#include <cstring>
#include <vector>
#include <map>
#include <memory>
#include <cmath>

namespace opdet {

  class OpticalFEM : public art::EDProducer{
  public:
    
    OpticalFEM(const fhicl::ParameterSet&);
    virtual ~OpticalFEM() {};
    
    // This method reads in any parameters from the .fcl files. 
    void reconfigure(fhicl::ParameterSet const& pset);
    
    void produce(art::Event&);
      
  private:  
    // The parameters we'll read from the .fcl file. There are two numbers associated with
    // most parameters: the low-gain and the high-gain; these are stored in two-element vectors.
    typedef std::vector<optdata::TimeSlice_t> timeVector_t;
    typedef std::vector<unsigned int>         adcVector_t; // ART can't handle reading vectors of ADC_Count_t
    std::string          fBeamModule;           // module that created input simulated beam gate
    std::string          fFakeBeamModule;       // module that created input "fake" beam gate (UBOpticalADCSim)
    std::string          fInputModule;          // module that created input ADC counts.
    timeVector_t         fm_delay0;             // delay for DIFF subtraction, in time slices
    std::vector<int>     fm_delay1;             // number of slices to go back relative to disc 0
    timeVector_t         fm_disc0window;        // disc 0 must have fired within this many slices for disc 1 to fire
    timeVector_t         fm_cosmicSlices;       // number of slices to write when disc 1 fires
    adcVector_t          fm_beamThreshold;      // total ADC counts among selected channels to satisfy PMT beam trigger
    adcVector_t          fm_cosmicThreshold;    // total ADC counts among selected channels to satisfy PMT cosmic trigger
    std::vector<short>   fm_beamMultiplicity;   // number of channels that must satisify conditions for PMT beam trigger
    std::vector<short>   fm_cosmicMultiplicity; // number of channels that must satisify conditions for PMT cosmic trigger
    adcVector_t          fm_threshold0;         // lower limit for disc 0, in ADC counts
    adcVector_t          fm_threshold1;         // lower limit for disc 1, in ADC counts
    adcVector_t          fm_threshold3;         // lower limit for disc 3, in ADC counts
    timeVector_t         fm_disc0quiet;         // quiet interval required between successive disc 0, in time slices
    timeVector_t         fm_disc1deadtime;      // dead time for disc 1, in time slices
    timeVector_t         fm_disc3deadtime;      // dead time for disc 3, in time slices
    timeVector_t         fm_disc1width;         // pulse width for disc 1, in time slices
    timeVector_t         fm_disc3width;         // pulse width for disc 3, in time slices
    timeVector_t         fm_beamWordsBNB;       // number of time slices to save for each BNB beam gate
    timeVector_t         fm_beamWordsNuMI;      // number of time slices to save for each NuMI beam gate
    timeVector_t         fm_beamDelayBNB;       // number of time slices to include prior to start of BNB beam gate
    timeVector_t         fm_beamDelayNuMI;      // number of time slices to include prior to start of NuMI beam gate
    timeVector_t         fm_triggerDeadtime;    // minimum number of time slices between successive PMT triggers
    std::vector<short>   fm_triggerFEMSlot;     // slot number for FEM modules to be used for generating PMT triggers
    bool                 fm_hist;               // if true, generate megatons of diagnostic histograms

    void FillBeamTimingVectors(const std::vector<sim::BeamGateInfo> &beamGates,
			       const optdata::TimeSlice_t readout_first_tdc,
			       const optdata::TimeSlice_t readout_size);
    // Determine the "begin" and "end" bin of the beam-gate windows relative to the first channel time slice.
    std::vector< std::vector< optdata::TimeSlice_t > > _beamBeginBin;
    std::vector< std::vector< optdata::TimeSlice_t > > _beamEndBin;
    // It's also useful to compute the frame number of the start of each beam gate.
    std::vector< std::vector< optdata::Frame_t     > > _gateFrame;
    // The time of first slice in the gate window
    std::vector< std::vector< optdata::TimeSlice_t > > _gateWindowTime;

  };
} // namespace opdet


// Required for any LArSoft module.
namespace opdet {
  DEFINE_ART_MODULE(OpticalFEM)
} 


// Implementation
namespace opdet {
  
  OpticalFEM::OpticalFEM(fhicl::ParameterSet const& parameterSet)
  {
    // Describe the data products we can write.
    produces< std::vector< optdata::FIFOChannel > >();
    produces< std::vector< optdata::PMTTrigger > >();

    // Read in the parameters from the .fcl file.
    this->reconfigure(parameterSet);

  }
  
  void OpticalFEM::reconfigure(fhicl::ParameterSet const& p)
  {
    // Read parameters from the .fcl file.
    fInputModule          = p.get< std::string >         ("OpticalDigitizationModule");
    fBeamModule           = p.get< std::string >         ("BeamGateModule");
    fFakeBeamModule       = p.get< std::string >         ("FakeBeamGateModule");
    fm_hist               = p.get< bool >                ("VerboseHistograms");
    fm_delay0             = p.get< timeVector_t >        ("PMTDelay0");
    fm_delay1             = p.get< std::vector<int> >    ("PMTDelay1");
    fm_disc0window        = p.get< timeVector_t >        ("Discriminator0Window");
    fm_cosmicSlices       = p.get< timeVector_t >        ("PMTWords");
    fm_beamWordsBNB       = p.get< timeVector_t >        ("BeamWordsBNB");
    fm_beamWordsNuMI      = p.get< timeVector_t >        ("BeamWordsNuMI");
    fm_beamDelayBNB       = p.get< timeVector_t >        ("BeamDelayBNB");
    fm_beamDelayNuMI      = p.get< timeVector_t >        ("BeamDelayNuMI");
    fm_beamThreshold      = p.get< timeVector_t >        ("BeamThreshold");
    fm_cosmicThreshold    = p.get< timeVector_t >        ("CosmicThreshold");
    fm_beamMultiplicity   = p.get< std::vector<short> >  ("BeamMultiplicity");
    fm_cosmicMultiplicity = p.get< std::vector<short> >  ("CosmicMultiplicity");
    fm_threshold0         = p.get< adcVector_t >         ("DiscriminatorThreshold0");
    fm_threshold1         = p.get< adcVector_t >         ("DiscriminatorThreshold1");
    fm_threshold3         = p.get< adcVector_t >         ("DiscriminatorThreshold3");
    fm_disc0quiet         = p.get< timeVector_t >        ("PMTPrecount");
    fm_disc1deadtime      = p.get< timeVector_t >        ("Discriminator1DeadTime");
    fm_disc3deadtime      = p.get< timeVector_t >        ("Discriminator3DeadTime");
    fm_disc1width         = p.get< timeVector_t >        ("Discriminator1Width");
    fm_disc3width         = p.get< timeVector_t >        ("Discriminator3Width");
    fm_triggerDeadtime    = p.get< timeVector_t >        ("PMTTriggerDeadtime");
    fm_triggerFEMSlot     = p.get< std::vector<short> >  ("TriggerFEMSlot");

    return;
  }

  //-------------------------------------------------
  void OpticalFEM::FillBeamTimingVectors(const std::vector<sim::BeamGateInfo> &beamGates,
					 const optdata::TimeSlice_t readout_start_tdc,
					 const optdata::TimeSlice_t readout_size)
  {
    
    // Obtain optical clock to be used for sample/frame number generation
    auto const* ts = lar::providerFrom<detinfo::DetectorClocksService>();
    ::detinfo::ElecClock clock = ts->OpticalClock();
    size_t numberOfGates = beamGates.size();

    // Determine the "begin" and "end" bin of the beam-gate
    // windows relative to the first channel time slice.
    _beamBeginBin.clear();
    _beamEndBin.clear();
    _gateFrame.clear();
    _gateWindowTime.clear();

    //std::cout << "FillBeamTimingVectors: " << "beam_bnb_delay=" << fm_beamDelayBNB.size() << " ngates=" << numberOfGates << std::endl;

    for(size_t gain_index=0; gain_index < fm_beamDelayBNB.size(); ++gain_index){

      _beamBeginBin.push_back(std::vector<optdata::TimeSlice_t>(numberOfGates,0));
      _beamEndBin.push_back(std::vector<optdata::TimeSlice_t>(numberOfGates,0));
      _gateFrame.push_back(std::vector<optdata::TimeSlice_t>(numberOfGates,0));
      _gateWindowTime.push_back(std::vector<optdata::TimeSlice_t>(numberOfGates,0));

      auto& beamBeginBin = _beamBeginBin.back();
      auto& beamEndBin   = _beamEndBin.back();
      auto& gateFrame    = _gateFrame.back();
      auto& gateWindowTime   = _gateWindowTime.back();

      for ( size_t gateIndex = 0; gateIndex < 2 && gateIndex < numberOfGates; ++gateIndex ) {
	
	// Fetch the start and width of the beam gate.
	const sim::BeamGateInfo& beamGateInfo = beamGates.at(gateIndex);
	
	if(ts->G4ToElecTime(beamGateInfo.Start()) < 0)
	  
	  throw cet::exception("OpticalFEM") 
	    << "\033[93m"
	    << " Found BeamGateInfo @ " << beamGateInfo.Start() << " [ns] (G4 time)"
	    << " which is " << ts->G4ToElecTime(beamGateInfo.Start()) << " [us] (Elec. time)"
	    << " ... aborting."
	    << "\033[00m" << std::endl;
	
	optdata::TimeSlice_t gateTime = ts->OpticalG4Time2TDC(beamGateInfo.Start());
	//optdata::TimeSlice_t gateWidth = clock.Ticks(beamGateInfo.Width());

	if(gateTime < readout_start_tdc) {
	  //std::cout << "FillBeamTimingVectors: gateTime < readout_start_tdc (" << gateTime << " < " << readout_start_tdc << ")" << std::endl;
	  continue;
	}

	// Figure out the first bin we should start to save.
	//optdata::TimeSlice_t firstSlice = (*channelDataHandle).TimeSlice();
	beamBeginBin[gateIndex] = 0;
	optdata::TimeSlice_t beam_delay = 0;
	optdata::TimeSlice_t beam_words = 0;
	switch( beamGateInfo.BeamType()) {
	case ::sim::kBNB:
	  beam_delay = fm_beamDelayBNB.at(gain_index);
	  beam_words = fm_beamWordsBNB.at(gain_index);
	  break;
	case ::sim::kNuMI:
	  beam_delay = fm_beamDelayNuMI.at(gain_index);
	  beam_words = fm_beamWordsNuMI.at(gain_index);
	  break;
	default:
	  throw cet::exception("OpticalFEM") 
	    << "\033[93m"
	    << "Unsupported Beam Type: " << beamGateInfo.BeamType()
	    << "\033[00m" << std::endl;
	}
	
	if( gateTime < beam_delay )
	  
	  throw cet::exception("OpticalFEM") 
	    << "\033[93m"
	    << " Beam delay " << beam_delay 
	    << " is larger than BeamGateStart time slice " << gateTime
	    << " ... aborting."
	    << "\033[00m" << std::endl;
	
	//beamBeginBin[gateIndex] = ( gateTime + beam_delay - readout_start_tdc); // tmw: is this a bug?
	if ( (int)gateTime-(int)beam_delay-(int)readout_start_tdc<0 )
	  beamBeginBin[gateIndex] = 0;
	else
	  beamBeginBin[gateIndex] = gateTime - beam_delay - readout_start_tdc;

	beamEndBin[gateIndex] = beamBeginBin[gateIndex] + beam_words;
	if(beamEndBin[gateIndex] > readout_size)
	  beamEndBin[gateIndex] = readout_size;
	
	// Compute the frame number of the first slice to be saved. 
	gateFrame[gateIndex] = ( beamBeginBin[gateIndex] / clock.FrameTicks() );
	// The time of the "beginBin[gateIndex]" slice within the gateFrame.
	gateWindowTime[gateIndex] = ( beamBeginBin[gateIndex] % clock.FrameTicks() ); 
      
	LOG_DEBUG("OpticalFEM") 
	//std::cout << "Beam gate #" << gateIndex
	  << " begin beam gate bin to save = " << beamBeginBin[gateIndex]
	  << "; end beam gate bin to save = " << beamEndBin[gateIndex]
	  << "; beam gate frame = " << gateFrame[gateIndex]
	  << "; starts at sample = " << gateWindowTime[gateIndex]
	  << "; slices to save = " << beam_words;
      }
    }

  }

  //-------------------------------------------------

  void OpticalFEM::produce(art::Event& event)
  {
    // Obtain optical clock to be used for sample/frame number generation
    auto const* ts = lar::providerFrom<detinfo::DetectorClocksService>();
    ::detinfo::ElecClock clock = ts->OpticalClock();

    // The collection of channels we'll write in response to beam
    // gates and cosmic signals.
    std::unique_ptr< std::vector<optdata::FIFOChannel> > 
      channelCollection (new std::vector<optdata::FIFOChannel>);

    // The collection of triggers we'll write if the beam/cosmic PMT
    // trigger conditions are satisified.
    std::unique_ptr< std::vector<optdata::PMTTrigger> > 
      triggerCollection (new std::vector<optdata::PMTTrigger>);

    //
    // Read in optdata::ChannelDataGroup. Without these, there's nothing
    // for this module to do, so use getByLabel. 
    //
    art::Handle< optdata::ChannelDataGroup > channelDataHandle;
    event.getByLabel(fInputModule, channelDataHandle);

    // If no channel data group found, provide warning an return
    if(!channelDataHandle.isValid()){

      std::cout << std::endl
		<<"\033[93m" 
		<< "<<" << __PRETTY_FUNCTION__ << ">>" << "  No channel data found. Skipping...."
		<< "\033[00m" 
		<< std::endl
		<< std::endl;
      return;
    }

    
    // Create a work vector for the channel diff subtraction. To save
    // on execution time and memory use, make a rough estimate to the
    // size of this vector and allocate it now. That estimate will be
    // the size of the the first channel of the first group.
    auto const sizeFirstChannel = ((*channelDataHandle).at(0)).size();
    std::vector< optdata::ADC_Count_t > diffVector;
    diffVector.reserve( sizeFirstChannel );
    auto const firstSlice = (*channelDataHandle).TimeSlice();
    auto const firstFrame = (*channelDataHandle).Frame();
    auto const firstTDC   = firstSlice + firstFrame * clock.FrameTicks();
    //
    // std::cout << " Figure out beamgate timings " << std::endl;
    //
    art::Handle< std::vector<sim::BeamGateInfo> > beamGates;
    event.getByLabel(fBeamModule, beamGates);
    // Did we actually read in any BeamGateInfo objects?
    size_t numberOfGates = 0;
    if ( beamGates.isValid() ) {
      if( !fFakeBeamModule.empty() ) {
	std::cout<< "\033[95m[ERROR]\033[00m Found both BeamGateModule and FakeBeamGateModule provided!" << std::endl;
	throw std::exception();
      }
      numberOfGates = beamGates->size();
      FillBeamTimingVectors(*beamGates,
			    firstTDC,
			    (firstSlice + sizeFirstChannel - 1) );
    }
    
    else if(!fFakeBeamModule.empty()) {
      event.getByLabel(fFakeBeamModule, beamGates);
      if( beamGates.isValid() ) {
	numberOfGates = beamGates->size();
	FillBeamTimingVectors(*beamGates,
			      firstTDC,
			      (firstSlice + sizeFirstChannel -1) );
      }
    }
  
    //std::cout << "Number of Gates: " << numberOfGates << std::endl;

    // Do the same for the vectors that will accumulate the sums of
    // the ADC and multiplicity counts for the PMT trigger processing.
    std::vector< std::vector< optdata::ADC_Count_t > > maxADCSum1;
    std::vector< std::vector< optdata::ADC_Count_t > > maxADCSum3;
    std::vector< std::vector< short                > > multiplicitySum1;
    std::vector< std::vector< short                > > multiplicitySum3;
    std::vector< size_t > slot_to_gain_index;

    //
    // Channel map algorithm
    //
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<geo::UBOpReadoutMap> ch_map;
    art::ServiceHandle<opdet::UBOpticalChConfig> ch_conf;

    // This two-dimensional map contains, for each
    // [channel][discriminator] entry, a list of time slices for
    // which the corresponding discriminator has fired.
    // "Discriminator" can be 0, 1, or 3; let's save a teeny
    // amount of space by giving it a type of "short"
    typedef std::vector< optdata::TimeSlice_t > fireList_t;
    std::map< optdata::Channel_t, std::map< short, fireList_t > > fireMap;
    
    // For every disc 1 or disc 3 entry in the above map, the
    // following map contains the max ADC count value during the
    // discriminator width.
    typedef std::vector< optdata::TimeSlice_t > maxADCList_t;
    std::map< optdata::Channel_t, std::map< short, maxADCList_t > > maxADCMap;

    // A list of the type and time slices of PMT triggers issued
    // by the FEM.
    typedef std::pair< optdata::Optical_Category_t, optdata::TimeSlice_t> trigger_t;
    typedef std::vector< trigger_t > triggers_t;
    triggers_t triggers;

    // For each ChannelDataGroup object in this event:
    //for ( auto const& channelDataGroup : (*channelDataHandle) ) {

    if ( numberOfGates>0 ) {
      for ( auto const& channelData : (*channelDataHandle) ) {

	optdata::Channel_t channel = channelData.ChannelNumber();

	::opdet::UBOpticalChannelType_t gain_type = ch_map->GetChannelType(channelData.ChannelNumber());
	size_t gain_index = 0;
	// Determine the gain category.
	switch(gain_type) {
	case ::opdet::LowGain:
	  gain_index=0; break;
	case ::opdet::HighGain:
	case ::opdet::LogicChannel:
	  gain_index=1; break;
	  //gain_index=2; break;
	default:
	  mf::LogError("OpticalFEM") 
	    << "Unknown channel data category = " <<  gain_type
	    << "; skipping channel data group";
	  continue; // skip to the next ChannelDataGroup
	}

	auto const& gateFrame      = _gateFrame.at      (gain_index);
	auto const& gateWindowTime = _gateWindowTime.at (gain_index);
	auto const& beamBeginBin   = _beamBeginBin.at   (gain_index);
	auto const& beamEndBin     = _beamEndBin.at     (gain_index);

	// For each beam gate...
	for ( size_t gateIndex = 0; 
	      //(gain_type == ::opdet::LowGain || gain_type == ::opdet::HighGain) && gateIndex < numberOfGates; 
	      gateIndex < numberOfGates;
	      ++gateIndex ) {
		  
	  // Create a new FIFOChannel, copying the channel
	  // number from the input channel, with length
	  // beam_words.

	  if( (beamEndBin[gateIndex] - beamBeginBin[gateIndex]) < 1 ) continue;

	  optdata::Optical_Category_t category = optdata::kFEMBeamLowGain;
	  if( gain_index == 1 ) category = optdata::kFEMBeamHighGain;
	  
	  optdata::FIFOChannel 
	    beamFIFOChannel( category, 
			     gateWindowTime[gateIndex], 
			     gateFrame[gateIndex],
			     channel,
			     beamEndBin[gateIndex] - beamBeginBin[gateIndex]);

	  mf::LogDebug("OpticalFEM")
	  //std::cout
	    << "Writing beam gate FIFO entry: chanel=" << channel
	    << " at frame=" << gateFrame[gateIndex]
	    << " slice="    << gateWindowTime[gateIndex]
	    << " beginBin=" << beamBeginBin[gateIndex]
	    << " endBin="   << beamEndBin[gateIndex] << std::endl;
	
	  // Copy the time slices.
	  for ( optdata::TimeSlice_t t = beamBeginBin[gateIndex]; 
		t != beamEndBin[gateIndex]; ++t )
	    beamFIFOChannel.push_back( channelData[t] );
	
	  // Dump the FIFO channels as histograms. 
	  if (fm_hist) {
	    art::ServiceHandle<art::TFileService> tfs;
	    std::ostringstream hname;
	    hname << "BFIFO_R" << event.run()
		  << "E" << event.id().event()
		  << "G" << gain_index
		  << "C" << channel
		  << "F" << beamFIFOChannel.Frame()
		  << "S" << beamFIFOChannel.TimeSlice();
	    std::ostringstream htitle;
	    htitle << ";Beam FIFO ADC counts for Run " << event.run()
		   << " Event " << event.id().event()
		   << " Gain " << gain_index 
		   << " Channel " << channel 
		   << " Frame " << beamFIFOChannel.Frame()
		   << " Sample " << beamFIFOChannel.TimeSlice()
		   << ";";
	    TH1* fifoHist = tfs->make<TH1S>(hname.str().c_str(),
					    htitle.str().c_str(),
					    beamFIFOChannel.size(), 
					    0, beamFIFOChannel.size() );
	    // Reminder: The first bin in a histogram is bin 1, NOT bin 0!
	    for ( size_t i = 0; i != beamFIFOChannel.size(); ++i )
	      fifoHist->SetBinContent(i+1,beamFIFOChannel[i]);
	    // The DIFF isn't written at all, but it's fun to look at. 
	    std::ostringstream dname;
	    dname << "BDIFF_R" << event.run()
		  << "E" << event.id().event()
		  << "G" << gain_index
		  << "C" << channel
		  << "F" << beamFIFOChannel.Frame()
		  << "S" << beamFIFOChannel.TimeSlice();
	    std::ostringstream dtitle;
	    dtitle << ";Beam FIFO DIFF for Run " << event.run()
		   << " Event " << event.id().event()
		   << " Gain " << gain_index 
		   << " Channel " << channel 
		   << " Frame " << beamFIFOChannel.Frame()
		   << " Sample " << beamFIFOChannel.TimeSlice()
		   << ";";
	    TH1* diffHist = tfs->make<TH1S>(dname.str().c_str(),
					    dtitle.str().c_str(),
					    beamFIFOChannel.size(), 
					    0, beamFIFOChannel.size() );
	    for ( size_t i = beamBeginBin[gateIndex], b = 1; 
		  i != beamEndBin[gateIndex]; ++i, ++b ) {
	      if ( i >= fm_delay0[gain_index] )
		diffHist->
		  SetBinContent(b,std::max(0,(int)channelData[i] 
					   - (int)channelData[i - fm_delay0[gain_index]]));
	    } // fill diffHist
	  } // if fm_hist
	
	  // Include this beam-gate channel in the output.
	  channelCollection->push_back( std::move(beamFIFOChannel) );
		      
	} // for each input channel
      } // for each beam gate.
    } // if number of beam gates > 0

    // Save on typing: define a type for the size of the
    // "diff vector": the ADC channel subtracted from
    // itself with a time delay.
    typedef optdata::ChannelData::size_type diffSize_t;
    
    // For now, the diff vector is going to be the same
    // length as the length of the ADC channel. (As we tweak
    // the algorithm, this might change.)
    diffSize_t diffSize = (*channelDataHandle).at(0).size();
    diffVector.resize( diffSize );
    
    // Discriminator processing. Go through each channel in this group.
    for ( auto const& channelData : *channelDataHandle ) {
      // Get readout channel number
      ::optdata::Channel_t channel = channelData.ChannelNumber();
      // Get FEM hardware info
      unsigned int crate, slot, femch;
      ch_map->GetCrateSlotFEMChFromReadoutChannel(channel,crate,slot,femch);
      // Get gain type
      ::opdet::UBOpticalChannelType_t gain_type = ch_map->GetChannelType(channel);
      bool is_logic_channel = false;
      if(gain_type == ::opdet::LogicChannel) {
	// Logic channel also belong to a specific FEM w/ gain setting.
	// Our simulation code currently organizes, however, logic channel type
	// as in the same grouping structure as high & low gain (i.e. cannot be
	// a logic channel AND high/low gain type). Here we hack to figure out the
	// gain type. --Kazu July 15 2015
	unsigned int gain_ref_channel = ch_map->GetChannelNumberFromCrateSlotFEMCh(crate,slot,0);
	gain_type = ch_map->GetChannelType( gain_ref_channel );
	is_logic_channel = true;
      }
      size_t gain_index = 0;
      // Determine the gain category.
      switch( ::opdet::UBOpticalChannelType_t(gain_type) ) {
      case ::opdet::LowGain:
	gain_index=0; break;
      case ::opdet::HighGain:
	gain_index=1; break;
      case ::opdet::LogicChannel:
      default:
	// If reaching this point, we do not understand the cause and it is an error.
	mf::LogError("OpticalFEM") 
	  << "Could not find gain setting for a channel: " << channel
	  << std::endl
	  << "Corresponding hardware address: Crate = " << crate << " Slot = " << slot << " FEMCh = " << femch
	  << std::endl;
	throw cet::exception("OpticalFEM");	
      }
    
      // Fill information for PMT Trigger generation here if relevant
      bool slot_for_trigger = false;
      for(auto const& trigger_slot : fm_triggerFEMSlot) {
	if(slot == (unsigned int)(trigger_slot)) {
	  slot_for_trigger = true;
	  break;
	}
      }
      
      if(slot_for_trigger && slot >= maxADCSum1.size()) {
	maxADCSum1.resize(slot+1);
	maxADCSum3.resize(slot+1);
	multiplicitySum1.resize(slot+1);
	multiplicitySum3.resize(slot+1);

	maxADCSum1[slot].resize(diffVector.size(),0);
	maxADCSum3[slot].resize(diffVector.size(),0);
	multiplicitySum1[slot].resize(diffVector.size(),0);
	multiplicitySum3[slot].resize(diffVector.size(),0);
	
	slot_to_gain_index.resize(slot+1,2);
	slot_to_gain_index[slot] = gain_index;
      }
      
      // The lists of previous discriminiator firings
      // and ADC maximums for this channel.
      fireList_t& fire0 = fireMap[channel][0];
      fireList_t& fire1 = fireMap[channel][1];
      fireList_t& fire3 = fireMap[channel][3];
      maxADCList_t& maxADC1 = maxADCMap[channel][1];
      maxADCList_t& maxADC3 = maxADCMap[channel][3];
      
      // Subtract the channel from a delayed version of itself
      // to get the diff vector. This eliminates pedestals, and
      // will be used for the discriminator and PMT trigger tests.
      for ( diffSize_t i = 0; i < diffSize; ++i ) {
	// For the first "fm_delay0[gain]" slices, we don't have
	// a delayed signal.
	if ( i < fm_delay0[gain_index] ) diffVector[i] = 0;
	else
	  // Make sure that negative results would get chopped
	  // to zero, and that unsigned ints won't "wrap around"
	  // to large positive values during the subtraction.
	  diffVector[i] 
	    = (optdata::ADC_Count_t)std::max(0,(int)channelData[ i ] 
					     - (int)channelData[ i - fm_delay0[gain_index] ]);
      }
      
      // Dump the channels as histograms. 
      if (fm_hist) {
	art::ServiceHandle<art::TFileService> tfs;
	std::ostringstream hname;
	hname << "AR" << event.run()
	      << "E" << event.id().event()
	      << "G" << gain_index
	      << "C" << channel;
	std::ostringstream htitle;
	htitle << ";ADC Counts for Run " << event.run()
	       << " Event " << event.id().event()
	       << " Gain " << gain_index 
	       << " Channel " << channel << ";";		      
	TH1* chanHist = tfs->make<TH1S>(hname.str().c_str(),
					htitle.str().c_str(),
					channelData.size(), 
					0, channelData.size() );
	for ( diffSize_t i = 0; i != channelData.size(); ++i )
	  chanHist->SetBinContent(i+1,channelData[i]);
	std::ostringstream dname;
	dname << "DR" << event.run()
		<< "E" << event.id().event()
		<< "G" << gain_index
		<< "C" << channel;
	  std::ostringstream dtitle;
	  dtitle << ";DIFF for Run " << event.run()
		 << " Event " << event.id().event()
		 << " Gain " << gain_index 
		 << " Channel " << channel << ";";  
	  TH1* diffHist = tfs->make<TH1S>(dname.str().c_str(),
					  dtitle.str().c_str(),
					  diffVector.size(), 
					  0, diffVector.size() );
	  for ( diffSize_t i = 0; i != diffVector.size(); ++i )
	    diffHist->SetBinContent(i+1,diffVector[i]);
      } // if fm_hist
      
      // Scan through the diff vector, testing for the
      // criteria for the discriminators to fire and
      // accumulating trigger sums.  We have to ignore the
      // first fm_delay0[gain] slices, since they haven't been
      // subtracted.
      for ( diffSize_t slice = fm_delay0[gain_index]; slice < diffSize; ++slice ) {
	
	// Check if we're outside all the beam gates AND this is not a logic pulse channel
	bool outsideBeamGates = true;
	//if( channel < kLogicStartChannel ) {
	if ( ch_map->GetChannelType( channel )!=::opdet::LogicChannel && numberOfGates>0 ) {
	  for ( size_t b = 0; b != numberOfGates; ++b )
	    if ( slice >= _beamBeginBin.at(gain_index)[b] && slice < _beamEndBin.at(gain_index)[b] )
	      outsideBeamGates = false;
	}
	
	// Check if discriminator 0 fired.
	if ( diffVector[slice] >= fm_threshold0[gain_index] )
	  {
	    // Did the previous discriminators (if any) fire too
	    // soon before this one?
	    if ( ( fire0.empty()  ||
		   fire0.back() + fm_disc0quiet[gain_index] < slice ) &&
		 ( fire1.empty()  ||
		   fire1.back() + fm_disc1deadtime[gain_index] < slice ) &&
		 ( fire3.empty()  ||
		   fire3.back() + fm_disc3deadtime[gain_index] < slice ) 
		 ) {
	      // No, so discriminator 0 fires. Add this
	      // slice to the list of slices for disc 0.
	      fire0.push_back( slice );
	      
	      LOG_DEBUG("OpticalFEM")
		<< "Disc 0 fires, channel=" << channel
		<< " at frame=" 
		<< (*channelDataHandle).Frame() + slice / clock.FrameTicks()
		<< " slice=" << slice % clock.FrameTicks();
	    }
	  } // threshold0 satisfied
	
	// Check if we're outside the beam gate and we've over the
	// discriminator 1 threshold.
	if ( outsideBeamGates && diffVector[slice] >= fm_threshold1[gain_index] ) {
	  
	  // See if the most recent firing of disc 0 occurred
	  // recently enough to fall within the disc 1 window, but
	  // outside the dead time from the last disc 1.
	  if ( ! fire0.empty()                              &&
	       slice - fire0.back() < fm_disc0window[gain_index]  &&
	       ( fire1.empty()  || 
		 fire1.back() + fm_disc1deadtime[gain_index] < slice ) ) {
	    
	    // Discriminator 1 fires. Add this slice to the list of
	    // slices for disc 1.
	    fire1.push_back( slice );
	    
	    // Look ahead in the diff vector to find the maximum ADC
	    // count within the discriminator width.
	    optdata::ADC_Count_t maxADC = 0;
	    optdata::TimeSlice_t endWidth 
	      = std::min( diffSize, slice + fm_disc1width[gain_index] );
	    for ( optdata::TimeSlice_t s = slice;
		  s != endWidth; ++s )
	      maxADC = std::max( maxADC, diffVector[s] );
	    
	    // Save this value for PMT trigger tests IF not logic pulse channel
	    if( !is_logic_channel ) maxADC1.push_back( maxADC );
	    
	    // Go back (if negative) or forward (if positive) from
	    // the point of the last disc0 firing to start saving
	    // slices.
	    optdata::TimeSlice_t saveSlice = slice + fm_delay1[gain_index];
	    // Make sure we don't go "off the end" of our data.
	    if ( (int)slice + (int)fm_delay1[gain_index] < 0 ) saveSlice = 0;
	    if ( saveSlice + fm_cosmicSlices[gain_index] > diffSize )
	      saveSlice = diffSize - fm_cosmicSlices[gain_index];
	    
	    // Time information for this FIFO channel.
	    
	    optdata::Frame_t cosmicFrame 
	      = firstFrame + (firstSlice + saveSlice) / clock.FrameTicks();
	    optdata::TimeSlice_t cosmicTime 
	      = (firstSlice + saveSlice) % clock.FrameTicks();
	    
	    LOG_DEBUG("OpticalFEM")
	      << "Disc 1 fires, Writing cosmic channel=" << channel
	      << " at frame=" << cosmicFrame
	      << " slice=" << cosmicTime
	      << " begin=" << saveSlice
	      << " end=" << saveSlice+fm_cosmicSlices[gain_index]
	      << " max ADC=" << maxADC; 

	    // Create a new FIFO channel, copying the channel
	    // number from the input: wrong categories.
	    
	    optdata::Optical_Category_t category = optdata::kFEMCosmicLowGain;
	    if ( gain_index == 1 ) category = optdata::kFEMCosmicHighGain;
	    else if ( gain_index == 2 ) category = optdata::kFEMCosmicLogicPulse;
	    
	    optdata::FIFOChannel
	      cosmicChannel( category, 
			     cosmicTime, 
			     cosmicFrame,
			     channel,
			     fm_cosmicSlices[gain_index] );
	    
	    // Copy the time slices.
	    for ( optdata::TimeSlice_t t = saveSlice; 
		  t != saveSlice + fm_cosmicSlices[gain_index]; ++t )
	      cosmicChannel.push_back( channelData[t] );
	    
	    if (fm_hist) {
	      // Dump the FIFO channels as histograms. 
	      art::ServiceHandle<art::TFileService> tfs;
	      std::ostringstream hname;
	      hname << "CFIFO_R" << event.run()
		    << "E" << event.id().event()
		    << "G" << gain_index
		    << "C" << channel
		    << "F" << cosmicChannel.Frame()
		    << "S" << cosmicChannel.TimeSlice();
	      std::ostringstream htitle;
	      htitle << ";Cosmic FIFO ADC counts for Run " << event.run()
		     << " Event " << event.id().event()
		     << " Gain " << gain_index
		     << " Channel " << channel 
		     << " Frame " << cosmicChannel.Frame()
		     << " Sample " << cosmicChannel.TimeSlice()
		     << ";";
	      TH1* fifoHist = tfs->make<TH1S>(hname.str().c_str(),
					      htitle.str().c_str(),
					      cosmicChannel.size(), 
					      0, cosmicChannel.size() );
	      for ( size_t i = 0; i != cosmicChannel.size(); ++i )
		fifoHist->SetBinContent(i+1,cosmicChannel[i]);
	      // The DIFF vector is not actually output, but it's fun to look at. 
	      std::ostringstream dname;
	      dname << "CDIFF_R" << event.run()
		    << "E" << event.id().event()
		    << "G" << gain_index
		    << "C" << channel
		    << "F" << cosmicChannel.Frame()
		    << "S" << cosmicChannel.TimeSlice();
	      std::ostringstream dtitle;
	      dtitle << ";Cosmic FIFO DIFF for Run " << event.run()
		     << " Event " << event.id().event()
		     << " Gain " << gain_index 
		     << " Channel " << channel 
		     << " Frame " << cosmicChannel.Frame()
		     << " Sample " << cosmicChannel.TimeSlice()
		     << ";";
	      TH1* diffHist = tfs->make<TH1S>(dname.str().c_str(),
					      dtitle.str().c_str(),
					      cosmicChannel.size(), 
					      0, cosmicChannel.size() );
	      for ( size_t i = saveSlice, b = 1; 
		    i != saveSlice + fm_cosmicSlices[gain_index]; ++i, ++b )
		diffHist->SetBinContent(b,diffVector[i]);
	    } // if fm_hist
	    
	    // Include this beam-gate channel in the output.
	    channelCollection->push_back( std::move(cosmicChannel) );
	    
	  } // disc 1 fired
	} // outside beam gate and threshold1 satisfied
	
	// If we're inside the beam gate, check if we've crossed
	// the discriminator 3 threshold.
	if ( (! outsideBeamGates)  &&  diffVector[slice] >= fm_threshold3[gain_index] ) {
	  // See if the most recent firing of disc 0 occurred recently
	  // enough to fall within the disc 3 window, but outside the
	  // dead time from the last disc 3.
	  if ( ! fire0.empty()                              &&
	       slice - fire0.back() < fm_disc0window[gain_index]  &&
	       ( fire3.empty()  || 
		 fire3.back() + fm_disc3deadtime[gain_index] < slice ) ) {
	    // Discriminator 3 fires. Add this slice to the list
	    // of slices for disc 3.
	    fire3.push_back( slice );
	    
	    // Look ahead in the diff vector to find the maximum
	    // ADC count within the discriminator width.
	    optdata::ADC_Count_t maxADC = 0;
	    optdata::TimeSlice_t endWidth 
	      = std::min( diffSize, slice + fm_disc3width[gain_index] );
	    for ( optdata::TimeSlice_t s = slice;
		  s != endWidth; ++s )
	      maxADC = std::max( maxADC, diffVector[s] );
	    
	    // Save this value for PMT trigger tests if not logic channel.
	    if( !is_logic_channel ) maxADC3.push_back( maxADC );
	    
	    LOG_DEBUG("OpticalFEM")
	      << "Disc 3 fires, channel=" << channel
	      << " at frame=" 
	      << (*channelDataHandle).Frame() + slice / clock.FrameTicks()
	      << " slice=" << slice % clock.FrameTicks()
	      << " max ADC=" << maxADC;

	    
	  } // disc 3 fired
	} // inside beam gate and threshold3 satisfied

	if(slot_for_trigger && !is_logic_channel) {
	  // Accumulate the PMT trigger sums. For each of
	  // discriminator {1,3}, see if we're within the width of
	  // that discriminator since the last time it fired. If so,
	  // sum the max ADC count and add to the multiplicity count.
	  // Only if this channel is NOT logic pulse channel
	  if ( ! fire1.empty ()   &&
	       slice < fire1.back() + fm_disc1deadtime[gain_index] ) {
	    maxADCSum1[slot][slice] += maxADC1.back();
	    ++multiplicitySum1[slot][slice];
	  }
	  if ( ! fire3.empty ()   &&
	       slice < fire3.back() + fm_disc3deadtime[gain_index] ) {
	    maxADCSum3[slot][slice] += maxADC3.back();
	    ++multiplicitySum3[slot][slice];
	  }
	}	
      } // diffVector time slice
    } // for each channel
    
    // Dump the trigger sums as histograms. 
    if (fm_hist) {
      for(size_t slot=0; slot<maxADCSum1.size(); ++slot) {
	if(!maxADCSum1[slot].size()) continue;
	art::ServiceHandle<art::TFileService> tfs;
	std::ostringstream h1name;
	h1name << "ADC1_R" << event.run()
	       << "E" << event.id().event()
	       << "G" << 0;
	std::ostringstream h1title;
	h1title << ";Sum Max ADC counts for Run " << event.run()
		<< " Event " << event.id().event()
		<< " Gain " << 0
		<< " Discr 1"
		<< ";";
	TH1* adc1Hist = tfs->make<TH1S>(h1name.str().c_str(),
					h1title.str().c_str(),
					maxADCSum1[slot].size(), 
					0, maxADCSum1[slot].size() );
	for ( size_t i = 0; i != maxADCSum1[slot].size(); ++i )
	  adc1Hist->SetBinContent(i+1,maxADCSum1[slot][i]);
	std::ostringstream m1name;
	m1name << "Mult1_R" << event.run()
	       << "E" << event.id().event()
	       << "G" << 0;
	std::ostringstream m1title;
	m1title << ";Sum Multiplicities for Run " << event.run()
		<< " Event " << event.id().event()
		<< " Gain " << 0	  
		<< " Discr 1"
		<< ";";
	TH1* mul1Hist = tfs->make<TH1S>(m1name.str().c_str(),
					m1title.str().c_str(),
					multiplicitySum1[slot].size(), 
				      0, multiplicitySum1[slot].size() );
      for ( size_t i = 0; i != multiplicitySum1[slot].size(); ++i )
	mul1Hist->SetBinContent(i+1,multiplicitySum1[slot][i]);
      std::ostringstream h3name;
      h3name << "ADC3_R" << event.run()
	     << "E" << event.id().event()
	     << "G" << 1;
      std::ostringstream h3title;
      h3title << ";Sum Max ADC counts for Run " << event.run()
	      << " Event " << event.id().event()
	      << " Gain " << 1
	      << " Discr 3"
	      << ";";
      TH1* adc3Hist = tfs->make<TH1S>(h3name.str().c_str(),
				      h3title.str().c_str(),
				      maxADCSum3[slot].size(), 
				      0, maxADCSum3[slot].size() );
      for ( size_t i = 0; i != maxADCSum3[slot].size(); ++i )
	adc3Hist->SetBinContent(i+1,maxADCSum3[slot][i]);
      std::ostringstream m3name;
      m3name << "Mult3_R" << event.run()
	     << "E" << event.id().event()
	     << "G" << 1;
      std::ostringstream m3title;
      m3title << ";Sum Multiplicities for Run " << event.run()
	      << " Event " << event.id().event()
	      << " Gain " << 1
	      << " Discr 3"
	      << ";";
      TH1* mul3Hist = tfs->make<TH1S>(m3name.str().c_str(),
				      m3title.str().c_str(),
				      multiplicitySum3[slot].size(), 
				      0, multiplicitySum3[slot].size() );
      for ( size_t i = 0; i != multiplicitySum3[slot].size(); ++i )
	mul3Hist->SetBinContent(i+1,multiplicitySum3[slot][i]);
      }
    } // if fm_hist
    
    // std::cout <<  "PMT Trigger processing. We only do this for the high-gain FEM" << std::endl;

    for(auto const& slot : fm_triggerFEMSlot) {
      auto gain_index = slot_to_gain_index[slot];
      // Go through the time slices again, testing if the PMT trigger
      // fires. 
      for ( diffSize_t slice = 0; slice < diffSize; ++slice ) {

	// If we're still within the deadtime of the last
	// trigger, skip the trigger tests.
	if ( ! triggers.empty()   &&
	     slice < triggers.back().second + fm_triggerDeadtime[gain_index] )
	  continue;
	
	if ( maxADCSum1[slot][slice]       > fm_cosmicThreshold[gain_index]  &&
	     multiplicitySum1[slot][slice] > fm_cosmicMultiplicity[gain_index] ) {
	  // Create a trigger record and save it. 
	  optdata::Frame_t frame 
	    = firstFrame + (firstSlice + slice) / clock.FrameTicks();
	  optdata::TimeSlice_t sample 
	    = (firstSlice + slice) % clock.FrameTicks();
	  optdata::PMTTrigger cosmicTrigger( optdata::kCosmicPMTTrigger,
					     sample, frame );
	  triggerCollection->push_back( std::move(cosmicTrigger) );
	  
	  // Save the trigger in our list.
	  triggers.push_back( trigger_t( optdata::kCosmicPMTTrigger, slice ) );
	  
	  LOG_DEBUG("OpticalFEM")
	    << "Cosmic PMT Trigger"
	    << " at frame=" << frame
	    << " slice=" << sample
	    << " max ADC sum=" << maxADCSum1[slot][slice]
	    << " multiplicity sum=" << multiplicitySum1[slot][slice]
	    << std::endl;

	  if (fm_hist) {
	    // Dump some trigger sums as histograms. 
	    // Pick some arbitrary display range for the trigger info.
	    size_t hBegin = slice + fm_delay1[gain_index];
	    size_t hEnd = slice + fm_disc1deadtime[gain_index];
	    size_t hSize = fm_disc1deadtime[gain_index] - fm_delay1[gain_index];
	    art::ServiceHandle<art::TFileService> tfs;
	    std::ostringstream hname;
	    hname << "CADC1_R" << event.run()
		  << "E" << event.id().event()
		  << "G" << gain_index
		  << "S" << slice;
	    std::ostringstream htitle;
	    htitle << ";Cosmic Trigger Sum Max ADC counts for Run " << event.run()
		   << " Event " << event.id().event()
		   << " Gain " << gain_index 
		   << " Sample " << slice
		   << " Discr 1;";
	    TH1* adcHist = tfs->make<TH1S>(hname.str().c_str(),
					   htitle.str().c_str(),
					   hSize, hBegin, hEnd );
	    for ( size_t i = hBegin, b = 1; i != hEnd; ++i, ++b )
	      adcHist->SetBinContent(b,maxADCSum1[slot][i]);
	    std::ostringstream mname;
	    mname << "CMult1_R" << event.run()
		  << "E" << event.id().event()
		  << "G" << gain_index
		  << "S" << slice;
	    std::ostringstream mtitle;
	    mtitle << ";Cosmic Trigger Sum Multiplicities for Run " << event.run()
		   << " Event " << event.id().event()
		   << " Gain " << gain_index 
		   << " Sample " << slice
		   << " Discr 1;";
	    TH1* mulHist = tfs->make<TH1S>(mname.str().c_str(),
					   mtitle.str().c_str(),
					   hSize,hBegin, hEnd );
	    for ( size_t i = hBegin, b = 1; i != hEnd; ++i, ++b )
	      mulHist->SetBinContent(b,multiplicitySum1[slot][i]);
	  } // if fm_hist
	  
	} // Cosmic PMT trigger
	      
	if ( maxADCSum3[slot][slice]       > fm_beamThreshold[gain_index]  &&
	     multiplicitySum3[slot][slice] > fm_beamMultiplicity[gain_index] ) {
	  // Create a trigger record and save it. 
	  optdata::Frame_t frame 
	    = firstFrame + (firstSlice + slice) / clock.FrameTicks();
	  optdata::TimeSlice_t sample 
	    = (firstSlice + slice) % clock.FrameTicks();
	  optdata::PMTTrigger beamTrigger( optdata::kBeamPMTTrigger,
					   sample, frame );
	  triggerCollection->push_back( std::move(beamTrigger) );
	  
	  // Save the trigger in our list.
	  triggers.push_back( trigger_t( optdata::kCosmicPMTTrigger, slice ) );
	  
	  LOG_DEBUG("OpticalFEM")
	    << "Beam PMT Trigger"
	    << " at frame=" << frame
	    << " slice=" << sample
	    << " max ADC sum=" << maxADCSum3[slot][slice]
	    << " multiplicity sum=" << multiplicitySum3[slot][slice];

	  if (fm_hist) {
	    // Dump some trigger sums as histograms. 
	    // Pick some arbitrary display range for the trigger info.
	    size_t hBegin = slice + fm_delay1[gain_index];
	    size_t hEnd = slice + fm_disc3deadtime[gain_index];
	    size_t hSize = fm_disc1deadtime[gain_index] - fm_delay1[gain_index];
	    art::ServiceHandle<art::TFileService> tfs;
	    std::ostringstream hname;
	    hname << "BADC3_R" << event.run()
		  << "E" << event.id().event()
		  << "G" << gain_index
		  << "S" << slice;
	    std::ostringstream htitle;
	    htitle << ";Beam Trigger Sum Max ADC counts for Run " << event.run()
		   << " Event " << event.id().event()
		   << " Gain " << gain_index 
		   << " Sample " << slice
		   << " Discr 3;";
	    TH1* adcHist = tfs->make<TH1S>(hname.str().c_str(),
					   htitle.str().c_str(),
					   hSize, hBegin, hEnd );
	    for ( size_t i = hBegin, b = 1; i != hEnd; ++i, ++b )
	      adcHist->SetBinContent(b,maxADCSum3[slot][i]);
	    std::ostringstream mname;
	    mname << "BMult3_R" << event.run()
		  << "E" << event.id().event()
		  << "G" << gain_index
		  << "S" << slice;
	    std::ostringstream mtitle;
	    mtitle << ";Beam Trigger Sum Multiplicities for Run " << event.run()
		   << " Event " << event.id().event()
		   << " Gain " << gain_index 
		   << " Sample " << slice
		   << " Discr 3;";
	    TH1* mulHist = tfs->make<TH1S>(mname.str().c_str(),
					   mtitle.str().c_str(),
					   hSize,hBegin, hEnd );
	    for ( size_t i = hBegin, b = 1; i != hEnd; ++i, ++b )
	      mulHist->SetBinContent(b,multiplicitySum3[slot][i]);
	  } // if fm_hist
	  
	} // Beam PMT Trigger
      } // for each time slice
    } // trigger processing for high-gain channels

    // Write out all the channels and triggers,
    event.put( std::move( channelCollection ) );
    event.put( std::move( triggerCollection ) );    
  }


} // namespace opdet
