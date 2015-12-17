////////////////////////////////////////////////////////////////////////
// Class:       OverlayRawDataDetailMicroBooNE
// Module Type: producer
// File:        OverlayRawDataDetailMicroBooNE_module.cc
//
// This borrows a lot from the Mu2e mixing module:
//      EventMixing/src/MixMCEvents_module.cc
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Modules/MixFilter.h"
#include "art/Framework/IO/ProductMix/MixHelper.h"
#include "art/Framework/Core/PtrRemapper.h"
#include "art/Persistency/Common/CollectionUtilities.h"
#include "art/Utilities/InputTag.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <exception>
#include <sstream>

#include "CalibrationDBI/Interface/IChannelStatusService.h"
#include "CalibrationDBI/Interface/IChannelStatusProvider.h"

#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"

#include "Simulation/SimChannel.h"
#include "Simulation/SimPhotons.h"
#include "Simulation/AuxDetSimChannel.h"

#include "MCBase/MCTrack.h"
#include "MCBase/MCShower.h"

#include "DataOverlay/RawDigitMixer.h"
#include "RawData/RawDigit.h"
#include "SimpleTypesAndConstants/RawTypes.h"

#include "DataOverlay/OpDetWaveformMixer.h"
#include "RawData/OpDetWaveform.h"

#include "DataOverlayProducts/EventMixingSummary.h"

namespace mix {
  class OverlayRawDataDetailMicroBooNE;
  typedef art::MixFilter<OverlayRawDataDetailMicroBooNE> OverlayRawDataMicroBooNE;
}

class mix::OverlayRawDataDetailMicroBooNE : public boost::noncopyable {
public:

  OverlayRawDataDetailMicroBooNE(fhicl::ParameterSet const& p,
				 art::MixHelper &helper);

  void startEvent(const art::Event&);  //called at the start of every event
  void finalizeEvent(art::Event &);    //called at the end of every event
  
  size_t nSecondaries() { return fEventsToMix; } 

  void processEventIDs(art::EventIDSequence const& seq); //bookkepping for event IDs

  // Mixing Functions

  // For now, allow exactly one  input. Assume MC inputs have been merged
  // previously and one detsim output created if needed. This could be changed
  // but would require mixing functions for MC here.

  //a lot of MC collections are just simple copies of the collections...
  template<typename T>
  bool MixSimpleCopy( std::vector< std::vector<T> const*> const& inputs,
		      std::vector< T > & output,
		      art::PtrRemapper const &);
  
  bool MixRawDigits( std::vector< std::vector<raw::RawDigit> const* > const& inputs,
  		     std::vector<raw::RawDigit> & output,
  		     art::PtrRemapper const &);
  
  bool MixOpDetWaveforms_HighGain( std::vector< std::vector<raw::OpDetWaveform> const* > const& inputs,
				   std::vector<raw::OpDetWaveform> & output,
				   art::PtrRemapper const &);
  bool MixOpDetWaveforms_LowGain( std::vector< std::vector<raw::OpDetWaveform> const* > const& inputs,
				  std::vector<raw::OpDetWaveform> & output,
				  art::PtrRemapper const &);

		   
  
private:

  // Declare member data here.
  RawDigitMixer              fRDMixer;
  OpDetWaveformMixer         fODMixer;
  
  fhicl::ParameterSet  fpset;
  short                fDefaultRawDigitSatPoint;
  short                fDefaultOpDetSatPoint;
  size_t               fOpDetMinSampleSize;
  bool                 fInputFileIsData;

  std::string          fRawDigitDataModuleLabel;
  std::string          fOpDetDataModuleLabel;
  std::string          fRawDigitMCModuleLabel;
  std::string          fOpDetMCModuleLabel;

  std::string          fRawDigitInputSourceModuleLabel;
  std::string          fOpDetInputSourceModuleLabel;
  std::string          fRawDigitMixerSourceModuleLabel;
  std::string          fOpDetMixerSourceModuleLabel;

  std::string          fG4InputModuleLabel;
  std::string          fGeneratorInputModuleLabel;

  bool                 fDoMCReco;
  std::string          fMCRecoInputModuleLabel;

  size_t               fEventsToMix;
  float                fDefaultMCRawDigitScale;
  float                fDefaultMCOpDetScale;
  
  art::Handle< std::vector<raw::RawDigit> > inputDigitHandle;
  art::Handle< std::vector<raw::OpDetWaveform> > inputOpDetHandle_HighGain;
  art::Handle< std::vector<raw::OpDetWaveform> > inputOpDetHandle_LowGain;

  void GenerateMCRawDigitScaleMap(std::vector<raw::RawDigit> const&);
  std::unordered_map<raw::ChannelID_t,float> fMCRawDigitScaleMap;

  void GenerateMCOpDetHighGainScaleMap(std::vector<raw::OpDetWaveform> const&);
  std::unordered_map<raw::Channel_t,float> fMCOpDetHighGainScaleMap;

  void GenerateMCOpDetLowGainScaleMap(std::vector<raw::OpDetWaveform> const&);
  std::unordered_map<raw::Channel_t,float> fMCOpDetLowGainScaleMap;

  std::unique_ptr< std::vector<mix::EventMixingSummary> > fEventMixingSummary;
  
};


mix::OverlayRawDataDetailMicroBooNE::OverlayRawDataDetailMicroBooNE(fhicl::ParameterSet const& p,
								    art::MixHelper &helper)
  :
  fRDMixer(false), //print warnings turned off
  fODMixer(false), //print warnings turned off
  fpset(p.get<fhicl::ParameterSet>("detail")),
  fDefaultRawDigitSatPoint(fpset.get<short>("DefaultRawDigitSaturationPoint",4096)),
  fDefaultOpDetSatPoint(fpset.get<short>("DefaultOpDetSaturationPoint",4096)),
  fOpDetMinSampleSize(fpset.get<size_t>("OpDetMinSampleSize",100)),
  fInputFileIsData(fpset.get<bool>("InputFileIsData")),
  fRawDigitDataModuleLabel(fpset.get<std::string>("RawDigitDataModuleLabel")),
  fOpDetDataModuleLabel(fpset.get<std::string>("OpDetMCModuleLabel")),
  fRawDigitMCModuleLabel(fpset.get<std::string>("RawDigitDataModuleLabel")),
  fOpDetMCModuleLabel(fpset.get<std::string>("OpDetMCModuleLabel")),
  fEventsToMix(fpset.get<size_t>("EventsToMix",1)),
  fDefaultMCRawDigitScale(fpset.get<float>("DefaultMCRawDigitScale",1)),
  fDefaultMCOpDetScale(fpset.get<float>("DefaultMCOpDetScale",1)),
  fEventMixingSummary(nullptr)
{
  
  if(fEventsToMix!=1){
    std::stringstream err_str;
    err_str << "ERROR! Really sorry, but we can only do mixing for one collection right now! ";
    err_str << "\nYep. We're gonna throw an exception now. You should change your fcl to set 'EventsToMix' to 1";
    throw cet::exception("OverlayRawDataMicroBooNE") << err_str.str() << std::endl;;
  }

  if(fInputFileIsData){
    fRawDigitInputSourceModuleLabel = fRawDigitDataModuleLabel;
    fOpDetInputSourceModuleLabel    = fOpDetDataModuleLabel;
    fRawDigitMixerSourceModuleLabel = fRawDigitMCModuleLabel;
    fOpDetMixerSourceModuleLabel    = fOpDetMCModuleLabel;
  }
  else if(!fInputFileIsData){
    fRawDigitInputSourceModuleLabel = fRawDigitMCModuleLabel;
    fOpDetInputSourceModuleLabel    = fOpDetMCModuleLabel;
    fRawDigitMixerSourceModuleLabel = fRawDigitDataModuleLabel;
    fOpDetMixerSourceModuleLabel    = fOpDetDataModuleLabel;
  }
  
  if(fInputFileIsData){
    fDoMCReco = fpset.get_if_present<std::string>("MCRecoInputModuleLabel",fMCRecoInputModuleLabel);
    fG4InputModuleLabel = fpset.get<std::string>("G4InputModuleLabel");
    fGeneratorInputModuleLabel = fpset.get<std::string>("GeneratorInputModuleLabel");
    
    //MC generator info is a simple copy
    helper.declareMixOp( art::InputTag(fGeneratorInputModuleLabel),
			 &OverlayRawDataDetailMicroBooNE::MixSimpleCopy<simb::MCTruth>,
			 *this );
    
    //Simple copies of G4 SimPhotons, MCParticles, SimChannels, and SimAuxDetChannel
    helper.declareMixOp( art::InputTag(fG4InputModuleLabel),
			 &OverlayRawDataDetailMicroBooNE::MixSimpleCopy<simb::MCParticle>,
			 *this );
    helper.declareMixOp( art::InputTag(fG4InputModuleLabel),
			 &OverlayRawDataDetailMicroBooNE::MixSimpleCopy<sim::SimPhotons>,
			 *this );
    helper.declareMixOp( art::InputTag(fG4InputModuleLabel),
			 &OverlayRawDataDetailMicroBooNE::MixSimpleCopy<sim::SimChannel>,
			 *this );
    helper.declareMixOp( art::InputTag(fG4InputModuleLabel),
			 &OverlayRawDataDetailMicroBooNE::MixSimpleCopy<sim::AuxDetSimChannel>,
			 *this );
    /*
    //Associations of MCParticles to MCTruth...hopefully a simple copy is enough
    helper.declareMixOp( art::InputTag(fG4InputModuleLabel),
    &OverlayRawDataDetailMicroBooNE::MixSimpleCopy
		       < art::Assns<simb::MCTruth,simb::MCParticle,void> >,
		       *this );
		       */
    
    //Copies of MCShower and MCTrack
    if(fDoMCReco){
      helper.declareMixOp( art::InputTag(fMCRecoInputModuleLabel),
			 &OverlayRawDataDetailMicroBooNE::MixSimpleCopy<sim::MCShower>,
			   *this );
      helper.declareMixOp( art::InputTag(fMCRecoInputModuleLabel),
			   &OverlayRawDataDetailMicroBooNE::MixSimpleCopy<sim::MCTrack>,
			   *this );
    }
  }//end if file is input data

  helper.declareMixOp( art::InputTag(fRawDigitMixerSourceModuleLabel),
		       &OverlayRawDataDetailMicroBooNE::MixRawDigits,
		       *this );
  
  helper.declareMixOp( art::InputTag(fOpDetMixerSourceModuleLabel,"OpdetBeamHighGain"),
		       &OverlayRawDataDetailMicroBooNE::MixOpDetWaveforms_HighGain,
		       *this );
  helper.declareMixOp( art::InputTag(fOpDetMixerSourceModuleLabel,"OpdetBeamLowGain"),
		       &OverlayRawDataDetailMicroBooNE::MixOpDetWaveforms_LowGain,
		       *this );

  //If it produces something on its own, declare it here
  helper.produces< std::vector<mix::EventMixingSummary> >();
}

//Initialize for each event
void mix::OverlayRawDataDetailMicroBooNE::startEvent(const art::Event& event) {

  if(!( (event.isRealData() && fInputFileIsData) || (!event.isRealData() && !fInputFileIsData)))
    throw cet::exception("OverlayRawDataMicroBooNE") << "Input file claimed to be data/not data, but it's not." << std::endl;;
    
  
  event.getByLabel(fRawDigitInputSourceModuleLabel,inputDigitHandle);
  if(!inputDigitHandle.isValid())
    throw cet::exception("OverlayRawDataMicroBooNE") << "Bad input digit handle." << std::endl;;
  fRDMixer.SetSaturationPoint(fDefaultRawDigitSatPoint);

  event.getByLabel(fOpDetInputSourceModuleLabel,"OpdetBeamLowGain",inputOpDetHandle_LowGain);
  if(!inputOpDetHandle_LowGain.isValid())
    throw cet::exception("OverlayRawDataMicroBooNE") << "Bad input opdet lowgain handle." << std::endl;;

  event.getByLabel(fOpDetInputSourceModuleLabel,"OpdetBeamHighGain",inputOpDetHandle_HighGain);
  if(!inputOpDetHandle_HighGain.isValid())
    throw cet::exception("OverlayRawDataMicroBooNE") << "Bad input opdet highgain handle." << std::endl;;

  fODMixer.SetSaturationPoint(fDefaultOpDetSatPoint);
  fODMixer.SetMinSampleSize(fOpDetMinSampleSize);


  fEventMixingSummary.reset(new std::vector<mix::EventMixingSummary>);
}

//For each of the mixed in events...bookkepping for event IDs
void mix::OverlayRawDataDetailMicroBooNE::processEventIDs(art::EventIDSequence const& seq){
  for (auto const& id : seq)
    fEventMixingSummary->emplace_back(id.event(),id.subRun(),id.run());
}


//End each event
void mix::OverlayRawDataDetailMicroBooNE::finalizeEvent(art::Event& event) {
  event.put(std::move(fEventMixingSummary));
}

template<typename T>
bool mix::OverlayRawDataDetailMicroBooNE::MixSimpleCopy( std::vector< std::vector<T> const*> const& inputs,
							 std::vector< T > & output,
							 art::PtrRemapper const &){
  art::flattenCollections(inputs,output);
  return true;
}

void mix::OverlayRawDataDetailMicroBooNE::GenerateMCRawDigitScaleMap(std::vector<raw::RawDigit> const& dataDigitVector){
  //right now, assume the number of channels is the number in the collection
  //and, loop through the channels one by one to get the right channel number
  //note: we will put here access to the channel database to determine dead channels
  fMCRawDigitScaleMap.clear();

  const lariov::IChannelStatusProvider& chanStatus = art::ServiceHandle<lariov::IChannelStatusService>()->GetProvider();
  
  for(auto const& d : dataDigitVector){
    if(chanStatus.IsBad(d.Channel()))
      fMCRawDigitScaleMap[d.Channel()] = 0.0;
    else
      fMCRawDigitScaleMap[d.Channel()] = fDefaultMCRawDigitScale;
  }
}

bool mix::OverlayRawDataDetailMicroBooNE::MixRawDigits( std::vector< std::vector<raw::RawDigit> const* > const& inputs,
							std::vector<raw::RawDigit> & output,
							art::PtrRemapper const & remap) {
  
  //make sure we only have two collections for now
  if(inputs.size()!=fEventsToMix || (inputs.size()!=1 && !fInputFileIsData)){
    std::stringstream err_str;
    err_str << "ERROR! We have the wrong number of collections of raw digits we are adding! " << inputs.size();
    throw std::runtime_error(err_str.str());
  }


  if(fInputFileIsData){
    GenerateMCRawDigitScaleMap(*inputDigitHandle);  
    fRDMixer.DeclareData(*inputDigitHandle);
    for(auto const& icol : inputs)
      fRDMixer.Mix(*icol,fMCRawDigitScaleMap);
  }
  else if(!fInputFileIsData){
    GenerateMCRawDigitScaleMap(*(inputs[0]));  
    fRDMixer.DeclareData(*(inputs[0]));
    fRDMixer.Mix(*inputDigitHandle,fMCRawDigitScaleMap);
  }
  
  fRDMixer.FillRawDigitOutput(output);
  
  return true;
}

void mix::OverlayRawDataDetailMicroBooNE::GenerateMCOpDetHighGainScaleMap(std::vector<raw::OpDetWaveform> const& dataVector){
  //right now, assume the number of channels is the number in the collection
  //and, loop through the channels one by one to get the right channel number
  //note: we will put here access to the channel database to determine dead channels
  fMCOpDetHighGainScaleMap.clear();
  for(auto const& d : dataVector)
    fMCOpDetHighGainScaleMap[d.ChannelNumber()] = fDefaultMCOpDetScale;
}

void mix::OverlayRawDataDetailMicroBooNE::GenerateMCOpDetLowGainScaleMap(std::vector<raw::OpDetWaveform> const& dataVector){
  //right now, assume the number of channels is the number in the collection
  //and, loop through the channels one by one to get the right channel number
  //note: we will put here access to the channel database to determine dead channels
  fMCOpDetLowGainScaleMap.clear();
  for(auto const& d : dataVector)
    fMCOpDetLowGainScaleMap[d.ChannelNumber()] = fDefaultMCOpDetScale;
}

bool mix::OverlayRawDataDetailMicroBooNE::MixOpDetWaveforms_HighGain( std::vector< std::vector<raw::OpDetWaveform> const* > const& inputs,
								      std::vector<raw::OpDetWaveform> & output,
								      art::PtrRemapper const & remap) {
  
  //make sure we only have two collections for now
  if(inputs.size()!=fEventsToMix || (inputs.size()!=1 && !fInputFileIsData)){
    std::stringstream err_str;
    err_str << "ERROR! We have the wrong number of collections of raw digits we are adding! " << inputs.size();
    throw std::runtime_error(err_str.str());
  }


  if(fInputFileIsData){
    GenerateMCOpDetHighGainScaleMap(*inputOpDetHandle_HighGain);  
    fODMixer.DeclareData(*inputOpDetHandle_HighGain,output);
    for(auto const& icol : inputs)
      fODMixer.Mix(*icol,fMCOpDetHighGainScaleMap,output);
  }
  else if(!fInputFileIsData){
    GenerateMCOpDetHighGainScaleMap(*(inputs[0]));  
    fODMixer.DeclareData(*(inputs[0]),output);
    fODMixer.Mix(*inputOpDetHandle_HighGain,fMCOpDetHighGainScaleMap,output);
  }
  
  return true;
}

bool mix::OverlayRawDataDetailMicroBooNE::MixOpDetWaveforms_LowGain( std::vector< std::vector<raw::OpDetWaveform> const* > const& inputs,
								     std::vector<raw::OpDetWaveform> & output,
								     art::PtrRemapper const & remap) {
  
  //make sure we only have two collections for now
  if(inputs.size()!=fEventsToMix || (inputs.size()!=1 && !fInputFileIsData)){
    std::stringstream err_str;
    err_str << "ERROR! We have the wrong number of collections of raw digits we are adding! " << inputs.size();
    throw std::runtime_error(err_str.str());
  }


  if(fInputFileIsData){
    GenerateMCOpDetLowGainScaleMap(*inputOpDetHandle_LowGain);  
    fODMixer.DeclareData(*inputOpDetHandle_LowGain,output);
    for(auto const& icol : inputs)
      fODMixer.Mix(*icol,fMCOpDetLowGainScaleMap,output);
  }
  else if(!fInputFileIsData){
    GenerateMCOpDetLowGainScaleMap(*(inputs[0]));  
    fODMixer.DeclareData(*(inputs[0]),output);
    fODMixer.Mix(*inputOpDetHandle_LowGain,fMCOpDetLowGainScaleMap,output);
  }
  
  return true;
}

DEFINE_ART_MODULE(mix::OverlayRawDataMicroBooNE)
