///////////////////////////////////////////////////////////////////////
/// \file  LArRawInputDriverUBooNE.cxx
/// \brief Source to convert raw binary files to root files
/// \Original Authors
/// \version $Id: LArRawInputDriver.h,v 1.7 2010/01/14 19:20:33 brebel Exp $
/// \author  brebel@fnal.gov, soderber@fnal.gov
/// \MicroBooNE Author: jasaadi@fnal.gov, zarko@fnal.gov (with much help from Wes and Eric)
////////////////////////////////////////////////////////////////////////

//LArSoft 
#include "uboone/RawData/utils/LArRawInputDriverUBooNE.h"
#include "lardata/RawData/RawDigit.h"
#include "lardata/RawData/TriggerData.h"
#include "lardata/RawData/DAQHeader.h"
#include "lardata/RawData/BeamInfo.h"
#include "lardata/RawData/OpDetWaveform.h"
#include "larcore/SummaryData/RunData.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfo/ElecClock.h"
#include "lardata/OpticalDetectorData/OpticalTypes.h" // I want to move the enums we use back to UBooNE as they are UBooNE-specific
#include "uboone/TriggerSim/UBTriggerTypes.h"

//ART, ...
#include "art/Framework/IO/Sources/put_product_in_principal.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Utilities/Exception.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//uboone datatypes

// uboonecode
#include "uboone/Geometry/UBOpChannelTypes.h"
#include "uboone/Geometry/UBOpReadoutMap.h"

#include "datatypes/raw_data_access.h"

//boost
//#include <boost/archive/binary_iarchive.hpp>
#include <boost/algorithm/string.hpp>

//root
#include "TH1D.h"
#include "TTree.h"

//other
#include <libpq-fe.h>

extern "C" {
#include <sys/types.h>
}

#include <bitset>
#include <iostream>

namespace ubdaq=gov::fnal::uboone::datatypes;

// +++ Blatant stealing from LongBo
namespace lris {

  unsigned int LArRawInputDriverUBooNE::RollOver(unsigned int ref,
						 unsigned int subject,
						 unsigned int nbits) 
  {
    // Return "ref" which lower "nbits" are replaced by that of "subject"                                      
    // Takes care of roll over effect.                                                                         
    // For speed purpose we only accept pre-defined nbits values.                                              
  
    unsigned int diff=0; // max diff should be (2^(nbits)-2)/2                                                       
    unsigned int mask=0; // mask to extract lower nbits from subject ... should be 2^(nbits)-1                       
    if      (nbits==3) {diff = 3; mask = 0x7;}
    else if (nbits==4) {diff = 0x7; mask = 0xf;}
    //else if (nbits==4) {diff = 0x7; mask = 0xf;}                                                             
    //    else if (nbits==4) {nbits=3; diff = 0x7; mask = 0x7;}                                                
    else {
      std::cerr<<"\033[93m<<ERROR>>\033[00m"<<" Only supported for nbits = [3,4]!"<<std::endl;
      throw std::exception();
    }
  
    subject = ( (ref>>nbits) << nbits) + (subject & mask);
    //std::cout<<ref<<" : "<<subject<<" : "<<nbits<< " : "<<diff<<std::endl;                                   
    // If exactly same, return                                                                                 
    if(subject == ref) return subject;
  
    // If subject is bigger than ref by a predefined diff value, inspect difference                            
    else if ( subject > ref && (subject - ref) > diff) {
    
      // Throw an exception if difference is exactly diff+1                                                    
      if ( (subject - ref) == diff+1 ) {
	std::cerr<<"\033[93m<<ERROR>>\033[00m"<<Form(" Unexpected diff: ref=%d, subject=%d",ref,subject)<<std::endl;
	throw std::exception();
      }

      // Else we have to subtract (mask+1)                                                                     
      else{
	//std::cout<<Form("Correcting %d to %d",subject,(subject-(mask+1)))<<std::endl;                        
	subject = subject - (mask + 1);
      }
    
    }
    // If subject is smaller than ref by a predefined diff value, inspect difference                           
    else if ( subject < ref && (ref - subject) > diff) {
    
      // Throw an exception if difference is exactly diff+1                                                    
      if ( (ref - subject) == diff+1 ) {
	std::cerr<<"\033[93m<<ERROR>>\033[00m"<<Form(" Unexpected diff: ref=%d, subject=%d",ref,subject)<<std::endl;
	throw std::exception();
      }
      else{
	//std::cout<<Form("Correcting %d to %d",subject,(subject + (mask+1)))<<std::endl;                      
	subject = subject + (mask + 1);
      }
    }
    return subject;
  }

  // ======================================================================
  LArRawInputDriverUBooNE::LArRawInputDriverUBooNE(fhicl::ParameterSet const & ps, 
                                                   art::ProductRegistryHelper &helper,
                                                   art::SourceHelper const &pm)
    :
    fSourceHelper(pm),
    fCurrentSubRunID(),
    fEventCounter(0),
    fFinalEventCutOff(30000000), //this value of 30 000 000 is the approximate size of an event.
    fPreviousPosition(0),
    fCompleteFile(true),
    fHuffmanDecode(ps.get<bool>("huffmanDecode",false)),
    fUseGPS(ps.get<bool>("UseGPS",false)),
    fUseNTP(ps.get<bool>("UseNTP",false)),
    fMaxEvents(-1),
    fSkipEvents(0)
  {
    ::peek_at_next_event<ub_TPC_CardData_v6>(false);
    ::peek_at_next_event<ub_PMT_CardData_v6>(false);
    ::handle_missing_words<ub_TPC_CardData_v6>(true);
    ::handle_missing_words<ub_PMT_CardData_v6>(true);

    helper.reconstitutes<raw::DAQHeader,                 art::InEvent>("daq");
    helper.reconstitutes<std::vector<raw::RawDigit>,     art::InEvent>("daq");
    helper.reconstitutes<raw::BeamInfo,                  art::InEvent>("daq");
    helper.reconstitutes<std::vector<raw::Trigger>,      art::InEvent>("daq");
    helper.reconstitutes<raw::ubdaqSoftwareTriggerData, art::InEvent>("daq");
    registerOpticalData( helper ); //helper.reconstitutes<std::vector<raw::OpDetWaveform>,art::InEvent>("daq");
    fDataTakingTime                    = ps.get< int  >("DataTakingTime", -1);
    fSwizzlingTime                     = ps.get< int  >("SwizzlingTime", -1);

    fSwizzleTPC = ps.get<bool>("swizzleTPC",true);
    fSwizzlePMT = ps.get<bool>("swizzlePMT",true);
    fSwizzlePMT_init = ps.get<bool>("swizzlePMT",true);
    fSwizzleTrigger = ps.get<bool>("swizzleTrigger",true);
    fSwizzleTriggerType = ps.get<std::string>("swizzleTriggerType"); // Only use ALL for this option, other options will not work
    fMaxEvents = ps.get<int>("maxEvents", -1);
    fSkipEvents = ps.get<int>("skipEvents", 0);

    //temporary kazuTestSwizzleTrigger
    kazuTestSwizzleTrigger = ps.get<bool>("kazuTestSwizzleTrigger",true);
   
    //if ( fHuffmanDecode )
    tpc_crate_data_t::doDissect(true); // setup for decoding

    art::ServiceHandle<art::TFileService> tfs;
    //initialize beam histograms specified in fhicl file
    art::TFileDirectory tfbeamdir = tfs->mkdir( "Beam" );
    std::vector<std::string> beam_hist=ps.get<std::vector<std::string> >("beam_histograms");
    for ( auto it : beam_hist ) {
      std::vector<std::string> hist;
      boost::split(hist, it, boost::is_any_of(","));
      if (hist.size() != 4)
        mf::LogWarning("") << "Bad definition in fhicl file for histogram "<<hist.at(0)<<". Ignoring it.";
      else {
        TH1D* h=tfbeamdir.make<TH1D>(hist[0].c_str(),hist[0].c_str(),
				     atoi(hist[1].c_str()),atof(hist[2].c_str()),atof(hist[3].c_str()));
        std::pair<std::string, TH1D*> p(hist[0],h);
        fHistMapBeam.insert(p);
      }
    }

    art::TFileDirectory tfdebugdir = tfs->mkdir( "Debug" );
    tMyTree = tfdebugdir.make<TTree>("tMyTree", "tree");
    tMyTree->Branch("event",&event,"event/I");
    tMyTree->Branch("triggerFrame",&triggerFrame,"triggerFrame/I");
    tMyTree->Branch("triggerSample",&triggerSample,"triggerSample/I");
    tMyTree->Branch("triggerTime",&triggerTime,"triggerTime/D");
    tMyTree->Branch("triggerActive",&triggerActive,"triggerActive/I");
    tMyTree->Branch("triggerBitBNB",&triggerBitBNB,"triggerBitBNB/I");
    tMyTree->Branch("triggerBitNuMI",&triggerBitNuMI,"triggerBitNuMI/I");
    tMyTree->Branch("triggerBitEXT",&triggerBitEXT,"triggerBitEXT/I");
    tMyTree->Branch("triggerBitPMTBeam",&triggerBitPMTBeam,"triggerBitPMTBeam/I");
    tMyTree->Branch("triggerBitPMTCosmic",&triggerBitPMTCosmic,"triggerBitPMTCosmic/I");
    tMyTree->Branch("FEM5triggerFrame",&FEM5triggerFrame,"FEM5triggerFrame/I");
    tMyTree->Branch("FEM5triggerSample",&FEM5triggerSample,"FEM5triggerSample/I");
    tMyTree->Branch("FEM6triggerFrame",&FEM6triggerFrame,"FEM6triggerFrame/I");
    tMyTree->Branch("FEM6triggerSample",&FEM6triggerSample,"FEM6triggerSample/I");
    tMyTree->Branch("FEM5triggerTime",&FEM5triggerTime,"FEM5triggerTime/D");
    tMyTree->Branch("FEM6triggerTime",&FEM6triggerTime,"FEM6triggerTime/D");

    tMyTree->Branch("RO_BNBtriggerFrame",&RO_BNBtriggerFrame,"RO_BNBtriggerFrame/I");
    tMyTree->Branch("RO_NuMItriggerFrame",&RO_NuMItriggerFrame,"RO_NuMItriggerFrame/I");
    tMyTree->Branch("RO_EXTtriggerFrame",&RO_EXTtriggerFrame,"RO_EXTtriggerFrame/I");
    tMyTree->Branch("RO_RWMtriggerFrame",&RO_RWMtriggerFrame,"RO_RWMtriggerFrame/I");
    tMyTree->Branch("RO_BNBtriggerSample",&RO_BNBtriggerSample,"RO_BNBtriggerSample/I");
    tMyTree->Branch("RO_NuMItriggerSample",&RO_NuMItriggerSample,"RO_NuMItriggerSample/I");
    tMyTree->Branch("RO_EXTtriggerSample",&RO_EXTtriggerSample,"RO_EXTtriggerSample/I");
    tMyTree->Branch("RO_RWMtriggerSample",&RO_RWMtriggerSample,"RO_RWMtriggerSample/I");
    tMyTree->Branch("RO_BNBtriggerTime",&RO_BNBtriggerTime,"RO_BNBtriggerTime/D");
    tMyTree->Branch("RO_NuMItriggerTime",&RO_NuMItriggerTime,"RO_NuMItriggerTime/D");
    tMyTree->Branch("RO_EXTtriggerTime",&RO_EXTtriggerTime,"RO_EXTtriggerTime/D");
    tMyTree->Branch("RO_RWMtriggerTime",&RO_RWMtriggerTime,"RO_RWMtriggerTime/D");

//    tMyTree->Branch("RO_Gate1Frame",&RO_Gate1Frame,"RO_Gate1Frame/I");
//    tMyTree->Branch("RO_Gate1Sample",&RO_Gate1Sample,"RO_Gate1Sample/I");
//    tMyTree->Branch("RO_Gate2Frame",&RO_Gate2Frame,"RO_Gate2Frame/I");
//    tMyTree->Branch("RO_Gate2Sample",&RO_Gate2Sample,"RO_Gate2Sample/I");

    tMyTree->Branch("TPC1triggerFrame",&TPC1triggerFrame,"TPC1triggerFrame/I");
    tMyTree->Branch("TPC1triggerSample",&TPC1triggerSample,"TPC1triggerSample/I");
    tMyTree->Branch("TPC2triggerFrame",&TPC2triggerFrame,"TPC2triggerFrame/I");
    tMyTree->Branch("TPC2triggerSample",&TPC2triggerSample,"TPC2triggerSample/I");
    tMyTree->Branch("TPC3triggerFrame",&TPC3triggerFrame,"TPC3triggerFrame/I");
    tMyTree->Branch("TPC3triggerSample",&TPC3triggerSample,"TPC3triggerSample/I");
    tMyTree->Branch("TPC4triggerFrame",&TPC4triggerFrame,"TPC4triggerFrame/I");
    tMyTree->Branch("TPC4triggerSample",&TPC4triggerSample,"TPC4triggerSample/I");
    tMyTree->Branch("TPC5triggerFrame",&TPC5triggerFrame,"TPC5triggerFrame/I");
    tMyTree->Branch("TPC5triggerSample",&TPC5triggerSample,"TPC5triggerSample/I");
    tMyTree->Branch("TPC6triggerFrame",&TPC6triggerFrame,"TPC6triggerFrame/I");
    tMyTree->Branch("TPC6triggerSample",&TPC6triggerSample,"TPC6triggerSample/I");
    tMyTree->Branch("TPC7triggerFrame",&TPC7triggerFrame,"TPC7triggerFrame/I");
    tMyTree->Branch("TPC7triggerSample",&TPC7triggerSample,"TPC7triggerSample/I");
    tMyTree->Branch("TPC8triggerFrame",&TPC8triggerFrame,"TPC8triggerFrame/I");
    tMyTree->Branch("TPC8triggerSample",&TPC8triggerSample,"TPC8triggerSample/I");
    tMyTree->Branch("TPC9triggerFrame",&TPC9triggerFrame,"TPC9triggerFrame/I");
    tMyTree->Branch("TPC9triggerSample",&TPC9triggerSample,"TPC9triggerSample/I");

    tMyTree->Branch("N_discriminators",N_discriminators,"N_discriminators[40]/I");
    tMyTree->Branch("discriminatorSample",discriminatorSample,"discriminatorSample[40][100]/I");
    tMyTree->Branch("discriminatorFrame",discriminatorFrame,"discriminatorFrame[40][100]/I");
    tMyTree->Branch("discriminatorType",discriminatorType,"discriminatorType[40][100]/I");
    
    tMyTree->Branch("N_PMT_waveforms",&N_PMT_waveforms,"N_PMT_waveforms/I");
    tMyTree->Branch("PMT_waveform_times",PMT_waveform_times,"PMT_waveform_times[N_PMT_waveforms]/D");

    tMyTree->Branch("ADCwords_crate0",&ADCwords_crate0,"ADCwords_crate0/I");
    tMyTree->Branch("ADCwords_crate1",&ADCwords_crate1,"ADCwords_crate1/I");
    tMyTree->Branch("ADCwords_crate2",&ADCwords_crate2,"ADCwords_crate2/I");
    tMyTree->Branch("ADCwords_crate3",&ADCwords_crate3,"ADCwords_crate3/I");
    tMyTree->Branch("ADCwords_crate4",&ADCwords_crate4,"ADCwords_crate4/I");
    tMyTree->Branch("ADCwords_crate5",&ADCwords_crate5,"ADCwords_crate5/I");
    tMyTree->Branch("ADCwords_crate6",&ADCwords_crate6,"ADCwords_crate6/I");
    tMyTree->Branch("ADCwords_crate7",&ADCwords_crate7,"ADCwords_crate7/I");
    tMyTree->Branch("ADCwords_crate8",&ADCwords_crate8,"ADCwords_crate8/I");
    tMyTree->Branch("ADCwords_crate9",&ADCwords_crate9,"ADCwords_crate9/I");
    tMyTree->Branch("NumWords_crate1",&NumWords_crate1,"NumWords_crate1/I");
    tMyTree->Branch("NumWords_crate2",&NumWords_crate2,"NumWords_crate2/I");
    tMyTree->Branch("NumWords_crate3",&NumWords_crate3,"NumWords_crate3/I");
    tMyTree->Branch("NumWords_crate4",&NumWords_crate4,"NumWords_crate4/I");
    tMyTree->Branch("NumWords_crate5",&NumWords_crate5,"NumWords_crate5/I");
    tMyTree->Branch("NumWords_crate6",&NumWords_crate6,"NumWords_crate6/I");
    tMyTree->Branch("NumWords_crate7",&NumWords_crate7,"NumWords_crate7/I");
    tMyTree->Branch("NumWords_crate8",&NumWords_crate8,"NumWords_crate8/I");
    tMyTree->Branch("NumWords_crate9",&NumWords_crate9,"NumWords_crate9/I");

    event = 0;

  }
  
  
  // ======================================================================
  void LArRawInputDriverUBooNE::closeCurrentFile()  
  {    
    mf::LogInfo(__FUNCTION__)<<"File boundary (processed "<<fEventCounter<<" events)"<<std::endl;
    fCurrentSubRunID.flushSubRun();
    fEventCounter=0;
    fNumberEventsInFile=0;
    fInputStream.close();
  }
    
  // ======================================================================
  void LArRawInputDriverUBooNE::readFile(std::string const &name,
                                         art::FileBlock* &fb)
  {
    // Fill and return a new Fileblock.
    fb = new art::FileBlock(art::FileFormatVersion(1, "LArRawInput 2011a"),
                            name);
    
    fInputStream.open(name.c_str(),std::ios_base::in | std::ios_base::binary);
    
    
    // Throwing an exception if the file fails to open
    if( !fInputStream.is_open() ) {
      throw art::Exception( art::errors::FileReadError )
        << "failed to open input file " << name << std::endl;
    }
      
    //seek to the end of file, check the end word, check number of events and sizes
    //read in the end of file word first
    if (fEventCounter==0) {
      uint16_t end_of_file_marker;
      uint32_t nevents;
      fPreviousPosition = fInputStream.tellg();
      fInputStream.seekg( -1*sizeof(uint16_t) , std::ios::end); //eof marker is 16 bits long so go 16 bits from the end.
      fInputStream.read( (char*)&end_of_file_marker , sizeof(uint16_t));
      if(end_of_file_marker != 0xe0f0){ //make sure that it is the correct marker
          fCompleteFile=false;
          fNumberEventsInFile = -1;
          mf::LogWarning("")<< "File "<<name<<" has incorrect end of file marker. Expected 0xe0f0 and instead got 0x"<< std::hex << end_of_file_marker<<std::endl;
          mf::LogWarning("")<< "The number of events in this file will be determined on the fly, and we will stop when there is less than " << fFinalEventCutOff << " bytes in the file." <<std::endl;
      } else {
          fCompleteFile=true;
          mf::LogInfo("")<<"Opened file "<<name<<" and found that is has a good end of file marker."<<std::endl;
      }
      
      if(fCompleteFile) {
          fInputStream.seekg( -3*sizeof(uint16_t) , std::ios::end); //need to go 48 bits from the end of the file (16 for eof marker and 32 for nevents
          fInputStream.read( (char*)&nevents, sizeof(uint32_t)); //read in the 32 bits that should be the event count
          if (nevents>0 && nevents <1E9) { //make sure that nevents is reasonable
              mf::LogInfo("")<<"Opened file "<<name<<" with "<< nevents <<" event(s)";
              fNumberEventsInFile = nevents;
          } else {
              throw art::Exception( art::errors::FileReadError )
              << "File "<<name<<" has incorrect number of events in trailer. "<< nevents<<std::endl;
          }
      }
    
      fInputStream.seekg(std::ios::beg);
    }
    
    return;

  }


  // =====================================================================
  void LArRawInputDriverUBooNE::registerOpticalData( art::ProductRegistryHelper &helper ) {
    // we make a data product for each category of channels
    fPMTdataProductNames.clear();
    for ( unsigned int cat=0; cat<(unsigned int)opdet::NumUBOpticalChannelCategories; cat++ ) {
      helper.reconstitutes<std::vector<raw::OpDetWaveform>,art::InEvent>( "pmtreadout", opdet::UBOpChannelEnumName( (opdet::UBOpticalChannelCategory_t)cat ) );
      fPMTdataProductNames.insert( std::make_pair( (opdet::UBOpticalChannelCategory_t)cat, opdet::UBOpChannelEnumName( (opdet::UBOpticalChannelCategory_t)cat ) ) );
    }
  }

  // =====================================================================
  void LArRawInputDriverUBooNE::putPMTDigitsIntoEvent( std::map< opdet::UBOpticalChannelCategory_t, std::unique_ptr< std::vector<raw::OpDetWaveform> > >& pmtdigitlist, 
						       art::EventPrincipal* &outE ) {
    for ( unsigned int cat=0; cat<(unsigned int)opdet::NumUBOpticalChannelCategories; cat++ ) {
      
      art::put_product_in_principal(std::move( pmtdigitlist[(opdet::UBOpticalChannelCategory_t)cat]  ),
      				    *outE,
				    "pmtreadout", // module
      				    fPMTdataProductNames[ (opdet::UBOpticalChannelCategory_t)cat ]); // instance
    }
    
  }
  
  // =====================================================================
  bool LArRawInputDriverUBooNE::readNext(art::RunPrincipal* const &/*inR*/,
                                         art::SubRunPrincipal* const &/*inSR*/,
                                         art::RunPrincipal* &outR,
                                         art::SubRunPrincipal* &outSR,
                                         art::EventPrincipal* &outE)
  {
      if (fCompleteFile) {
          if (fEventCounter==fNumberEventsInFile) {
              mf::LogInfo(__FUNCTION__)<<"Already read " << fEventCounter << " events, so checking end of file..." << std::endl;
              std::ios::streampos current_position = fInputStream.tellg(); //find out where in the file we're located
              std::ios::streampos data_offset = current_position - fPreviousPosition;
              mf::LogInfo(__FUNCTION__)<< "For event " << fEventCounter << ", the number of bytes read since the last event is " << data_offset << std::endl;
              fPreviousPosition = current_position;
              fInputStream.seekg(0,std::ios::end); //go to the end of the file
              std::ios::streampos file_length = fInputStream.tellg(); //get the location which will tell the size.
              fInputStream.seekg(current_position); //put the ifstream back to where it was before
              if ( ((uint8_t)file_length - (uint8_t)current_position) > (uint8_t)1000 ) {
                  throw art::Exception( art::errors::FileReadError ) << "We processed " << fEventCounter << " events from the file " <<  std::endl << "But there are still " << (file_length - current_position) << " bytes in the file" << std::endl;
              }
              mf::LogInfo(__FUNCTION__)<<"Completed reading file and closing output file." << std::endl;
              return false; //tells readNext that you're done reading all of the events in this file.
          } else {
              std::ios::streampos current_position = fInputStream.tellg(); //find out where in the file we're located
              std::ios::streampos data_offset = current_position - fPreviousPosition;
              mf::LogInfo(__FUNCTION__)<< "For event " << fEventCounter << ", the number of bytes read since the last event is " << data_offset << std::endl;
              fPreviousPosition = current_position;
	  }
      } else {
          mf::LogInfo(__FUNCTION__)<<"Read " << fEventCounter << " events from an incomplete file, and checking end of file..." << std::endl;
          std::ios::streampos current_position = fInputStream.tellg(); //find out where in the file we're located
          std::ios::streampos data_offset = current_position - fPreviousPosition;
          mf::LogInfo(__FUNCTION__)<< "For event " << fEventCounter << ", the number of bytes read since the last event is " << data_offset << std::endl;
          fPreviousPosition = current_position;
          fInputStream.seekg(0,std::ios::end); //go to the end of the file
          std::ios::streampos file_length = fInputStream.tellg(); //get the location which will tell the size.
          fInputStream.seekg(current_position); //put the ifstream back to where it was before
          if ( (file_length - current_position) < fFinalEventCutOff ) { //this value of 30000000 is the approximate size of an event.
              mf::LogInfo(__FUNCTION__) << "We processed " << fEventCounter << " events from the incomplete file " << std::endl << "But there are only " << (file_length - current_position) << " bytes in the file, so we're calling that the end of the road." << std::endl;
              return false; //tells readNext that you're done reading all of the events in this file.
          } else {
              mf::LogInfo(__FUNCTION__)<<"We processed " << fEventCounter << " events from the incomplete file. " << std::endl << "There are " << (file_length - current_position) << " bytes remaining to be read in the file, so we're gonna keep processing..." << std::endl;
          }
      }

    if (fMaxEvents > 0 && fEventCounter == unsigned(fMaxEvents))
      return false;

    mf::LogInfo(__FUNCTION__)<<"Attempting to read event: "<<fEventCounter<<std::endl;
    // Create empty result, then fill it from current file:
    std::unique_ptr<raw::DAQHeader> daq_header(new raw::DAQHeader);
    std::unique_ptr<std::vector<raw::RawDigit> >  tpc_raw_digits( new std::vector<raw::RawDigit>  );
    std::unique_ptr<raw::BeamInfo> beam_info(new raw::BeamInfo);
    std::unique_ptr<std::vector<raw::Trigger>> trig_info( new std::vector<raw::Trigger> );
    std::map< opdet::UBOpticalChannelCategory_t, std::unique_ptr< std::vector<raw::OpDetWaveform> > > pmt_raw_digits;
    std::unique_ptr<raw::ubdaqSoftwareTriggerData> sw_trig_info( new raw::ubdaqSoftwareTriggerData );
    //raw::ubdaqSoftwareTriggerData * sw_trig_info(new raw::ubdaqSoftwareTriggerData);

    for ( unsigned int opdetcat=0; opdetcat<(unsigned int)opdet::NumUBOpticalChannelCategories; opdetcat++ ) {
      pmt_raw_digits.insert( std::make_pair( (opdet::UBOpticalChannelCategory_t)opdetcat, std::unique_ptr< std::vector<raw::OpDetWaveform> >(  new std::vector<raw::OpDetWaveform> ) ) );
    }
    event  = fEventCounter;
    NumWords_crate1 = 0;
    NumWords_crate2 = 0;
    NumWords_crate3 = 0;
    NumWords_crate4 = 0;
    NumWords_crate5 = 0;
    NumWords_crate6 = 0;
    NumWords_crate7 = 0;
    NumWords_crate8 = 0;
    NumWords_crate9 = 0;
    ADCwords_crate0 = 0;
    ADCwords_crate1 = 0;
    ADCwords_crate2 = 0;
    ADCwords_crate3 = 0;
    ADCwords_crate4 = 0;
    ADCwords_crate5 = 0;
    ADCwords_crate6 = 0;
    ADCwords_crate7 = 0;
    ADCwords_crate8 = 0;
    ADCwords_crate9 = 0;

    uint32_t event_number = 0;

    bool res = false;
    bool done = false;

    while(!done) {
      res=processNextEvent(*tpc_raw_digits, pmt_raw_digits, *daq_header, *beam_info, *trig_info, *sw_trig_info, event_number, fSkipEvents > 0);
      if(fSkipEvents > 0)
	--fSkipEvents;
      else
	done = true;
    }

    if (res) {
      fEventCounter++;
      art::RunNumber_t rn = daq_header->GetRun();//+1;
      art::Timestamp tstamp = daq_header->GetTimeStamp();
      art::SubRunID newID(rn, daq_header->GetSubRun());
      if (fCurrentSubRunID.runID() != newID.runID()) { // New Run
        outR = fSourceHelper.makeRunPrincipal(rn, tstamp);
      }
      if (fCurrentSubRunID != newID) { // New SubRun
        outSR = fSourceHelper.makeSubRunPrincipal(rn,
                                                  daq_header->GetSubRun(),
                                                  tstamp);
        fCurrentSubRunID = newID;        
      }
      /*
	std::cout<<"\033[93mAbout to make a principal for run: " << fCurrentSubRunID.run()
	<<" subrun: " << fCurrentSubRunID.subRun()
	<<" event: " << event_number
	<<"\033[00m"<< std::endl;
      */
      outE = fSourceHelper.makeEventPrincipal(fCurrentSubRunID.run(),
					      fCurrentSubRunID.subRun(),
					      event_number,
					      tstamp);
      //std::cout<<"\033[93mDone\033[00m"<<std::endl;

      // Put products in the event.
      art::put_product_in_principal(std::move(tpc_raw_digits),
                                    *outE,
                                    "daq"); // Module label
      art::put_product_in_principal(std::move(daq_header),
                                    *outE,
                                    "daq"); // Module label
      art::put_product_in_principal(std::move(beam_info),
                                    *outE,
                                    "daq"); // Module label
      art::put_product_in_principal(std::move(trig_info),
				    *outE,
				    "daq"); // Module label
      art::put_product_in_principal(std::move(sw_trig_info),
				    *outE,
				    "daq"); // Module label
      putPMTDigitsIntoEvent( pmt_raw_digits, outE );
     
    }

    return res;

  }

  // =====================================================================
  bool LArRawInputDriverUBooNE::processNextEvent(std::vector<raw::RawDigit>& tpcDigitList,
                                                 std::map< opdet::UBOpticalChannelCategory_t, 
                                                 std::unique_ptr<std::vector<raw::OpDetWaveform>> >& pmtDigitList,
                                                 raw::DAQHeader& daqHeader,
                                                 raw::BeamInfo& beamInfo,
						 std::vector<raw::Trigger>& trigInfo,
						 raw::ubdaqSoftwareTriggerData& sw_trigInfo,
						 uint32_t& event_number,
						 bool skip)
  {  
     triggerFrame = -999;
     TPCframe = -999;
     PMTframe = -999;
     FEM5triggerSample=-999;
     FEM6triggerSample=-999;
     FEM5triggerFrame=-999;
     FEM6triggerFrame=-999;
     RO_BNBtriggerFrame=-999;
     RO_BNBtriggerSample=-999;
     RO_NuMItriggerFrame=-999;
     RO_NuMItriggerSample=-999;
     RO_EXTtriggerFrame=-999;
     RO_EXTtriggerSample=-999;
     RO_RWMtriggerFrame=-999;
     RO_RWMtriggerSample=-999;
    
     skipEvent = false;

//     RO_Gate1Frame=-999;
//     RO_Gate1Sample=-999;
//     RO_Gate2Frame=-999;
//     RO_Gate2Sample=-999;

    //try {
    boost::archive::binary_iarchive ia(fInputStream); 
    ubdaq::ub_EventRecord event_record;  
    ia >> event_record;
    if(skip)
      return false;
    //std::cout<<event_record.debugInfo()<<std::endl;
    //set granularity 
    //      event_record.updateIOMode(ubdaq::IO_GRANULARITY_CHANNEL);
    _trigger_beam_window_time = std::numeric_limits<double>::max();
    fillTriggerData(event_record, trigInfo);
    //if (skipEvent){return false;} // check that trigger data doesn't suggest we should skip event. // commented out because this doesn't work at the moment
    fillDAQHeaderData(event_record, daqHeader);
    fillTPCData(event_record, tpcDigitList);
    //please keep fillPMTData ahead of fillSWTriggerData in cases of events without any PMT data
    fillPMTData(event_record, pmtDigitList);
    fillBeamData(event_record, beamInfo);
    //please keep fillPMTData ahead of fillSWTriggerData in cases of events without any PMT data
    fillSWTriggerData(event_record, sw_trigInfo);
      
    event_number = event_record.getGlobalHeader().getEventNumber()+1;

    checkTimeStampConsistency();
      
    tMyTree->Fill();
      
    //Note that the fSwizzlePMT value needs to be reset every event, since it can be set to false if there is no PMT data
    //Setting it to false makes sure that the checkTimeStampConsistency() doesn't try to compare nonexistent data

    fSwizzlePMT = fSwizzlePMT_init;
      
    /*
      } catch (...) {
      //throw art::Exception( art::errors::FileReadError )
      std::cout<< "\033[93mFailed to read the event.\033[00m\n"<< std::endl;
      return false;
      }
    */  
    return true;
  }
  
  // =====================================================================
  void LArRawInputDriverUBooNE::fillDAQHeaderData(ubdaq::ub_EventRecord& event_record,
                                                  raw::DAQHeader& daqHeader)
  {
    ubdaq::ub_GlobalHeader global_header = event_record.getGlobalHeader();
    if(fUseGPS)
      global_header.useGPSTime();
    else if(fUseNTP)
      global_header.useLocalHostTime();

    // art::Timestamp is an unsigned long long. The conventional 
    // use is for the upper 32 bits to have the seconds since 1970 epoch 
    // and the lower 32 bits to be the number of nanoseconds within the 
    // current second.
    // (time_t is a 64 bit word)

      uint32_t seconds=global_header.getSeconds();
      uint32_t nano_seconds=global_header.getNanoSeconds()+
	                    global_header.getMicroSeconds()*1000;
      time_t mytime = ((time_t)seconds<<32) | nano_seconds;

    //\/      uint32_t subrun_num = global_header->getSubrunNumber();
      
    daqHeader.SetStatus(1);
    daqHeader.SetFileFormat(global_header.getRecordType());
    daqHeader.SetSoftwareVersion(global_header.DAQ_version_number);
    daqHeader.SetRun(global_header.getRunNumber());
    daqHeader.SetSubRun(global_header.getSubrunNumber());
      
    //\/ Add the subRun number too!
    daqHeader.SetEvent(global_header.getEventNumber()+1);
    daqHeader.SetTimeStamp(mytime);

    /// \todo: What is the "fixed word" ? Leaving it unset for now
    /// \todo: What is the "spare word" ? Leaving it unset for now
    //daqHeader.SetFixedWord(h1.header);
    //daqHeader.SetSpareWord(h1.spare);
  }

  // =====================================================================
  void LArRawInputDriverUBooNE::fillTPCData(ubdaq::ub_EventRecord& event_record,
                                            std::vector<raw::RawDigit>& tpcDigitList)
  {    
    //Channel map has changed each time the detector has been re-cabled.
    //Provide data-taking time as first argument. (integer epoch seconds) 
    //Optionally recover outdated mappings with 'swizzling time' second arg. (also integer epoch seconds)
    if (fDataTakingTime == -1)
      fChannelMap = art::ServiceHandle<util::DatabaseUtil>()->GetUBChannelMap(event_record.LocalHostTime().seb_time_sec, fSwizzlingTime); 
    else
      fChannelMap = art::ServiceHandle<util::DatabaseUtil>()->GetUBChannelMap(fDataTakingTime, fSwizzlingTime); 

    
    // ### Swizzling to get the number of channels...trying the method used in write_read.cpp
    // ### provided by Wes --- About the data:
    // ### The format of the data is in levels: crate, card, channel.
    // ### Level 1: The event record contains a map of (crateHeader,crateData) pairs.
    // ### Level 2: Each crateData object may contain a map of (cardHeader,cardData) pairs.
    // ### Level 3: Each cardData object may contain a map of (int,channelData) pairs. The int
    // ### is the channel number.
      
    //get the seb map, and do a loop over all sebs/crates

    for( auto const& seb_it : event_record.getTPCSEBMap()) {    // I think auto should be tpc_map_t::const_iterator  -NJT

          
      //get the crateHeader/crateData objects
      // ubdaq::crateHeader crate_header = seb_it->first;
      //        ubdaq::crateData crate_data = seb_it->second;
      //      int tpc_seb_num = seb_it.first;
      tpc_crate_data_t const& tpc_crate = seb_it.second; // (ub_TPC_CrateData_v6)
      int crate_number = seb_it.first; // confirmed this is correct. for now, crate's are not given their ID number to store and then retrieve.
      
      if ( !tpc_crate.wasDissected() ) {
	std::cerr << "Warning crate data corrupted! Skipping." << std::endl;
	tpc_crate.dissectionException().what();
	continue;
      }

      // Get Time information:
      //uint32_t sebTSec = crate_header.getSebTimeSec();
      //std::cout << "Seb Time (sec) : " << sebTSec << std::endl;
      //jmsj crate_header_t crHeader = crate_header.getCrateHeader();
      // GPStime in UNIX second/micro/nano info
      //jmsjgps_time_t GPStime = crHeader.gps_time;
      // DAQtime is time of last update of GPS time (in frame, sample, div)
      //tbclkub_t  DAQtime = crHeader.daqClock_time;
      //std::cout << "GPS Time seconds: " << GPStime.second << std::endl;
      //std::cout << "DAQ Frame: " << DAQtime.frame << "\tSample: " << DAQtime.sample << std::endl;

      //      auto const& tpc_crate_header = tpc_crate.header();    
      //      auto const& tpc_crate_trailer = tpc_crate.trailer();  

      //Special to the crate, there is a special header that the DAQ attaches. You can access this
      //like so. The type here is a unique ptr to a ub_CrateHeader_v6 struct. That has useful info
      //like the local host time, which may or may not be set properly right now...
      //auto const& tpc_crate_DAQ_header = tpc_crate.crateHeader(); // I think auto should be tpc_crate_data_t::ub_CrateHeader_t --NJT
      //     ub_LocalHostTime this_time = tpc_crate_DAQ_header->local_host_time;
      
      //The Crate Data is split up into Cards. You use the "getCards()" command to get access to
      //each of those. Note that calling this function will dissect the data if it has not already
      //been dissected (debugInfo() calls getCards()). You can do a look over the cards like so:
      for(auto const& card : tpc_crate.getCards()){  // This auto is tpc_crate_data_t::card_t

	if ( !card.wasDissected() ) {
	  std::cerr << "Warning card data corrupted! Skipping." << std::endl;
	  card.dissectionException().what();
	  continue;
	}

        //The format here is similar to the crate! There's a header (which is a ub_TPC_CardHeader_v*
        //object), and technically a trailer (though here it's empty!).
	auto const& tpc_card_header = card.header();   
        
        unsigned int frame = RollOver(tpc_card_header.getFrame(), tpc_card_header.getTrigFrameMod16(), 3);
        
        if (TPCframe == -999){TPCframe = frame;} // internal frame consistency checking
        if (abs(frame - TPCframe) > 1){ // if the frame doesn't match the other TPC frames here then we have a problem
          std::cerr << "ERROR!" << std::endl;
          std::cerr << "TPC card header trigger frames not within one frame of each other!!" << std::endl;
          throw std::exception();
        }

//        if (fSwizzleTrigger){
//          if (triggerFrame == -999){ // if we have swizzled trigger data, and we have the default triggerFrame, we have a problem
//            std::cerr << "ERROR!" << std::endl;
//            std::cerr << "Trigger data should have been swizzled, but triggerFrame remains as the default value of -999" << std::endl;
//            throw std::exception();
//          }
//          if (abs(frame - triggerFrame) > 1){ // if the triggerFrame and frame don't agree here then we have a problem
//            std::cerr << "ERROR!" << std::endl;
//            std::cerr << "TPC card header trigger frame not within one frame of trigger header frame!! Trigger header frame is " 
//                      << triggerFrame << " and TPC card header frame is " << frame << std::endl;
//            throw std::exception();
//          }
//        }
        // Output tree variables - for calculating compression
        if (crate_number == 1){
          NumWords_crate1 += tpc_card_header.getWordCount();
          TPC1triggerFrame = frame;
          TPC1triggerSample = tpc_card_header.getTrigSample();
        }
        if (crate_number == 2){
          NumWords_crate2 += tpc_card_header.getWordCount();
          TPC2triggerFrame = frame;
          TPC2triggerSample = tpc_card_header.getTrigSample();
        }
        if (crate_number == 3){
          NumWords_crate3 += tpc_card_header.getWordCount();
          TPC3triggerFrame = frame;
          TPC3triggerSample = tpc_card_header.getTrigSample();
        }
        if (crate_number == 4){
          NumWords_crate4 += tpc_card_header.getWordCount();
          TPC4triggerFrame = frame;
          TPC4triggerSample = tpc_card_header.getTrigSample();
        }
        if (crate_number == 5){
          NumWords_crate5 += tpc_card_header.getWordCount();
          TPC5triggerFrame = frame;
          TPC5triggerSample = tpc_card_header.getTrigSample();
        }
        if (crate_number == 6){
          NumWords_crate6 += tpc_card_header.getWordCount();
          TPC6triggerFrame = frame;
          TPC6triggerSample = tpc_card_header.getTrigSample();
        }
        if (crate_number == 7){
          NumWords_crate7 += tpc_card_header.getWordCount();
          TPC7triggerFrame = frame;
          TPC7triggerSample = tpc_card_header.getTrigSample();
        }
        if (crate_number == 8){
          NumWords_crate8 += tpc_card_header.getWordCount();
          TPC8triggerFrame = frame;
          TPC8triggerSample = tpc_card_header.getTrigSample();
        }
        if (crate_number == 9){
          NumWords_crate9 += tpc_card_header.getWordCount();
          TPC9triggerFrame = frame;
          TPC9triggerSample = tpc_card_header.getTrigSample();
        }
         

	//	auto const& tpc_card_trailer = card.trailer(); 

        //Of course, you can probe for information in the card header. You'll have to find the appropriate
        //header file to know what type you have, but again, these will follow typical practice. And, you
        //can always use debugInfo to not only print the info, but it tells you the type.
        // auto const this_event_number = card.getEvent(); /// auto are ints here
        // auto const this_frame_number = card.getFrame(); /// auto are ints here


        //And, you guessed it, the tpc card data is split up into one more level: by channel.
        for(auto const& channel : card.getChannels()){ // auto here tpc_crate_data_t::card_t::card_channel_type


	  if ( !channel.wasDissected() ) {
	    std::cerr << "Warning channel data corrupted! Skipping." << std::endl;
	    //channel.dissectionException().what();
	    continue;
	  }

	  //There's a header and trailer here. Remember these are just uint16_t, that contain the
	  //channel number.
	  // auto const& tpc_channel_header = channel.header();   // unused
	  // auto const& tpc_channel_trailer = channel.trailer(); // unsued
	  
	  //The channel object (ub_MarkedRawChannelData) has a method for returning the channel.
	  //You can look at the other objects too (like ub_MarkedRawCardData) and see methods of
	  //use there as well.
	  auto const tpc_channel_number = channel.getChannelNumber(); // auto is int here
                        

            // output:
            std::vector<short> adclist;
	    size_t chdsize(0); 

	    //Huffman decoding

	    channel.decompress(adclist); // All-in-one call.
	    
	    /* // Commented out as trailer check is now donw via swizzler
	       uint16_t frailer = channel.getChannelTrailerWord();
	       short lachadawin = adclist.at( adclist.size()-1 );
	       if ( (frailer>>12 != 0x5) || ( (frailer&0xfff) != tpc_channel_number ) ) {
	       std::vector<short> kazufix = decodeChannelTrailer( (unsigned short)lachadawin, (unsigned short)frailer );
	       for ( auto& it : kazufix )
	       adclist.emplace_back( it );
	       }
	    */
        // Output tree variables - for calculating compression
	    chdsize = adclist.size();
        if (crate_number == 1){
          ADCwords_crate1 += chdsize;
        }
        if (crate_number == 2){
          ADCwords_crate2 += chdsize;
        }
        if (crate_number == 3){
          ADCwords_crate3 += chdsize;
        }
        if (crate_number == 4){
          ADCwords_crate4 += chdsize;
        }
        if (crate_number == 5){
          ADCwords_crate5 += chdsize;
        }
        if (crate_number == 6){
          ADCwords_crate6 += chdsize;
        }
        if (crate_number == 7){
          ADCwords_crate7 += chdsize;
        }
        if (crate_number == 8){
          ADCwords_crate8 += chdsize;
        }
        if (crate_number == 9){
          ADCwords_crate9 += chdsize;
        }
        if (crate_number == 0){
          ADCwords_crate0 += chdsize; 
        }
	    const static size_t          fAdcList_size = chdsize;
	    if (fAdcList_size!=chdsize) {
	      throw art::Exception( art::errors::FileReadError ) 
		<< "Unexpected data: Number of words for channel: " 
		<< tpc_channel_number << " different than first waveform in the readout ("
		<< fAdcList_size << "!=" << chdsize << ") ... That's really bad!!!" << std::endl;
	    }
	      
	  /* else {
	    const ub_RawData& chD = channel.data(); 
	    // chdsize=(chD.getChannelDataSize()/sizeof(uint16_t));    
	    // chdsize = chD.size()/sizeof(uint16_t);    
	    chdsize = chD.size();
	    adclist.reserve(chD.size()); // optimize
	    for(ub_RawData::const_iterator it = chD.begin(); it!= chD.end(); it++) {
	    adclist.push_back(*it);
	    }
	    //              chD.decompress();
	    }*/
	    
	  //int crate_number = tpc_crate.crateHeader()->crate_number;
	  util::UBDaqID daqId( crate_number, card.getModule(), tpc_channel_number);

	  int ch=0;
	  auto it_chsearch = fChannelMap.find(daqId);
	  if ( it_chsearch!=fChannelMap.end() ){
	    ch=(*it_chsearch).second;
	    //              fChannelMap[daqId];
	    //              wire=fWireMap[daqId];
	    //              pl=fPlaneMap[daqId];
	  }
	  else {
	    if ( ( crate_number==1 && card.getModule()==8 && (tpc_channel_number>=32 && tpc_channel_number<64) ) ||
		 ( crate_number==9 && card.getModule()==5 && (tpc_channel_number>=32 && tpc_channel_number<64) ) ) {
	      // As of 6/22/2016: We expect these FEM channels to have no database entry.
	      continue; // do not write to data product
	    }
	    else {
	      // unexpected channels are missing. throw.
	      char warn[256];
	      sprintf( warn, "Warning DAQ ID not found ( %d, %d, %d )!", crate_number, card.getModule(), tpc_channel_number );
	      throw std::runtime_error( warn );
	    }
	  }
	  //\todo fix this once there is a proper channel table
	  // else{
	  //   //continue;
	  //   ch=10000*tpc_crate.crateHeader()->crate_number
	  //     +100*card.getModule()
	  //     +tpc_channel_number;
	  // }

	  //if (int(ch) >= 8254)
	  // continue;
	  //raw::Compress_t compression=raw::kHuffman;
	  //if (fHuffmanDecode) compression=raw::kNone;
	  raw::Compress_t compression=raw::kNone; // as of June 19,2015 compression not used by the DAQ. Data stored is uncompressed.
	  if ( adclist.size()!=9595 ) {
	    char warn[256];
	    sprintf( warn, "Error: Number of ADCs in (crate,slot,channel)=( %d, %d, %d ) does not equal 9595!", crate_number, card.getModule(), tpc_channel_number );
	    throw std::runtime_error( warn );
	  }
	  if (fSwizzleTPC){ // here is where we actually fill the output
  	    raw::RawDigit rd(ch,chdsize,adclist,compression);
  	    tpcDigitList.push_back(rd);
        }
	    
	  /*
            std::cout << ch << "\t"
	    << int(crate_header.getCrateNumber()) << "\t" 
	    << card_header.getModule() << "\t"
	    << channel_number << "\t"
	    << rms << std::endl;
	  */

	}//<--End channel_it for loop
      }//<---End card_it for loop
    }//<---End seb_it for loop

    if ( tpcDigitList.size()!=8256 ) {
      char warn[256];
      sprintf( warn, "Error: Number of channels saved (%d) did not match the expectation (8256)!", (int)tpcDigitList.size() );
      //throw std::runtime_error( warn );
    }

    mf::LogInfo("")<< "Got to end of fillTPCData().";
  }

  // =====================================================================
  void LArRawInputDriverUBooNE::fillPMTData(ubdaq::ub_EventRecord& event_record,
					    std::map< opdet::UBOpticalChannelCategory_t, std::unique_ptr<std::vector<raw::OpDetWaveform>> >& pmtDigitList )
  {
    //fill PMT data

    if (!fSwizzlePMT){
        std::cout << "Swizzling turned off for PMT Data so skipping that.." << std::endl;
        return;
    }

      
    // MODIFIED by Nathaniel Sat May 16, to use my new version of datatypes (v6_08, on branch master)
    
    //crate -> card -> channel -> window

    auto const* timeService = lar::providerFrom<detinfo::DetectorClocksService>();
    ::art::ServiceHandle<geo::UBOpReadoutMap> ub_pmt_channel_map;
    
    // pmt channel map is assumed to be time dependent. therefore we need event time to set correct map.
    ubdaq::ub_GlobalHeader global_header = event_record.getGlobalHeader();
    if(fUseGPS)
      global_header.useGPSTime();
    else if(fUseNTP)
      global_header.useLocalHostTime();
    uint32_t seconds=global_header.getSeconds();
    time_t mytime = (time_t)seconds;
    if ( mytime==0 ) {
      // some events seem o be missing time stamp. use run number in this case.
      std::cout << "[LArRawInputDriverUBooNE::fillPMTData] event epoch time 0 (!?). using run to set channel map" << std::endl;
      ub_pmt_channel_map->SetOpMapRun( global_header.getRunNumber() );
    }
    else
      ub_pmt_channel_map->SetOpMapTime( mytime );
    
    using namespace gov::fnal::uboone::datatypes;
    
    auto const seb_pmt_map = event_record.getPMTSEBMap();
    
    if (seb_pmt_map.empty()) {
        std::cerr << "Warning swizzler didn't find any PMT data in the event." << std::endl;
        std::cerr << "If this is a calibration or laser run, that's ok." << std::endl;
        fSwizzlePMT=false;
        return;
    }

    N_PMT_waveforms = 0;
    for(auto const& it:  seb_pmt_map) {
      pmt_crate_data_t const& crate_data = it.second;
      //      int crate_number = crate_data.crateHeader()->crate_number;
      
      if ( !crate_data.wasDissected() ) {
	std::cerr << "Warning PMT crate data corrupted! Skipping." << std::endl;
	continue;
      }

      //now get the card map (for the current crate), and do a loop over all cards
      std::vector<pmt_crate_data_t::card_t> const& cards = crate_data.getCards();
       
      for( pmt_crate_data_t::card_t const& card_data : cards ) {
        
	if ( !card_data.wasDissected() ) {
	  std::cerr << "Warning PMT card data corrupted! Skipping." << std::endl;
	  continue;
	}
            
// Frame and sample for trigger
        uint32_t frame = RollOver(card_data.getFrame(), card_data.getTrigFrameMod16(), 4);
        uint32_t sample = card_data.getTrigSample();

        if (PMTframe == -999){PMTframe = frame;} // internal frame consistency checking
        if (abs(frame - PMTframe) > 1){ // if the frame doesn't match the other PMT frames here then we have a problem
          std::cerr << "ERROR!" << std::endl;
          std::cerr << "PMT card header trigger frames not within one frame of each other!!" << std::endl;
          throw std::exception();
        }


//        // check if we have swizzled the trigger data
//        if (fSwizzleTrigger){
//          if (triggerFrame == -999){ // if we have swizzled trigger data, and we have the default triggerFrame, we have a problem
//            std::cerr << "ERROR!" << std::endl;
//            std::cerr << "Trigger data should have been swizzled, but triggerFrame remains as the default value of -999" << std::endl;
//            throw std::exception();
//          }
//          if (abs(frame - triggerFrame) > 1){ // if the triggerFrame and frame don't agree here then we have a problem
//            std::cerr << "ERROR!" << std::endl;
//            std::cerr << "PMT card header trigger frame not within one frame of trigger header frame!! Trigger header frame is " 
//                      << triggerFrame << " and PMT card header frame is " << frame << std::endl;
//            throw std::exception();
//          }
//        }
//         Filling output tree variables
        if (card_data.getModule() == 5){
          FEM5triggerFrame = frame;
          FEM5triggerSample = sample;
          FEM5triggerTime = timeService->OpticalClock().Time( sample, frame );
        }
        if (card_data.getModule() == 6){
          FEM6triggerFrame = frame;
          FEM6triggerSample = sample;
          FEM6triggerTime = timeService->OpticalClock().Time( sample, frame );
        }
	//        int card_number = card_data.getModule();
        
        // nathaniel's version of datatypes:
        for(auto const& channel_data : card_data.getChannels() ) { // auto here is pmt_crate_data_t::card_t::card_channel-type

	  // if ( !channel_data.wasDissected() ){
	  //   std::cerr << "Warning PMT channel data corrupted! Skipping." << std::endl;
	  //   continue;
	  // }

          int channel_number = channel_data.getChannelNumber();
          
          //now get the windows
          auto const& windows = channel_data.getWindows();  // auto here is std::vector<ub_PMT_WindowData_v6>
          for(const auto& window: windows ) {               // auto here is ub_PMT_WindowData_v6
            const auto& window_header = window.header();    // auto here is ub_PMT_WindowHeader_v6
            const ub_RawData& window_data = window.data();
            size_t win_data_size=window_data.size();

            
            // //\todo check category, time & frame
            // optdata::Optical_Category_t category = optdata::kUndefined;
            // if ((window_header.getDiscriminantor()&0x04)==0x04) {
            //   category=optdata::kBeamPMTTrigger;
            // } else {
            //   category=optdata::kCosmicPMTTrigger;
            // }
	    // tmw: In this new scheme, category is no longer needed (5/26/15)
            
            uint32_t sample=window_header.getSample();
	    uint32_t frame =RollOver(card_data.getFrame(),window_header.getFrame(),3);
	    //std::cout<<" FEM: " << card_data.getFrame() << " ... Channel: " << frame << " ... sample: " << sample << std::endl;
	    unsigned int data_product_ch_num = ub_pmt_channel_map->GetChannelNumberFromCrateSlotFEMCh( crate_data.crateHeader()->crate_number, card_data.getModule(), channel_number );
	    //int crate_number = crate_data.crateHeader()->crate_number; 
	    //std::cout << "fill (CSF): " << crate_number << ", " << card_data.getModule() << ", " << channel_number << " ==> Readout Channel " << data_product_ch_num << std::endl;
	    
	    // here we translate crate/card/daq channel to data product channel number
	    // also need to go from clock time to time stamp
	    opdet::UBOpticalChannelCategory_t ch_category = ub_pmt_channel_map->GetChannelCategory( data_product_ch_num );
	    double window_timestamp = timeService->OpticalClock().Time( sample, frame );

            raw::OpDetWaveform rd( window_timestamp, data_product_ch_num, win_data_size);
            rd.reserve(win_data_size); // Don't know if this compiles, but it is more efficient. push_back is terrible without it.

	    //std::cout << " into ReadoutCH=" << data_product_ch_num << " category=" << opdet::UBOpChannelEnumName( ch_category ) << std::endl;
	    short adc_max=0;
	    short adc_max_sample=0;
            for(ub_RawData::const_iterator it = window_data.begin(); it!= window_data.end(); it++){ 
              rd.push_back(*it & 0xfff);                
              if(adc_max < rd.back()){
                adc_max_sample = rd.size();
                adc_max = rd.back();
              }

            }
            // fill OpDetWaveform time
            double OpDetWaveForm_time = rd.TimeStamp();
            if (N_PMT_waveforms < 400 and channel_number < 36){
              PMT_waveform_times[N_PMT_waveforms] = OpDetWaveForm_time;
              N_PMT_waveforms += 1;
            }
            // Saving trigger readout stream variables to output file
            if(adc_max>2150) { // logic pulses have high ADC values

              if (channel_number == 39){
//                std::cout << "Found RWM signal!" << std::endl;
//                std::cout << "RWM signal at frame, sample " << RollOver(card_data.getFrame(),window_header.getFrame(),3) << ", " <<  window_header.getSample() + adc_max_sample << std::endl;
                RO_RWMtriggerFrame = RollOver(card_data.getFrame(),window_header.getFrame(),3);
                RO_RWMtriggerSample = window_header.getSample() + adc_max_sample;
                RO_RWMtriggerTime = timeService->OpticalClock().Time( RO_RWMtriggerSample, RO_RWMtriggerFrame);
//                std::cout << "window size = " << win_data_size << std::endl;
              }
              else if (channel_number == 38){
//                std::cout << "Found STROBE signal!" << std::endl;
//                std::cout << "STROBE signal at frame, sample " << RollOver(card_data.getFrame(),window_header.getFrame(),3) <<  ", " << window_header.getSample() + adc_max_sample << std::endl;
                RO_EXTtriggerFrame = RollOver(card_data.getFrame(),window_header.getFrame(),3);
                RO_EXTtriggerSample = window_header.getSample() + adc_max_sample;
                RO_EXTtriggerTime = timeService->OpticalClock().Time( RO_EXTtriggerSample, RO_EXTtriggerFrame);
//                std::cout << "window size = " << win_data_size << std::endl;
              }
              else if (channel_number == 37){
//                std::cout << "Found NuMI signal!" << std::endl;
//                std::cout << "NuMI signal at frame, sample " << RollOver(card_data.getFrame(),window_header.getFrame(),3) <<  ", " << window_header.getSample() << std::endl;
                RO_NuMItriggerFrame = RollOver(card_data.getFrame(),window_header.getFrame(),3);
                RO_NuMItriggerSample = window_header.getSample() + adc_max_sample;
                RO_NuMItriggerTime = timeService->OpticalClock().Time( RO_NuMItriggerSample, RO_NuMItriggerFrame);
//                std::cout << "window size = " << win_data_size << std::endl;
              }
              else if (channel_number == 36){
//                std::cout << "Found BNB signal!" << std::endl;
//                std::cout << "BNB signal at frame, sample " << RollOver(card_data.getFrame(),window_header.getFrame(),3) <<  ", " << window_header.getSample() + adc_max_sample << std::endl;
                RO_BNBtriggerFrame = RollOver(card_data.getFrame(),window_header.getFrame(),3);
                RO_BNBtriggerSample = window_header.getSample() + adc_max_sample;
                RO_BNBtriggerTime = timeService->OpticalClock().Time( RO_BNBtriggerSample, RO_BNBtriggerFrame);
//                std::cout << "window size = " << win_data_size << std::endl;
              }
//              else if (channel_number == 46){
//                std::cout << "Found GATE2 (NuMI) signal!" << std::endl;
//                std::cout << "NuMI gate at frame, sample " << RollOver(card_data.getFrame(),window_header.getFrame(),3) <<  ", " << window_header.getSample() << std::endl;
//                RO_Gate2Frame = RollOver(card_data.getFrame(),window_header.getFrame(),3);
//                RO_Gate2Sample = window_header.getSample();
//                std::cout << "window size = " << win_data_size << std::endl;
//              }
//              else if (channel_number == 47){
//                std::cout << "Found GATE1 (BNB) signal!" << std::endl;
//                std::cout << "BNB gate at frame, sample " << RollOver(card_data.getFrame(),window_header.getFrame(),3) <<  ", " << window_header.getSample() << std::endl;
//                RO_Gate1Frame = RollOver(card_data.getFrame(),window_header.getFrame(),3);
//                RO_Gate1Sample = window_header.getSample();
//                std::cout << "window size = " << win_data_size << std::endl;
//              }
//              else {
//                std::cout << "channel " << channel_number << ",regular PMT window size = " << win_data_size << std::endl;
//              }
            }
            if (fSwizzlePMT){
              pmtDigitList[ch_category]->emplace_back(rd);
            }
          }
        }//<--End channel_pmt_it for loop
      }//<---End card_pmt_it for loop
    }//<---End seb_pmt_it for loop
    
  }

  // =====================================================================
  void LArRawInputDriverUBooNE::fillBeamData(ubdaq::ub_EventRecord& event_record,
                                             raw::BeamInfo& beamInfo)
  {
    /*
      ubdaq::ub_BeamHeader bh=event_record.getBeamHeader();
      std::vector<ubdaq::ub_BeamData> bdv=event_record.getBeamDataVector();
      if (bdv.size()>0) {
      beamInfo.SetRecordType(bh.getRecordType());
      beamInfo.SetSeconds(bh.getSeconds());
      beamInfo.SetMilliSeconds(bh.getMilliSeconds());
      beamInfo.SetNumberOfDevices(bh.getNumberOfDevices());
      
      for (int i=0;i<bh.getNumberOfDevices();i++) {
      beamInfo.Set(bdv[i].getDeviceName(),bdv[i].getData());
      if (fHistMapBeam.find(bdv[i].getDeviceName())!=fHistMapBeam.end()) 
      fHistMapBeam[bdv[i].getDeviceName()]->Fill(bdv[i].getData()[0]);
      }
      }
    */
  }

  // =====================================================================
  void LArRawInputDriverUBooNE::fillTriggerData(gov::fnal::uboone::datatypes::ub_EventRecord &event_record,
						std::vector<raw::Trigger>& trigInfo)
  {

    auto const* timeService = lar::providerFrom<detinfo::DetectorClocksService>();

    for(auto const& it_trig_map : event_record.getTRIGSEBMap()){

      //int seb_num = it_trig_map.first;
      auto const& trig_crate = it_trig_map.second;  //  is typedef of ub_Trigger_CrateData_X
      auto const& trig_card  = trig_crate.getTriggerCardData(); // typedef of ub_Trigger_CardData_X
      auto const& trig_header = trig_crate.getTriggerHeader();   // ub_Trigger_HeaderData_X
      auto const& trig_data   = trig_crate.getTriggerData();     // ub_Trigger_
      
      // The following is to skip events if we don't care about that trigger type.  It does not currently work.
      if (fSwizzleTriggerType == "BNB" and trig_data.Trig_Gate2() == 0){
        skipEvent = true;
        std::cout << "skipping non BNB" << std::endl;
        return;
      } 
      else if (fSwizzleTriggerType == "NuMI" and not trig_data.Trig_Gate1()){
        skipEvent = true;
        std::cout << "skipping non NuMI" << std::endl;
        return;
      } 
      else if (fSwizzleTriggerType == "EXT" and not trig_data.Trig_EXT()){
        skipEvent = true;
        std::cout << "skipping non EXT" << std::endl;
        return;
      } 
      else if (fSwizzleTriggerType == "CALIB" and not trig_data.Trig_Calib()){
        skipEvent = true;
        std::cout << "skipping non CALIB" << std::endl;
        return;
      } 

      // Make a trigger clock 
      unsigned int sample_64MHz = (trig_header.get2MHzSampleNumber() * 32) + (trig_header.get16MHzRemainderNumber() * 4) + trig_data.getPhase();
      unsigned int frame = trig_header.getFrame();
      //std::cout << "Trigger frame: " << frame << " ... sample : " << sample_64MHz << std::endl;
      detinfo::ElecClock trig_clock = timeService->OpticalClock( sample_64MHz, frame);

      double trigger_time = trig_clock.Time();
      double beam_time = -1;
      if ( trig_data.Trig_Gate1() || trig_data.Trig_Gate2() ) // 1) NUMI : 2) BNB
	beam_time = trigger_time;
      uint32_t trig_bits = trig_data.getPMTTrigData();
      if( trig_data.Trig_PC()       ) trig_bits += ( 0x1 << ::trigger::kTriggerPC    ); 
      if( trig_data.Trig_EXT()      ) trig_bits += ( 0x1 << ::trigger::kTriggerEXT   );
      if( trig_data.Trig_Active()   ) trig_bits += ( 0x1 << ::trigger::kActive       );
      if( trig_data.Trig_Gate1()    ) trig_bits += ( 0x1 << ::trigger::kTriggerNuMI  );
      if( trig_data.Trig_Gate2()    ) trig_bits += ( 0x1 << ::trigger::kTriggerBNB   );
      if( trig_data.Trig_Veto()     ) trig_bits += ( 0x1 << ::trigger::kVeto         );
      if( trig_data.Trig_Calib()    ) trig_bits += ( 0x1 << ::trigger::kTriggerCalib );
      if( trig_data.Trig_GateFake() ) trig_bits += ( 0x1 << ::trigger::kFakeGate     );
      if( trig_data.Trig_BeamFake() ) trig_bits += ( 0x1 << ::trigger::kFakeBeam     );
      if( trig_data.Trig_Spare1()   ) trig_bits += ( 0x1 << ::trigger::kSpare        );
	  
      raw::Trigger swiz_trig( trig_card.getTrigNumber(),
			      trigger_time,
			      beam_time,
			      trig_bits );
      if (not kazuTestSwizzleTrigger){return;}
      trigInfo.emplace_back( swiz_trig );
      
      if (not fSwizzleTrigger){return;} // if we don't want to swizzle the trigger data, then stop here.
// variables saving to output tree
      triggerFrame = frame;
      triggerSample = sample_64MHz;
      triggerActive = trig_data.Trig_Active();
      triggerBitBNB = trig_data.Trig_Gate2(); 
      triggerBitNuMI = trig_data.Trig_Gate1(); 
      triggerBitEXT = trig_data.Trig_EXT(); 
      triggerBitPMTBeam = trig_bits & 0x1;
      triggerBitPMTCosmic = trig_bits & 0x2;
      triggerTime = trigger_time;

    }
  }

  // =====================================================================  
  void LArRawInputDriverUBooNE::fillSWTriggerData(gov::fnal::uboone::datatypes::ub_EventRecord &event_record,
						raw::ubdaqSoftwareTriggerData& trigInfo)
  {
    auto const* timeService = lar::providerFrom<detinfo::DetectorClocksService>();
    
    std::vector<ub_FEMBeamTriggerOutput> swTrig_vect;
    try {

      // Software trigger data pulled from DAQ software trigger
      swTrig_vect = event_record.getSWTriggerOutputVector();
    }
    catch(...){ // softwrare trigger data product not present in binary file (because it's too old probably).  Just set values to a default
      std::cout << "failed to obtain software trigger object from binary file - setting all values to default" << std::endl;
      return;
    }

    if (swTrig_vect.empty()){
        std::cout << "The SWTriggerOutputVector was empty, so setting all values to the default." << std::endl;
        return;
    }
    //
    // Figure out time w.r.t. Trigger - dirty but works.
    //
    auto const& opt_clock = timeService->OpticalClock();
    //auto const& trg_clock = timeService->TriggerClock();
    //auto const& tpc_clock = timeService->TPCClock();
    auto const& trig_map  = event_record.getTRIGSEBMap();
    auto const& trig_header = trig_map.begin()->second.getTriggerHeader();
    auto const& trig_data = trig_map.begin()->second.getTriggerCardData().getTriggerData();
    uint64_t trig_sample_number = trig_header.get2MHzSampleNumber() * 32;
    trig_sample_number += trig_header.get16MHzRemainderNumber() * 4;
    trig_sample_number += trig_data.getPhase();
    double trigger_time = (double)(trig_header.getFrame() * opt_clock.FramePeriod());
    trigger_time += ((double)trig_sample_number) * opt_clock.TickPeriod();

    uint64_t trig_tick = trig_sample_number + trig_header.getFrame() * opt_clock.FrameTicks();
     
    //NOTE that if there's no PMT Data, there isn't gonna be any SW Trigger information.
    //need to make this more complicated if there are ever triggers based on things other than PMTS

    auto const seb_pmt_map = event_record.getPMTSEBMap();
    if (seb_pmt_map.empty()) {
        std::cout << "The PMTSEBMap was empty and there was no PMT data, so setting all SWTrigger values to the default." << std::endl;
        return;
    }

    auto const& crate_data = event_record.getPMTSEBMap().begin()->second;
     
    //auto const& pmt_card    = crate_data.getCards().at(0);
    
    uint64_t beam_ro_tick = 0;
    auto const& card_data = crate_data.getCards().front();
    uint64_t min_dt = 1e12; //FIXME this should be set to max integer value from compiler
    // First search the target timing
    for(auto const& ch_data : card_data.getChannels()){
  
      for(auto const& window : ch_data.getWindows()) {
     
        if(window.header().getDiscriminantor()!=ub_PMT_DiscriminatorTypes_v6::BEAM && 
           window.header().getDiscriminantor()!=ub_PMT_DiscriminatorTypes_v6::BEAM_GATE){
          continue; //ignore non-BEAM signals
        }
      
        uint64_t window_time = RollOver(card_data.getFrame(), window.header().getFrame(), 3) * 102400;
        window_time += window.header().getSample();
        uint64_t window_trigger_dt = 
          ( window_time < trig_tick ? trig_tick - window_time : window_time - trig_tick );
        
        if( min_dt > window_trigger_dt ) {
          min_dt       = window_trigger_dt;
          beam_ro_tick = window_time;
        }
      }
    }

    if(beam_ro_tick > trig_tick){
      _trigger_beam_window_time = beam_ro_tick - trig_tick;
    }
    else{
      _trigger_beam_window_time = trig_tick - beam_ro_tick;
      _trigger_beam_window_time *= -1.;
    }
    
    for (unsigned int i(0); i < swTrig_vect.size(); ++i){ // loop through swtrigger algos filling info
      ub_FEMBeamTriggerOutput swTrig = swTrig_vect.at(i); // fetch algorithm
      
      trigInfo.addAlgorithm(swTrig.algo_instance_name, // add algorithm to art data product
                            swTrig.pass_algo, 
                            swTrig.pass_prescale, 
                            swTrig.amplitude, 
                            swTrig.multiplicity, 
                            swTrig.time, 
                            _trigger_beam_window_time + swTrig.time * opt_clock.TickPeriod(), 
                            swTrig.prescale_weight);
    } // end loop over swtrigger algos
  }
    

  // =====================================================================  
  std::vector<short> LArRawInputDriverUBooNE::decodeChannelTrailer(unsigned short last_adc, unsigned short data)
  {
    // bug fix for missing channel trailer in TPC Data.
    // undoes the hack that fixed the above where the last word is used as the trailer
    // we then use the fake trailer, or frailer, combine it with the last word in the channel data window, or lachadawin, 
    // to recover the end of the channel waveform.

    //std::vector<unsigned short> res;
    std::vector<short> res;
    if(data>>12 == 0x0) {
      //std::cout << "Non-huffman data word..." << std::endl;
      res.push_back(  (short) data & 0xfff);
      return res;
    }
    if(data>>14 == 0x2) {
      //std::cout << "Huffman data word..." << std::endl;
      size_t zero_count=0;
      for(int index=13; index>=0; --index) {
	if(!(data>>index & 0x1)) zero_count +=1;
	else {
	  switch(zero_count){
	      
	  case 0:
	    break;
	  case 1:
	    last_adc -= 1; break;
	  case 2:
	    last_adc += 1; break;
	  case 3:
	    last_adc -= 2; break;
	  case 4:
	    last_adc += 2; break;
	  case 5:
	    last_adc -= 3; break;
	  case 6:
	    last_adc += 3; break;
	  default:

	    std::cerr << "Unexpected 0-count for huffman word: "
		      << "\033[95m"
		      << zero_count << " zeros in the word 0x"
		      << std::hex
		      << data
		      << std::dec
		      << "\033[00m"
		      << std::endl;
	    std::cerr << "Binary representation of the whole word: "
		      << "\033[00m";
	    for(int i=15; i>=0; --i)
	      std::cout << ((data>>i) & 0x1);
	    std::cout << "\033[00m" << std::endl;
	    throw std::exception();
	  }
	  res.push_back((short)last_adc);
	  zero_count = 0;
	}
      }
      return res;
    }

    std::cerr << "\033[93mERROR\033[00m Unexpected upper 4 bit: 0x"
	      << std::hex
	      << ((data >> 12) & 0xf)
	      << std::dec
	      << std::endl;
    throw std::exception();
  }


  // =====================================================================  
  void LArRawInputDriverUBooNE::checkTimeStampConsistency(){

    //Note that fSwizzlePMT will be set to false if there is no PMT data in the event.
    //at the end of the event, it is reset to the initial value from the fhicl parameters.
    //but this means that the trig-TPC comparison will continue, but not the TPC-PMT or Trig-PMT
    if (fSwizzleTrigger && fSwizzlePMT ){ // trig-PMT comparison
      if (abs(PMTframe - triggerFrame)>1){
        std::cout << "ERROR!" << std::endl;
        std::cout << "trigger data and PMT data both read out, but frames disagree (by more than 1)!" << std::endl;
        std::cout << "trigger data frame = " << triggerFrame << ", PMT data frame = " << PMTframe << std::endl;
        throw std::exception();
      } // Done trig-PMT comparison
    }
    if (fSwizzleTrigger && fSwizzleTPC ){ // trig-TPC comparison
      if (abs(TPCframe - triggerFrame)>1){
        std::cout << "ERROR!" << std::endl;
        std::cout << "trigger data and TPC data both read out, but frames disagree (by more than 1)!" << std::endl;
        std::cout << "trigger data frame = " << triggerFrame << ", TPC data frame = " << TPCframe << std::endl;
        throw std::exception();
      } // Done trig-TPC comparison
    }
    if (fSwizzleTPC && fSwizzlePMT ){ // TPC-PMT comparison
      if (abs(PMTframe - TPCframe)>1){
        std::cout << "ERROR!" << std::endl;
        std::cout << "TPC data and PMT data both read out, but frames disagree (by more than 1)!" << std::endl;
        std::cout << "TPC data frame = " << TPCframe << ", PMT data frame = " << PMTframe << std::endl;
        throw std::exception();
      }
    } // Done TPC-PMT comparison

  }
  // =====================================================================  

}//<---Endlris
