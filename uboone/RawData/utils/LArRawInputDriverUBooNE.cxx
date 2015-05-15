///////////////////////////////////////////////////////////////////////
/// \file  LArRawInputDriverUBooNE.cxx
/// \brief Source to convert raw binary files to root files
/// \Original Authors
/// \version $Id: LArRawInputDriver.h,v 1.7 2010/01/14 19:20:33 brebel Exp $
/// \author  brebel@fnal.gov, soderber@fnal.gov
/// \MicroBooNE Author: jasaadi@fnal.gov, zsrko@fnal.gov (with much help from Wes and Eric)
////////////////////////////////////////////////////////////////////////

//LArSoft 
#include "uboone/RawData/utils/LArRawInputDriverUBooNE.h"
#include "RawData/RawDigit.h"
#include "RawData/DAQHeader.h"
#include "RawData/BeamInfo.h"
#include "OpticalDetectorData/FIFOChannel.h"
#include "Geometry/Geometry.h"
#include "SummaryData/RunData.h"

//ART, ...
#include "art/Framework/IO/Sources/put_product_in_principal.h"
#include "art/Utilities/Exception.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//uboone datatypes


/*
#include "datatypes/eventRecord.h"
#include "datatypes/globalHeader.h"
#include "datatypes/triggerData.h"
#include "datatypes/crateHeader.h"
#include "datatypes/crateData.h"
#include "datatypes/cardHeader.h"
#include "datatypes/cardData.h"
#include "datatypes/channelData.h"
#include "datatypes/crateDataPMT.h"
#include "datatypes/cardHeaderPMT.h"
#include "datatypes/cardDataPMT.h"
#include "datatypes/channelDataPMT.h"
#include "datatypes/windowHeaderPMT.h"
#include "datatypes/windowDataPMT.h"
*/
#include "datatypes/ub_BeamHeader.h"
#include "datatypes/ub_BeamData.h"



//boost
//#include <boost/archive/binary_iarchive.hpp>
#include <boost/algorithm/string.hpp>

//root
#include "TH1D.h"

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

  // ======================================================================
  LArRawInputDriverUBooNE::LArRawInputDriverUBooNE(fhicl::ParameterSet const & ps, 
						   art::ProductRegistryHelper &helper,
						   art::SourceHelper const &pm)
    :
    fSourceHelper(pm),
    fCurrentSubRunID(),
    fNumberOfEvents(-1),
    fEventCounter(0),
    fHuffmanDecode(ps.get<bool>("huffmanDecode",false))
  {
    helper.reconstitutes<raw::DAQHeader,              art::InEvent>("daq");
    helper.reconstitutes<std::vector<raw::RawDigit>,  art::InEvent>("daq");
    helper.reconstitutes<std::vector<optdata::FIFOChannel>,art::InEvent>("daq");
    helper.reconstitutes<raw::BeamInfo,               art::InEvent>("daq");
    initChannelMap();

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
  }
  
  
  // ======================================================================
  void LArRawInputDriverUBooNE::closeCurrentFile()
  {
    fCurrentSubRunID.flushSubRun();
    fNumberOfEvents=-1;
    fEventCounter=0;
    fChannelMap.clear();
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
    uint16_t end_of_file_marker;
    fInputStream.seekg( -1*sizeof(uint16_t) , std::ios::end);
    fInputStream.read( (char*)&end_of_file_marker , sizeof(uint16_t));

    if(end_of_file_marker != 0xe0f0){
      throw art::Exception( art::errors::FileReadError )
	<< "File "<<name<<" has incorrect end of file marker. "<< end_of_file_marker<<std::endl;
    }
    
    //get number of events from word at end of file
    fInputStream.seekg( -1*(sizeof(uint16_t)+sizeof(uint32_t)), std::ios::end);
    fInputStream.read( (char*)&fNumberOfEvents , sizeof(uint32_t));

    //now get all of the event sizes, 
    uint32_t tmp_event_size;
    std::streampos count_event_size=0;
    fInputStream.seekg( -1*(sizeof(uint16_t)+(fNumberOfEvents+1)*sizeof(uint32_t)), std::ios::end);
    for(uint32_t i=0; i<fNumberOfEvents; i++){
      // since we want the beginning, push back the event size before incrementing it
      fEventLocation.push_back(count_event_size);
      fInputStream.read( (char*)&tmp_event_size , sizeof(uint32_t));
      count_event_size += tmp_event_size;
    }
    fEventLocation.push_back(count_event_size);
    fInputStream.seekg(std::ios::beg);

    mf::LogInfo("")<<"Opened file "<<name<<" with "<<fNumberOfEvents<<" event(s)";
  }

  // ======================================================================  
  void LArRawInputDriverUBooNE::initChannelMap()
  {

    fChannelMap.clear();
    mf::LogInfo("")<<"Fetching channel map from DB";

    PGconn *conn = PQconnectdb("host=fnalpgsdev.fnal.gov port=5436 dbname=uboonedaq_dev user=uboonedaq_web password=argon!uBooNE");	
    
    if(PQstatus(conn)!=CONNECTION_OK) {
      mf::LogError("") << "Couldn't open connection to postgresql interface";
      PQfinish(conn);
      throw art::Exception( art::errors::FileReadError )
	<< "Failed to get channel map from DB." << std::endl;
    }

    PGresult *res  = PQexec(conn, "BEGIN");    
    if (PQresultStatus(res) != PGRES_COMMAND_OK) { 
      mf::LogError("")<< "postgresql BEGIN failed";
      PQclear(res);
      PQfinish(conn);
      throw art::Exception( art::errors::FileReadError )
	<< "postgresql BEGIN failed." << std::endl;
    }

    // Andrzej's changed database

    PQclear(res);
    res = PQexec(conn,
		 "SELECT crate_id, slot, wireplane, larsoft_channel, channel_id "
		 " FROM channels NATURAL JOIN asics NATURAL JOIN motherboards NATURAL JOIN coldcables NATURAL JOIN motherboard_mapping NATURAL JOIN intermediateamplifiers_copy NATURAL JOIN servicecables NATURAL JOIN servicecards NATURAL JOIN warmcables_copy NATURAL JOIN ADCreceivers_copy_new NATURAL JOIN crates NATURAL JOIN fecards"
		 );

    // Standard & wrong Hoot database
    /*
    PQclear(res);
    res = PQexec(conn,
		 "SELECT crate_id, slot, wireplane, larsoft_channel, channel_id "
		 " FROM channels NATURAL JOIN asics NATURAL JOIN motherboards NATURAL JOIN coldcables NATURAL JOIN motherboard_mapping NATURAL JOIN intermediateamplifiers_copy NATURAL JOIN servicecables NATURAL JOIN servicecards NATURAL JOIN warmcables NATURAL JOIN ADCreceivers NATURAL JOIN crates NATURAL JOIN fecards"
		 );
    */
    if ((!res) || (PQresultStatus(res) != PGRES_TUPLES_OK))
    {
      mf::LogError("")<< "SELECT command did not return tuples properly";
      PQclear(res);
      PQfinish(conn);
      throw art::Exception( art::errors::FileReadError )
	<< "postgresql SELECT failed." << std::endl;
    }

    int num_records=PQntuples(res);
    for (int i=0;i<num_records;i++) {
      int crate_id     = atoi(PQgetvalue(res, i, 0));
      int slot         = atoi(PQgetvalue(res, i, 1));
      //auto const wPl   =      PQgetvalue(res, i, 2);
      int larsoft_chan = atoi(PQgetvalue(res, i, 3));
      int channel_id   = atoi(PQgetvalue(res, i, 4));

      int boardChan = channel_id%64;


      if (crate_id==9 && slot==5)
	boardChan = (channel_id-32)%64;

      /*
      std::cout << "Looking up in DB: [Crate, Card, Channel]: [" << crate_id << ", "
		<< slot << ", " << boardChan << "]";
      std::cout << "\tCh. Id (LArSoft): " << larsoft_chan  << std::endl;
      */
      daqid_t daq_id(crate_id,slot,boardChan);
      std::pair<daqid_t, int> p(daq_id,larsoft_chan);

      if (fChannelMap.find(daq_id) != fChannelMap.end()){
	std::cout << "Multiple entries!" << std::endl;
	mf::LogWarning("")<<"Multiple DB entries for same (crate,card,channel). "<<std::endl
		       <<"Redefining (crate,card,channel)=>id link ("
		       <<daq_id.crate<<", "<<daq_id.card<<", "<<daq_id.channel<<")=>"
		       <<fChannelMap[daq_id];
      }
      fChannelMap.insert(p);
    }
  }

  // =====================================================================
  bool LArRawInputDriverUBooNE::readNext(art::RunPrincipal* const &/*inR*/,
					 art::SubRunPrincipal* const &/*inSR*/,
					 art::RunPrincipal* &outR,
					 art::SubRunPrincipal* &outSR,
					 art::EventPrincipal* &outE)
  {
    // Create empty result, then fill it from current file:
    std::unique_ptr<raw::DAQHeader> daq_header(new raw::DAQHeader);
    std::unique_ptr<std::vector<raw::RawDigit> >  tpc_raw_digits( new std::vector<raw::RawDigit>  );
    std::unique_ptr<std::vector<optdata::FIFOChannel> >  pmt_raw_digits( new std::vector<optdata::FIFOChannel>  );
    std::unique_ptr<raw::BeamInfo> beam_info(new raw::BeamInfo);

    bool res=false;

    if (fEventCounter < fNumberOfEvents) 
      res=processNextEvent(*tpc_raw_digits, *pmt_raw_digits, *daq_header, *beam_info );
    
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

      outE = fSourceHelper.makeEventPrincipal(fCurrentSubRunID.run(),
						fCurrentSubRunID.subRun(),
						daq_header->GetEvent(),
						tstamp);
    

      // Put products in the event.
      art::put_product_in_principal(std::move(tpc_raw_digits),
				    *outE,
				    "daq"); // Module label
      art::put_product_in_principal(std::move(pmt_raw_digits),
				    *outE,
				    "daq"); // Module label
      art::put_product_in_principal(std::move(daq_header),
				    *outE,
				    "daq"); // Module label
      art::put_product_in_principal(std::move(beam_info),
				    *outE,
				    "daq"); // Module label
     
    }

    return res;

  }

  // =====================================================================
  bool LArRawInputDriverUBooNE::processNextEvent(std::vector<raw::RawDigit>& tpcDigitList,
						 std::vector<optdata::FIFOChannel>& pmtDigitList,
						 raw::DAQHeader& daqHeader,
						 raw::BeamInfo& beamInfo)
  {       
    fInputStream.seekg(fEventLocation[fEventCounter]);
    try {
      boost::archive::binary_iarchive ia(fInputStream); 
      ubdaq::ub_EventRecord event_record;  
      ia >> event_record;
      //set granularity 
      //      event_record.updateIOMode(ubdaq::IO_GRANULARITY_CHANNEL);

      fillDAQHeaderData(event_record, daqHeader);
      fillTPCData(event_record, tpcDigitList);
      fillPMTData(event_record, pmtDigitList);
      fillBeamData(event_record, beamInfo);

    } catch ( art::Exception const& e ) {
      throw art::Exception( art::errors::FileReadError )
	<< "Failed to read the event."<<e.what() << std::endl;
    }
  
    return true;
  }
  
  // =====================================================================
  void LArRawInputDriverUBooNE::fillDAQHeaderData(ubdaq::ub_EventRecord& event_record,
						  raw::DAQHeader& daqHeader)
  {
      ubdaq::ub_GlobalHeader global_header = event_record.getGlobalHeader();
      
      // art::Timestamp is an unsigned long long. The conventional 
      // use is for the upper 32 bits to have the seconds since 1970 epoch 
      // and the lower 32 bits to be the number of nanoseconds within the 
      // current second.
      // (time_t is a 64 bit word)

      uint32_t seconds=global_header.getSeconds();
      uint32_t nano_seconds=global_header.getNanoSeconds();
      time_t mytime = ((time_t)seconds<<32) | nano_seconds;

      //\/      uint32_t subrun_num = global_header->getSubrunNumber();
      
      daqHeader.SetStatus(1);
      daqHeader.SetFileFormat(global_header.getRecordType());
      daqHeader.SetSoftwareVersion(global_header.DAQ_version_number);
      daqHeader.SetRun(global_header.getRunNumber());
      daqHeader.SetSubRun(global_header.getSubrunNumber());
      
      //\/ Add the subRun number too!
      daqHeader.SetEvent(global_header.getEventNumber());
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
      // ### Swizzling to get the number of channels...trying the method used in write_read.cpp
      // ### provided by Wes --- About the data:
      // ### The format of the data is in levels: crate, card, channel.
      // ### Level 1: The event record contains a map of (crateHeader,crateData) pairs.
      // ### Level 2: Each crateData object may contain a map of (cardHeader,cardData) pairs.
      // ### Level 3: Each cardData object may contain a map of (int,channelData) pairs. The int
      // ### is the channel number.
      
      //get the seb map, and do a loop over all sebs/crates

    for( auto const& seb_it : event_record.getTPCSEBMap()) {
	
	//get the crateHeader/crateData objects
      //	ubdaq::crateHeader crate_header = seb_it->first;
      //	ubdaq::crateData crate_data = seb_it->second;

	// Get Time information:
	//uint32_t sebTSec = crate_header.getSebTimeSec();
	//std::cout << "Seb Time (sec) : " << sebTSec << std::endl;
	//crate_header_t crHeader = crate_header.getCrateHeader();
	// GPStime in UNIX second/micro/nano info
	//gps_time_t GPStime = crHeader.gps_time;
	// DAQtime is time of last update of GPS time (in frame, sample, div)
	//tbclkub_t  DAQtime = crHeader.daqClock_time;
	//std::cout << "GPS Time seconds: " << GPStime.second << std::endl;
	//std::cout << "DAQ Frame: " << DAQtime.frame << "\tSample: " << DAQtime.sample << std::endl;

      auto const& tpc_crate_header = tpc_crate.header();
      auto const& tpc_crate_trailer = tpc_crate.trailer();

      //Special to the crate, there is a special header that the DAQ attaches. You can access this
      //like so. The type here is a unique ptr to a ub_CrateHeader_v6 struct. That has useful info
      //like the local host time, which may or may not be set properly right now...
      auto const& tpc_crate_DAQ_header = tpc_crate.crateHeader();
      ub_LocalHostTime this_time = tpc_crate_DAQ_header->local_host_time;
      
      //The Crate Data is split up into Cards. You use the "getCards()" command to get access to
      //each of those. Note that calling this function will dissect the data if it has not already
      //been dissected (debugInfo() calls getCards()). You can do a look over the cards like so:
      for(auto const& card : tpc_crate.getCards()){

	//The format here is similar to the crate! There's a header (which is a ub_TPC_CardHeader_v*
	//object), and technically a trailer (though here it's empty!).
	auto const& tpc_card_header = card.header();
	auto const& tpc_card_trailer = card.trailer();

	//Of course, you can probe for information in the card header. You'll have to find the appropriate
	//header file to know what type you have, but again, these will follow typical practice. And, you
	//can always use debugInfo to not only print the info, but it tells you the type.
	auto const this_event_number = card.getEvent();
	auto const this_frame_number = card.getFrame();


	//And, you guessed it, the tpc card data is split up into one more level: by channel.
	for(auto const& channel : card.getChannels()){

	  //There's a header and trailer here. Remember these are just uint16_t, that contain the
	  //channel number.
	  auto const& tpc_channel_header = channel.header();
	  auto const& tpc_channel_trailer = channel.trailer();

	  //The channel object (ub_MarkedRawChannelData) has a method for returning the channel.
	  //You can look at the other objects too (like ub_MarkedRawCardData) and see methods of
	  //use there as well.
	  auto const tpc_channel_number = channel.getChannelNumber();
	    
	  ubdaq::channelData chD = channel.data;

	    //Huffman decoding
	  if (fHuffmanDecode) chD.decompress();

	    //\todo here fill in the detector RawData structures. something like this:
	    //\tode following code is from pulse_viewer.cpp, not sure what is the most 
	    //      elegant way to do this
	    //\todo adclist is a vector<short>, but need only 12bit so probably 
	    //      ok to cast uint16_t to short

	    std::unique_ptr<uint16_t> blk_chD(new uint16_t);
	    std::vector<short> adclist;
	    size_t chdsize=(chD.getChannelDataSize()/sizeof(uint16_t));
	    char* cdptr=chD.getChannelDataPtr();

	    for (unsigned int i=0;i<chdsize;i++) {
	      std::copy(cdptr+i*sizeof(uint16_t),
			cdptr+(i+1)*sizeof(uint16_t),
			(char*)blk_chD.get());

	      adclist.push_back(*blk_chD);	      
	    }

	    daqid_t daqId(crate_header.getCrateNumber(),
			  card_header.getModule(),
			  channel_number);

	    int ch=-1;
	    if (fChannelMap.find(daqId)!=fChannelMap.end()){
	      ch=fChannelMap[daqId];
	      //	      wire=fWireMap[daqId];
	      //	      pl=fPlaneMap[daqId];
	    }
	    //\todo fix this once there is a proper channel table
	    else{
	      //continue;
	      ch=10000*crate_header.getCrateNumber()
		+100*card_header.getModule()
		+channel_number;
	    }
	    //if (int(ch) >= 8254)
	    // continue;
	    raw::Compress_t compression=raw::kHuffman;
	    if (fHuffmanDecode) compression=raw::kNone;

	    raw::RawDigit rd(ch,chdsize,adclist,compression);

	    /*
	    std::cout << ch << "\t"
		      << int(crate_header.getCrateNumber()) << "\t" 
		      << card_header.getModule() << "\t"
		      << channel_number << "\t"
		      << rms << std::endl;
	    */

	    tpcDigitList.push_back(rd);
	  }//<--End channel_it for loop
	}//<---End card_it for loop
      }//<---End seb_it for loop
  }

  // =====================================================================
  void LArRawInputDriverUBooNE::fillPMTData(ubdaq::ub_EventRecord& event_record,
					    std::vector<optdata::FIFOChannel>& pmtDigitList)
  {
    //fill PMT data

    //crate -> card -> channel -> window
    std::map<ubdaq::crateHeader,ubdaq::crateDataPMT,ubdaq::compareCrateHeader> seb_pmt_map = event_record.getSEBPMTMap();
    std::map<ubdaq::crateHeader,ubdaq::crateDataPMT>::iterator seb_pmt_it;
    for( seb_pmt_it = seb_pmt_map.begin(); seb_pmt_it != seb_pmt_map.end(); seb_pmt_it++){
      //get the crateHeader/crateData objects
      //crateHeader crate_header = seb_pmt_it->first;
      ubdaq::crateDataPMT crate_pmt_data = seb_pmt_it->second;
      
      //now get the card map (for the current crate), and do a loop over all cards
      std::map<ubdaq::cardHeaderPMT,ubdaq::cardDataPMT>::iterator card_pmt_it;
      std::map<ubdaq::cardHeaderPMT,ubdaq::cardDataPMT,ubdaq::compareCardHeaderPMT> card_pmt_map = crate_pmt_data.getCardMap();
      for(card_pmt_it = card_pmt_map.begin(); card_pmt_it != card_pmt_map.end(); card_pmt_it++){
	//get the cardHeader/cardData objects
	//cardHeaderPMT card_pmt_header = card_pmt_it->first;
	ubdaq::cardDataPMT card_pmt_data = card_pmt_it->second;
	
	//now get the channel map (for the current card), and do a loop over all channels
	std::map<int,ubdaq::channelDataPMT> channel_pmt_map = card_pmt_data.getChannelMap();
	std::map<int,ubdaq::channelDataPMT>::iterator channel_pmt_it;
	for(channel_pmt_it = channel_pmt_map.begin(); channel_pmt_it != channel_pmt_map.end(); channel_pmt_it++){
	  //get the channel number and channelData
	  int channel_number = channel_pmt_it->first;
	  ubdaq::channelDataPMT channel_pmt_data = channel_pmt_it->second;
	  
	  //now get the windows
	  std::map<ubdaq::windowHeaderPMT,ubdaq::windowDataPMT>::iterator window_it;
	  std::map<ubdaq::windowHeaderPMT,ubdaq::windowDataPMT,ubdaq::compareWindowHeaderPMT> window_map = channel_pmt_data.getWindowMap();
	  for(window_it = window_map.begin(); window_it != window_map.end(); window_it++){
	    ubdaq::windowHeaderPMT winHeader=window_it->first;
	    ubdaq::windowDataPMT winData=window_it->second;
	    size_t win_data_size=winData.getWindowDataSize()/sizeof(uint16_t);
	    
	    //\todo check category, time & frame
	    optdata::Optical_Category_t category = optdata::kUndefined;
	    switch (winHeader.getDiscriminant()) {
	    case 1:
	      category=optdata::kCosmicPMTTrigger;
	      break;
	    case 3:
	      category=optdata::kBeamPMTTrigger;
	      break;
	    default:
	      category=optdata::kUndefined;
	     
	    }
	    optdata::TimeSlice_t time=winHeader.getSample();
	    optdata::Frame_t frame=winHeader.getFrame();
	    optdata::FIFOChannel rd(category, time, frame, channel_number,win_data_size);

	    std::unique_ptr<uint16_t> blk_winD(new uint16_t);	    
	    for (unsigned int i=0;i<win_data_size;i++) {
	      std::copy(winData.getWindowDataPtr()+i*sizeof(uint16_t),
			winData.getWindowDataPtr()+(i+1)*sizeof(uint16_t),
			(char*)blk_winD.get());
	      
	      rd.push_back(*blk_winD & 0xfff);
	    }
	    pmtDigitList.push_back(rd);
	  }
	}//<--End channel_pmt_it for loop
      }//<---End card_pmt_it for loop
    }//<---End seb_pmt_it for loop
    
  }

  // =====================================================================
  void LArRawInputDriverUBooNE::fillBeamData(ubdaq::ub_EventRecord& event_record,
					     raw::BeamInfo& beamInfo)
  {
    ubdaq::beamHeader bh=event_record.getBeamHeader();
    std::vector<ubdaq::beamData> bdv=event_record.getBeamDataVector();
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
  }
}//<---Endlris

