////////////////////////////////////////////////////////////////////////
/// \file  LArRawInputDriverUBooNE.cxx
/// \brief Source to convert raw binary files to root files
/// \Original Authors
/// \version $Id: LArRawInputDriver.h,v 1.7 2010/01/14 19:20:33 brebel Exp $
/// \author  brebel@fnal.gov, soderber@fnal.gov
/// \MicroBooNE Author: jasaadi@fnal.gov (with much help from Wes and Eric)
////////////////////////////////////////////////////////////////////////

#include "RawData/utils/LArRawInputDriverUBooNE.h"

#include "RawData/RawDigit.h"
#include "RawData/DAQHeader.h"
#include "RawData/BeamInfo.h"
#include "RawData/OpDetPulse.h"
#include "Geometry/Geometry.h"
#include "SummaryData/RunData.h"

#include "art/Framework/IO/Sources/put_product_in_principal.h"
#include "art/Utilities/Exception.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

//uboone datatypes
#include "RawData/uboone_datatypes/beamData.h"
#include "RawData/uboone_datatypes/beamHeader.h"
#include "RawData/uboone_datatypes/cardData.h"
#include "RawData/uboone_datatypes/cardHeader.h"
#include "RawData/uboone_datatypes/channelData.h"
#include "RawData/uboone_datatypes/crateData.h"
#include "RawData/uboone_datatypes/crateHeader.h"
#include "RawData/uboone_datatypes/eventHeaderTrailer.h"
#include "RawData/uboone_datatypes/eventRecord.h"
#include "RawData/uboone_datatypes/globalHeader.h"
#include "RawData/uboone_datatypes/gps.h"
#include "RawData/uboone_datatypes/triggerData.h"

extern "C" {
#include <sys/types.h>
#include <iostream>
#include <dirent.h>
#include <cstdlib>
#include <libpq-fe.h>
}
#include <bitset>

using namespace gov::fnal::uboone::datatypes;

// ======================================================================
// UBooNE data interface, adapted from code by Rebel/Soderberg:
//  modified J. Asaadi March 20th, 2013
//  modified Z. Pavlovic September 23rd, 2013  
// ======================================================================

// +++ Blatant stealing from LongBo
namespace lris {
  // ======================================================================
  // class c'tor/d'tor:
  LArRawInputDriverUBooNE::LArRawInputDriverUBooNE(fhicl::ParameterSet const & ps, 
						   art::ProductRegistryHelper &helper,
						   art::PrincipalMaker const &pm)
    :
    fPrincipalMaker(pm),
    fCurrentSubRunID(),
    fNumberOfEvents(-1),
    fEventCounter(0),
    fHuffmanDecode(ps.get<bool>("huffmanDecode",false))
  {
    helper.reconstitutes<raw::DAQHeader,              art::InEvent>("daq");
    helper.reconstitutes<std::vector<raw::RawDigit>,  art::InEvent>("daq");
    helper.reconstitutes<std::vector<raw::OpDetPulse>,art::InEvent>("daq");
    helper.reconstitutes<raw::BeamInfo,               art::InEvent>("daq");
    initChannelMap();
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
    //\todo? perhaps could use these locations when reading in the file
    
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

    PQclear(res);
    res = PQexec(conn,
		 "SELECT crate_id, slot, wireplane, wirenum, channel_id "
		 " FROM channels NATURAL JOIN asics NATURAL JOIN motherboards NATURAL JOIN coldcables NATURAL JOIN motherboard_mapping NATURAL JOIN intermediateamplifiers NATURAL JOIN servicecables NATURAL JOIN servicecards NATURAL JOIN warmcables NATURAL JOIN ADCreceivers NATURAL JOIN crates NATURAL JOIN fecards "
		 );

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
      int crate_id = atoi(PQgetvalue(res, i, 0));
      int slot     = atoi(PQgetvalue(res, i, 1));
      //char wireplane = *PQgetvalue(res, i, 2);
      int wirenum  = *PQgetvalue(res, i, 3);
      int channel_id = atoi(PQgetvalue(res, i, 4));

      daqid_t daq_id(crate_id,slot,wirenum);
      std::pair<daqid_t, int> p(daq_id,channel_id);

      if (fChannelMap.find(daq_id) != fChannelMap.end())
	mf::LogInfo("")<<"Multiple DB entries for same (crate,card,channel)"<<std::endl
		       <<"Redefining (crate,card,channel)=>id link ("
		       <<daq_id.crate<<", "<<daq_id.card<<", "<<daq_id.channel<<")=>"
		       <<fChannelMap[daq_id];
    
      fChannelMap.insert(p);
    }

    /*
    std::stringstream ss;
    ss<<"Acquired channel mapping from DB"<<std::endl;
    ss<<"     (  crate,   card, channel) => ID"<<std::endl;
    for (auto it=fChannelMap.begin();it!=fChannelMap.end();it++) {
      ss << "     ("
	 << std::setw(7) << it->first.crate  <<","
	 << std::setw(7) << it->first.card   <<", "
	 << std::setw(7) << it->first.channel<<") => "<< it->second<<std::endl;
    }    
    mf::LogInfo("")<<ss.str();    
    */
  }

  // =====================================================================
  bool LArRawInputDriverUBooNE::readNext(art::RunPrincipal* const &inR,
					 art::SubRunPrincipal* const &inSR,
					 art::RunPrincipal* &outR,
					 art::SubRunPrincipal* &outSR,
					 art::EventPrincipal* &outE)
  {
    // Create empty result, then fill it from current filename:
    std::unique_ptr<raw::DAQHeader> daq_header(new raw::DAQHeader);
    std::unique_ptr<std::vector<raw::RawDigit> >  raw_digits( new std::vector<raw::RawDigit>  );
    std::unique_ptr<std::vector<raw::OpDetPulse> >  pmt_raw_digits( new std::vector<raw::OpDetPulse>  );
    std::unique_ptr<raw::BeamInfo> beam_info(new raw::BeamInfo);

    bool res=false;

    if (fEventCounter < fNumberOfEvents) 
      res=processNextEvent(*raw_digits, *pmt_raw_digits, *daq_header, *beam_info );
    
    if (res) {
      fEventCounter++;
      art::RunNumber_t rn = daq_header->GetRun();
      art::Timestamp tstamp = daq_header->GetTimeStamp();
      art::SubRunID newID(rn, 1); //daq_header has no subrun information
      
      if (fCurrentSubRunID.runID() != newID.runID()) { // New Run
	outR = fPrincipalMaker.makeRunPrincipal(rn, tstamp);
      }
      if (fCurrentSubRunID != newID) { // New SubRun
	outSR = fPrincipalMaker.makeSubRunPrincipal(rn,
						    1,
						    tstamp);
	fCurrentSubRunID = newID;	
      }

      outE = fPrincipalMaker.makeEventPrincipal(fCurrentSubRunID.run(),
						fCurrentSubRunID.subRun(),
						daq_header->GetEvent(),
						tstamp);
    

      // Put products in the event.
      art::put_product_in_principal(std::move(raw_digits),
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
  bool LArRawInputDriverUBooNE::processNextEvent(std::vector<raw::RawDigit>& digitList,
						 std::vector<raw::OpDetPulse>& pmtDigitList,
						 raw::DAQHeader& daqHeader,
						 raw::BeamInfo& beamInfo)
  {    
    
    fInputStream.seekg(fEventLocation[fEventCounter]);
    try {
      boost::archive::binary_iarchive ia(fInputStream); 
      eventRecord event_record;  // Declare an eventRecord object. This is yours.
      ia >> event_record;

      //      mf::LogInfo("") <<"Processing event "<<fEventCounter;
      //set granularity 
      event_record.updateIOMode(IO_GRANULARITY_CHANNEL);

      //Here's how to get info from the global header in the event record.
      //Other header-type objects will be very, very similar.
      globalHeader *global_header = event_record.getGlobalHeaderPtr();
      
      // ########################################################
      // ### Back to swizzling: using global_header-> methods ###
      // ########################################################
      time_t mytime = global_header->getSeconds();
      
      // I suspect this is antiquated...but might need to confirm?
      
      mytime = mytime << 32;//Nov. 2, 2010 - "time_t" is a 64-bit word on many 64-bit machines
      //so we had to change types in header struct to read in the correct
      //number of bits.  Once we have the 32-bit timestamp from the binary
      //data, shift it up to the upper half of the 64-bit timestamp.  - Mitch
      
      daqHeader.SetStatus(1);
      daqHeader.SetFileFormat(global_header->getRecordType());
      daqHeader.SetSoftwareVersion(global_header->DAQ_version_number);
      daqHeader.SetRun(global_header->getRunNumber());
      daqHeader.SetEvent(global_header->getEventNumber());
      daqHeader.SetTimeStamp(mytime);

      /// \todo: What is the "fixed word" ? Leaving it unset for now
      /// \todo: What is the "spare word" ? Leaving it unset for now
      //daqHeader.SetFixedWord(h1.header);
      //daqHeader.SetSpareWord(h1.spare);
      
      // ### Swizzling to get the number of channels...trying the method used in write_read.cpp
      // ### provided by Wes --- About the data:
      // ### The format of the data is in levels: crate, card, channel.
      // ### Level 1: The event record contains a map of (crateHeader,crateData) pairs.
      // ### Level 2: Each crateData object may contain a map of (cardHeader,cardData) pairs.
      // ### Level 3: Each cardData object may contain a map of (int,channelData) pairs. The int
      // ### is the channel number.
      
      //get the seb map, and do a loop over all sebs/crates
      std::map<crateHeader,crateData,compareCrateHeader> seb_map = event_record.getSEBMap();
      std::map<crateHeader,crateData>::iterator seb_it;
      for( seb_it = seb_map.begin(); seb_it != seb_map.end(); seb_it++){
	
	//get the crateHeader/crateData objects
	crateHeader crate_header = seb_it->first;
	crateData crate_data = seb_it->second;

	//now get the card map (for the current crate), and do a loop over all cards
	std::map<cardHeader,cardData>::iterator card_it;
	std::map<cardHeader,cardData,compareCardHeader> card_map = crate_data.getCardMap();
	for(card_it = card_map.begin(); card_it != card_map.end(); card_it++){
	  
	  //get the cardHeader/cardData objects
	  cardHeader card_header = card_it->first;
	  cardData card_data = card_it->second;
	  
	  //now get the channel map (for the current card), and do a loop over all channels
	  std::map<int,channelData> channel_map = card_data.getChannelMap();
	  std::map<int,channelData>::iterator channel_it;
	  for(channel_it = channel_map.begin(); channel_it != channel_map.end(); channel_it++){
	    
	    //get the channel number and channelData
	    int channel_number = channel_it->first;
	    channelData chD = channel_it->second;

	    //Huffman decoding
	    if (fHuffmanDecode) chD.decompress();

	    //\todo here fill in the detector RawData structures. something like this:
	    //\tode following code is from pulse_viewer.cpp, not sure what is the most 
	    //      elegant way to do this
	    //\todo adclist is a vector<short> can't really go from uint16_t to short
	    std::unique_ptr<uint16_t> blk_chD(new uint16_t);
	    std::vector<short> adclist;

	    size_t chdsize=(chD.getChannelDataSize()/sizeof(uint16_t));
	    char* cdptr=chD.getChannelDataPtr();
	    
	    for (uint i=0;i<chdsize;i++) {
	      std::copy(cdptr+i*sizeof(uint16_t),
			cdptr+(i+1)*sizeof(uint16_t),
			(char*)blk_chD.get());
	    
	      adclist.push_back(*blk_chD);
	    }
	    
	    daqid_t daqId(crate_header.getCrateNumber(),
			  card_header.getModule(),
			  channel_number);

	    int ch=-1;
	    if (fChannelMap.find(daqId)!=fChannelMap.end()) ch=fChannelMap[daqId];
	    //\todo fix this once there is a proper channel table
	    else ch=10000*crate_header.getCrateNumber()
	             +100*card_header.getModule()
	                 +channel_number;
	    raw::Compress_t compression=raw::kHuffman;
	    if (fHuffmanDecode) compression=raw::kNone;

	    raw::RawDigit rd(ch,chdsize,adclist,compression);
	    digitList.push_back(rd);
	  }//<--End channel_it for loop
	}//<---End card_it for loop
      }//<---End seb_it for loop
    
      //fill beam data
      beamHeader bh=event_record.getBeamHeader();
      std::vector<beamData> bdv=event_record.getBeamDataVector();
      
      if (bdv.size()>0) {
	beamInfo.SetRecordType(bh.getRecordType());
	beamInfo.SetSeconds(bh.getSeconds());
	beamInfo.SetMilliSeconds(bh.getMilliSeconds());
	beamInfo.SetNumberOfDevices(bh.getNumberOfDevices());
      
	for (int i=0;i<bh.getNumberOfDevices();i++) 
	  beamInfo.Set(bdv[i].getDeviceName(),bdv[i].getData());
      }

      //fill PMT data
      std::map<crateHeader,crateDataPMT,compareCrateHeader> seb_pmt_map = event_record.getSEBPMTMap();
      std::map<crateHeader,crateDataPMT>::iterator seb_pmt_it;
      for( seb_pmt_it = seb_pmt_map.begin(); seb_pmt_it != seb_pmt_map.end(); seb_pmt_it++){
	//get the crateHeader/crateData objects
	//crateHeader crate_header = seb_pmt_it->first;
	crateDataPMT crate_pmt_data = seb_pmt_it->second;
	
	//now get the card map (for the current crate), and do a loop over all cards
	std::map<cardHeaderPMT,cardDataPMT>::iterator card_pmt_it;
	std::map<cardHeaderPMT,cardDataPMT,compareCardHeaderPMT> card_pmt_map = crate_pmt_data.getCardMap();
	for(card_pmt_it = card_pmt_map.begin(); card_pmt_it != card_pmt_map.end(); card_pmt_it++){
	  //get the cardHeader/cardData objects
	  //cardHeaderPMT card_pmt_header = card_pmt_it->first;
	  cardDataPMT card_pmt_data = card_pmt_it->second;
	  
	  //now get the channel map (for the current card), and do a loop over all channels
	  std::map<int,channelDataPMT> channel_pmt_map = card_pmt_data.getChannelMap();
	  std::map<int,channelDataPMT>::iterator channel_pmt_it;
	  for(channel_pmt_it = channel_pmt_map.begin(); channel_pmt_it != channel_pmt_map.end(); channel_pmt_it++){
	    //get the channel number and channelData
	    int channel_number = channel_pmt_it->first;
	    channelDataPMT channel_pmt_data = channel_pmt_it->second;

	    //now get the windows
	    std::map<windowHeaderPMT,windowDataPMT>::iterator window_it;
	    std::map<windowHeaderPMT,windowDataPMT,compareWindowHeaderPMT> window_map = channel_pmt_data.getWindowMap();
	    for(window_it = window_map.begin(); window_it != window_map.end(); window_it++){
	      windowHeaderPMT winHeader=window_it->first;
	      windowDataPMT winData=window_it->second;

	      //\todo fix the uint16_t vs short 
	      std::unique_ptr<uint16_t> blk_winD(new uint16_t);
	      std::vector<short> adclist;
	      
	      size_t win_data_size=winData.getWindowDataSize()/sizeof(uint16_t);
	      for (uint i=0;i<win_data_size;i++) {
		std::copy(winData.getWindowDataPtr()+i*sizeof(uint16_t),
			  winData.getWindowDataPtr()+(i+1)*sizeof(uint16_t),
			  (char*)blk_winD.get());

		adclist.push_back(*blk_winD);
	      }
	     
	      unsigned int pmt_frame=winHeader.getFrame();
	      unsigned int first_sample=0;
	      raw::OpDetPulse rd(channel_number,adclist,pmt_frame,first_sample);
	      pmtDigitList.push_back(rd);
	    }
	  }//<--End channel_pmt_it for loop
	}//<---End card_pmt_it for loop
      }//<---End seb_pmt_it for loop
    
    } catch ( ... ) {
      throw art::Exception( art::errors::FileReadError )
	<< "Failed to read the event." << std::endl;
    }
  
    return true;

  }
}//<---Endlris

