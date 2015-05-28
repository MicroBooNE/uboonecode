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
#include "RawData/RawDigit.h"
#include "RawData/TriggerData.h"
#include "RawData/DAQHeader.h"
#include "RawData/BeamInfo.h"
#include "RawData/OpDetWaveform.h"
#include "Geometry/Geometry.h"
#include "Geometry/ExptGeoHelperInterface.h"
#include "SummaryData/RunData.h"
#include "Utilities/TimeService.h"
#include "OpticalDetectorData/OpticalTypes.h" // lardata -- I want to move the enums we use back to UBooNE as they are UBooNE-specific

//ART, ...
#include "art/Framework/IO/Sources/put_product_in_principal.h"
#include "art/Utilities/Exception.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//uboone datatypes

// uboonecode
#include "uboone/Geometry/ChannelMapUBooNEAlg.h"
#include "uboone/Geometry/UBOpChannelTypes.h"

#include "datatypes/raw_data_access.h"




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
    fEventCounter(0),
    fHuffmanDecode(ps.get<bool>("huffmanDecode",false))
  {
    helper.reconstitutes<raw::DAQHeader,              art::InEvent>("daq");
    helper.reconstitutes<std::vector<raw::RawDigit>,  art::InEvent>("daq");
    helper.reconstitutes<raw::BeamInfo,               art::InEvent>("daq");
    registerOpticalData( helper ); //helper.reconstitutes<std::vector<raw::OpDetWaveform>,art::InEvent>("daq");
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
    mf::LogInfo(__FUNCTION__)<<"File boundary (processed "<<fEventCounter<<" events)"<<std::endl;
    fCurrentSubRunID.flushSubRun();
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

    return;

    //seek to the end of file, check the end word, check number of events and sizes
    //read in the end of file word first
    uint16_t end_of_file_marker;
    fInputStream.seekg( -1*sizeof(uint16_t) , std::ios::end);
    fInputStream.read( (char*)&end_of_file_marker , sizeof(uint16_t));

    if(end_of_file_marker != 0xe0f0){
      //throw art::Exception( art::errors::FileReadError )
      std::cout << "File "<<name<<" has incorrect end of file marker. "<< end_of_file_marker<<std::endl;
    }
    fInputStream.seekg(std::ios::beg);
    
    return;
    /*
    //get number of events from word at end of file
    fInputStream.seekg( -1*(sizeof(uint16_t)+sizeof(uint32_t)), std::ios::end);
    fInputStream.read( (char*)&fNumberOfEvents , sizeof(uint32_t));

    fNumberOfEvents = 10;
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
    */
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
  void LArRawInputDriverUBooNE::registerOpticalData( art::ProductRegistryHelper &helper ) {
    // we make a data product for each category of channels
    fPMTdataProductNames.clear();
    for ( unsigned int cat=0; cat<(unsigned int)opdet::NumUBOpticalChannelCategories; cat++ ) {
      std::stringstream ss;
      ss << "pmtreadout" << opdet::UBOpChannelEnumName( (opdet::UBOpticalChannelCategory_t)cat );
      helper.reconstitutes<std::vector<raw::OpDetWaveform>,art::InEvent>(ss.str()); 
      fPMTdataProductNames.insert( std::make_pair( (opdet::UBOpticalChannelCategory_t)cat, ss.str() ) );
    }
  }

  // =====================================================================
  void LArRawInputDriverUBooNE::putPMTDigitsIntoEvent( std::map< opdet::UBOpticalChannelCategory_t, std::unique_ptr< std::vector<raw::OpDetWaveform> > >& pmtdigitlist, 
						       art::EventPrincipal* &outE ) {
    for ( unsigned int cat=0; cat<(unsigned int)opdet::NumUBOpticalChannelCategories; cat++ ) {
      
      art::put_product_in_principal(std::move( pmtdigitlist[(opdet::UBOpticalChannelCategory_t)cat]  ),
				    *outE,
				    fPMTdataProductNames[ (opdet::UBOpticalChannelCategory_t)cat ]); // Module label
    }
    
  }
  
  // =====================================================================
  bool LArRawInputDriverUBooNE::readNext(art::RunPrincipal* const &/*inR*/,
                                         art::SubRunPrincipal* const &/*inSR*/,
                                         art::RunPrincipal* &outR,
                                         art::SubRunPrincipal* &outSR,
                                         art::EventPrincipal* &outE)
  {
    mf::LogInfo(__FUNCTION__)<<"Attempting to read event: "<<fEventCounter<<std::endl;
    // Create empty result, then fill it from current file:
    std::unique_ptr<raw::DAQHeader> daq_header(new raw::DAQHeader);
    std::unique_ptr<std::vector<raw::RawDigit> >  tpc_raw_digits( new std::vector<raw::RawDigit>  );
    std::unique_ptr<raw::BeamInfo> beam_info(new raw::BeamInfo);
    std::unique_ptr<raw::Trigger> trig_info(new raw::Trigger);
    std::map< opdet::UBOpticalChannelCategory_t, std::unique_ptr< std::vector<raw::OpDetWaveform> > > pmt_raw_digits;
    for ( unsigned int opdetcat=0; opdetcat<(unsigned int)opdet::NumUBOpticalChannelCategories; opdetcat++ ) {
      pmt_raw_digits.insert( std::make_pair( (opdet::UBOpticalChannelCategory_t)opdetcat, std::unique_ptr< std::vector<raw::OpDetWaveform> >(  new std::vector<raw::OpDetWaveform> ) ) );
    }

    bool res=false;

    res=processNextEvent(*tpc_raw_digits, pmt_raw_digits, *daq_header, *beam_info, *trig_info );

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
	       <<" event: " << daq_header->GetEvent()
	       <<"\033[00m"<< std::endl;
      */
      outE = fSourceHelper.makeEventPrincipal(fCurrentSubRunID.run(),
					      fCurrentSubRunID.subRun(),
					      daq_header->GetEvent(),
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
      putPMTDigitsIntoEvent( pmt_raw_digits, outE );
      // art::put_product_in_principal(std::move(pmt_raw_digits),
      //                               *outE,
      //                               "daq"); // Module label
     
    }

    return res;

  }

  // =====================================================================
  bool LArRawInputDriverUBooNE::processNextEvent(std::vector<raw::RawDigit>& tpcDigitList,
                                                 std::map< opdet::UBOpticalChannelCategory_t, std::unique_ptr<std::vector<raw::OpDetWaveform>> >& pmtDigitList,
                                                 raw::DAQHeader& daqHeader,
                                                 raw::BeamInfo& beamInfo,
						 raw::Trigger& trigInfo)
  {       
    //try {
      boost::archive::binary_iarchive ia(fInputStream); 
      ubdaq::ub_EventRecord event_record;  
      ia >> event_record;
      //std::cout<<event_record.debugInfo()<<std::endl;
      //set granularity 
      //      event_record.updateIOMode(ubdaq::IO_GRANULARITY_CHANNEL);
      
      fillDAQHeaderData(event_record, daqHeader);
      fillTPCData(event_record, tpcDigitList);
      fillPMTData(event_record, pmtDigitList);
      fillBeamData(event_record, beamInfo);
      fillTriggerData(event_record, trigInfo);
      //std::cout<<"Done ProcessNextEvent..."<<std::endl;
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
      //        ubdaq::crateHeader crate_header = seb_it->first;
      //        ubdaq::crateData crate_data = seb_it->second;
      //      int tpc_seb_num = seb_it.first;
      tpc_crate_data_t const& tpc_crate = seb_it.second;


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

      //      auto const& tpc_crate_header = tpc_crate.header();    
      //      auto const& tpc_crate_trailer = tpc_crate.trailer();  

      //Special to the crate, there is a special header that the DAQ attaches. You can access this
      //like so. The type here is a unique ptr to a ub_CrateHeader_v6 struct. That has useful info
      //like the local host time, which may or may not be set properly right now...
      auto const& tpc_crate_DAQ_header = tpc_crate.crateHeader(); // I think auto should be tpc_crate_data_t::ub_CrateHeader_t --NJT
      //     ub_LocalHostTime this_time = tpc_crate_DAQ_header->local_host_time;
      
      //The Crate Data is split up into Cards. You use the "getCards()" command to get access to
      //each of those. Note that calling this function will dissect the data if it has not already
      //been dissected (debugInfo() calls getCards()). You can do a look over the cards like so:
      for(auto const& card : tpc_crate.getCards()){  // This auto is tpc_crate_data_t::card_t

        //The format here is similar to the crate! There's a header (which is a ub_TPC_CardHeader_v*
        //object), and technically a trailer (though here it's empty!).
	//	auto const& tpc_card_header = card.header();   
	//	auto const& tpc_card_trailer = card.trailer(); 

        //Of course, you can probe for information in the card header. You'll have to find the appropriate
        //header file to know what type you have, but again, these will follow typical practice. And, you
        //can always use debugInfo to not only print the info, but it tells you the type.
        // auto const this_event_number = card.getEvent(); /// auto are ints here
        // auto const this_frame_number = card.getFrame(); /// auto are ints here


        //And, you guessed it, the tpc card data is split up into one more level: by channel.
        for(auto const& channel : card.getChannels()){ // auto here tpc_crate_data_t::card_t::card_channel_type

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
	    if (fHuffmanDecode) {
              channel.decompress(adclist); // All-in-one call.
            } else {
              const ub_RawData& chD = channel.data(); 
	      //	      chdsize=(chD.getChannelDataSize()/sizeof(uint16_t));    
	      chdsize = chD.size()/sizeof(uint16_t);    
              adclist.reserve(chD.size()); // optimize
              for(ub_RawData::const_iterator it = chD.begin(); it!= chD.end(); it++) {
                adclist.push_back(*it);
              }
	      //              chD.decompress();
            }

            daqid_t daqId(tpc_crate.crateHeader()->crate_number,
                          card.getModule(),
                          tpc_channel_number);

            int ch=-1;
            if (fChannelMap.find(daqId)!=fChannelMap.end()){
              ch=fChannelMap[daqId];
              //              wire=fWireMap[daqId];
              //              pl=fPlaneMap[daqId];
            }
            //\todo fix this once there is a proper channel table
            else{
              //continue;
              ch=10000*tpc_crate.crateHeader()->crate_number
                +100*card.getModule()
                +tpc_channel_number;
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
					    std::map< opdet::UBOpticalChannelCategory_t, std::unique_ptr<std::vector<raw::OpDetWaveform>> >& pmtDigitList )
  {
    //fill PMT data

      // MODIFIED by Nathaniel Sat May 16, to use my new version of datatypes (v6_08, on branch master)
    
    //crate -> card -> channel -> window

    ::art::ServiceHandle<geo::Geometry> geom;
    ::art::ServiceHandle< util::TimeService > timeService;
    std::shared_ptr<const geo::ChannelMapUBooNEAlg> ub_pmt_channel_map = std::dynamic_pointer_cast< const geo::ChannelMapUBooNEAlg >( geom->GetChannelMapAlg() ); // With UB-specific methods we need
    
    using namespace gov::fnal::uboone::datatypes;
    
    auto const seb_pmt_map = event_record.getPMTSEBMap();
    
    for(auto const& it:  seb_pmt_map) {
      pmt_crate_data_t const& crate_data = it.second;
      //      int crate_number = crate_data.crateHeader()->crate_number;
      
      //now get the card map (for the current crate), and do a loop over all cards
      std::vector<pmt_crate_data_t::card_t> const& cards = crate_data.getCards();
       
      for( pmt_crate_data_t::card_t const& card_data : cards ) {
        
	//        int card_number = card_data.getModule();
        
        // nathaniel's version of datatypes:
        for(auto const& channel_data : card_data.getChannels() ) { // auto here is pmt_crate_data_t::card_t::card_channel-type
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
            
            optdata::TimeSlice_t time=window_header.getSample();
            optdata::Frame_t frame=window_header.getFrame();
	    //int crate_number = crate_data.crateHeader()->crate_number; 
	    //std::cout << "fill (CSF): " << crate_number << ", " << card_data.getModule() << ", " << channel_number;
	    
	    // here we translate crate/card/daq channel to data product channel number
	    // also need to go from clock time to time stamp
	    unsigned int data_product_ch_num = ub_pmt_channel_map->GetChannelNumberFromCrateSlotFEMCh( crate_data.crateHeader()->crate_number, card_data.getModule(), channel_number );
	    opdet::UBOpticalChannelCategory_t ch_category = ub_pmt_channel_map->GetChannelType( data_product_ch_num );
	    double window_timestamp = timeService->OpticalClock().Time( time, frame );
            raw::OpDetWaveform rd( window_timestamp, channel_number,win_data_size);
            rd.reserve(win_data_size); // Don't know if this compiles, but it is more efficient. push_back is terrible without it.

	    //std::cout << " into ReadoutCH=" << data_product_ch_num << " category=" << opdet::UBOpChannelEnumName( ch_category ) << std::endl;
	    
	    
            for(ub_RawData::const_iterator it = window_data.begin(); it!= window_data.end(); it++){ 
              rd.push_back(*it & 0xfff);                
            }
            pmtDigitList[ch_category]->emplace_back(rd);
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

  void LArRawInputDriverUBooNE::fillTriggerData(gov::fnal::uboone::datatypes::ub_EventRecord &event_record,
						raw::Trigger& trigInfo)
  {

    /*
      auto const& trig_data = event_record.triggerData();
    trigInfo = raw::Trigger( trig_data.getTrigEventNum(),
                              
    Trigger(unsigned int counter,
            double       trigger_time,
            double       beamgate_time,
            uint32_t     bits)
      : fTriggerNumber       ( counter           ),
      fTriggerTime         ( trigger_time      ),
      fBeamGateTime        ( beamgate_time     ),
      fTriggerBits         ( bits              )
      {}



    uint16_t getSampleNumber() const noexcept;
    uint16_t getSampleRemainder() const noexcept;
    uint16_t getSampleNumber_2MHz() const noexcept;
    uint16_t getSampleNumber_16MHz() const noexcept;
    uint16_t getSampleNumber_64MHz() const noexcept;
    bool     getBusy() const noexcept;
    uint32_t getFrame() const noexcept;
    uint32_t getTrigEventNum() const noexcept;
    uint16_t  getTriggerBits() const noexcept;
    bool     isPmtTrigger() const noexcept;
    bool     isExtTrigger() const noexcept;
    bool     isActiveTrigger() const noexcept;
    bool     isBnbTrigger()    const noexcept;
    bool     isNumiTrigger()   const noexcept;
    bool     isVetoTrigger()   const noexcept;
    bool     isCalibTrigger()  const noexcept;
    uint16_t getPhase64Mhz_1() const noexcept;
    uint16_t getPhase64Mhz_2() const noexcept;
    */
  }

}//<---Endlris

