#include "beamRun.h"
#include "httpResponse.h"
#include "beamDAQConfig.h"
#include "MWRData.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>

#include <messagefacility/MessageLogger/MessageLogger.h>

#include "boost/date_time/posix_time/posix_time.hpp"
#include <boost/algorithm/string.hpp>
#include "boost/date_time/gregorian/gregorian.hpp"
#include <boost/date_time/c_local_time_adjustor.hpp>

using namespace gov::fnal::uboone::beam;
using namespace gov::fnal::uboone::datatypes;
using namespace std;
using namespace boost::posix_time;

beamRun::beamRun()
{
  //get config
  beamDAQConfig* bdconfig=beamDAQConfig::GetInstance();

  //set parameters
  fOutputDirData = bdconfig->GetDataOutputDir();
  fOutputDirInfo = bdconfig->GetInfoOutputDir();
  fBeamLine      = bdconfig->GetBeamLineList();
  fBundle        = bdconfig->GetBundles();
  fIFDBURL       = bdconfig->GetIFDBURL();
  fIFDBLatency   = bdconfig->GetIFDBLatency();
  fEventTypeMap  = bdconfig->GetEventTypeMap();
  fTimeWindowMap = bdconfig->GetTimeWindowMap();
  fMWRTimeOffset = bdconfig->GetMWRTimeOffset();

  fZoneOffset=second_clock::local_time()-second_clock::universal_time();
}

beamRun::~beamRun()
{
  for (unsigned int i=0;i<fBeamLine.size();i++) {
    fBundle[fBeamLine[i]].clear();
  }
  fBundle.clear();
  fBeamLine.clear();
  fEventTypeMap.clear();
}

void beamRun::StartRun(beamRunHeader& rh, boost::posix_time::ptime tstart) 
{
  fRunHeader.fRun=rh.fRun;
  fRunHeader.fSubRun=rh.fSubRun;
  fRunHeader.fRunStart=tstart;

  mf::LogInfo("") <<"Starting new beam run "<<fRunHeader.fRun<<" subrun "<<fRunHeader.fSubRun<<" on "<<tstart;
  fOut=new ofstream[fBeamLine.size()];
  for (unsigned int i=0;i<fBeamLine.size();i++) {
    stringstream ss;
    ss<<fOutputDirData<<"/beam_"<<fBeamLine[i]<<"_"
      <<setfill('0')<<setw(7)<<fRunHeader.fRun<<"_"
      <<setfill('0')<<setw(5)<<fRunHeader.fSubRun
      <<".dat";
    mf::LogInfo("") <<"Writing beam "<<fBeamLine[i]<<" data output to "
		    <<ss.str();
    fOut[i].open(ss.str().c_str(), ios::out | ios::binary);
    // fOA[fBeamLine[i]]=new boost::archive::binary_oarchive(fOut[i]);
    fLastQueryTime[fBeamLine[i]]=tstart;
  }
}

void beamRun::Update(boost::posix_time::ptime tend)
{
  static char sbuf[1024];
 
  for (unsigned int ibeam=0;ibeam<fBeamLine.size();ibeam++) {
    bool isdone=false;
    beamdatamap_t data_map;
    ptime last_proc_time=tend-hours(10000);
    //now get data from db 10min at the time
    ptime t1 = (fLastQueryTime[fBeamLine[ibeam]]+minutes(10)<tend) ? 
	(fLastQueryTime[fBeamLine[ibeam]]+minutes(10)) : tend;

    while ( !isdone ) {
      string endtime=to_iso_extended_string(t1);
      if (microsec_clock::local_time()-t1<milliseconds(fTimeWindowMap[fBeamLine[ibeam]])) endtime="now";

      httpResponse* response=new httpResponse();
      for ( unsigned int i=0;i<fBundle[fBeamLine[ibeam]].size();i++ ) {
	if ( (fBundle[fBeamLine[ibeam]][i].find("MWR") == std::string::npos) &&
	     (fBundle[fBeamLine[ibeam]][i].find("MultiWire") == std::string::npos) ){
	  sprintf(sbuf, "%s/data?b=%s&t0=%s&t1=%s&f=csv", fIFDBURL.c_str(), fBundle[fBeamLine[ibeam]][i].c_str(),
		  to_iso_extended_string(fLastQueryTime[fBeamLine[ibeam]]).c_str(),
		  endtime.c_str());  
	} else {
	  mf::LogInfo("")<<"Pad time for multiwire bundle";
	  endtime=to_iso_extended_string(t1+minutes(1));
	  sprintf(sbuf, "%s/data?b=%s&t0=%s&t1=%s&f=csv", fIFDBURL.c_str(), fBundle[fBeamLine[ibeam]][i].c_str(),
		  to_iso_extended_string(fLastQueryTime[fBeamLine[ibeam]]-minutes(1)).c_str(),
		  endtime.c_str());  
	}
	mf::LogDebug("") <<"Query server:\n"<<sbuf;
	fIFDB.GetData(sbuf,response);	
      }
      ProcessResponse(response, data_map, fBeamLine[ibeam], t1);   
      if (data_map.size() > 0 ) 
	last_proc_time = ToPtime(data_map.rbegin()->first) + milliseconds(fTimeWindowMap[fBeamLine[ibeam]]);
      delete response;
    
      WriteData(fBeamLine[ibeam],data_map);
      data_map.clear();
      if (last_proc_time>fLastQueryTime[fBeamLine[ibeam]]) 
	fLastQueryTime[fBeamLine[ibeam]]=last_proc_time;
      else if ( t1 < microsec_clock::local_time()-minutes(fIFDBLatency) ) 
	fLastQueryTime[fBeamLine[ibeam]]=t1;
      
      if ( t1==tend ) isdone=true;
      t1 = (t1+minutes(10)<tend) ? (t1+minutes(10)) : tend;
    }
  }
}

void beamRun::EndRun(boost::posix_time::ptime tstop) 
{
  fRunHeader.fRunEnd=tstop;
  
  // in case run ends up in the future?
  while (microsec_clock::local_time() < tstop ) {
    Update(microsec_clock::local_time());
    sleep(60);
  }

  //now update until finished
  while ( 1 ) {
    //update beyond tstop to see data (data after 
    //fRunHeader.fRunEnd will not be written to file anyway)
    Update(tstop+minutes(fIFDBLatency)); 
    bool all_done=true;
    for (unsigned int i=0;i<fBeamLine.size();i++) {
      if (ToPtime(fLastBeamHeader[fBeamLine[i]]) < tstop ) all_done=false;
    }
    //finished if all beamlines already have 
    if (all_done) fRunHeader.fFinished=true;
    //or if we are beyond latency
    else if (tstop < microsec_clock::local_time()-minutes(fIFDBLatency)) fRunHeader.fFinished=true;
    
    if (fRunHeader.fFinished) break;
    sleep(60);
  }
  mf::LogInfo("")<<"Closing beam files.";
  for (unsigned int i=0;i<fBeamLine.size();i++) {
    fOut[i].close();
  }

  ofstream infofile;
  stringstream ss;
  ss<<fOutputDirInfo<<"/beam_"
      <<setfill('0')<<setw(7)<<fRunHeader.fRun<<"_"
      <<setfill('0')<<setw(4)<<fRunHeader.fSubRun
      <<".info";
  mf::LogInfo("") <<"Writing beam run summary to "<<ss.str();
  infofile.open(ss.str().c_str());
  infofile<<fRunHeader;
  infofile.close();
}

void beamRun::ProcessResponse(httpResponse* response, beamdatamap_t &data_map, std::string beamline, boost::posix_time::ptime tend)
{
  /****************************************************************
   * Function fills the data_map with data contained in httpResponse.  
   * Returns the timestamp for last processed event. 
   ***************************************************************/

  std::vector<std::string> all_rows;
  boost::split(all_rows, response->memory, boost::is_any_of("\n")); //split into lines
  all_rows.pop_back(); //remove last blank line

  std::map<uint64_t, std::vector<ub_BeamData> > bd_map;

  //need to unpack multiwire data, so will append all_rows which will change the size within the loop =P
  for (unsigned int i=1;i<all_rows.size();i++) {
    //parse single row
    
    //skip rows that are not data (like timestamp,name,units,value(s))
    if (all_rows[i].find("timestamp,")!=std::string::npos ||
	all_rows[i].find("units,")!=std::string::npos ||
	all_rows[i].find("value(s)")!=std::string::npos ) continue;

    std::vector<std::string> row(0);
    boost::split(row, all_rows[i], boost::is_any_of(","));
    uint64_t timestamp;
    if (all_rows[i].find("RAW")!=std::string::npos) {
      MWRData mwrdata;
      std::vector<std::string> unpackedData=mwrdata.unpackMWR(all_rows[i],fMWRTimeOffset);
      all_rows.insert(all_rows.end(),unpackedData.begin(),unpackedData.end());
      continue;
    }
    std::stringstream ss(row[0]);
    ss>>timestamp;
    
    uint32_t secs=uint32_t(timestamp/1000);
    uint16_t msecs=uint16_t(timestamp-(timestamp/1000)*1000);
    if (from_time_t(secs)+hours(fZoneOffset.hours())+microseconds(msecs)<fLastQueryTime[beamline]) continue;
    if (from_time_t(secs)+hours(fZoneOffset.hours())+microseconds(msecs)>tend) continue;

    ub_BeamData dt;
    dt.setDeviceName(row[1]);
    dt.setUnits(row[2]);
    dt.setSeconds(uint32_t(timestamp/1000));
    dt.setMilliSeconds(uint16_t(timestamp-(timestamp/1000)*1000));
    for (unsigned int j=3;j<row.size();j++) dt.pushData(atof(row[j].c_str()));
    
    if (bd_map.find(timestamp) == bd_map.end() ) {
      std::vector<ub_BeamData> tmp;
      tmp.push_back(dt);
      bd_map[timestamp]=tmp;
    } else {
      bd_map[timestamp].push_back(dt);
    }
  }

  uint32_t number_of_bytes_in_record=0;
  uint16_t number_of_devices=0;

  uint32_t last_secs=0;
  uint16_t last_msecs=0;

  std::vector<ub_BeamData> vdt;

  bool is_first=true;
  
  for (auto itr=bd_map.begin(); itr!=bd_map.end();itr++) {
    uint64_t timestamp=itr->first;
    uint32_t secs=uint32_t(timestamp/1000);
    uint16_t msecs=uint16_t(timestamp-(timestamp/1000)*1000);

    //collect all data within fTimeWindowMap[beamline]
    //once we step over fTimeWindowMap, create header and insert into data_map
    //clear vectors and start again
    if ( float(secs-last_secs)+float(msecs-last_msecs)/1000.0 > float(fTimeWindowMap[beamline])/1000. && !is_first) {    
      ub_BeamHeader hdr;
      hdr.setRecordType(8);
      hdr.setEventSignal(std::to_string(fEventTypeMap[beamline]));
      hdr.setSeconds(last_secs);
      hdr.setMilliSeconds(last_msecs);
      hdr.setNumberOfBytesInRecord(number_of_bytes_in_record+sizeof(hdr));
      hdr.setNumberOfDevices(number_of_devices);
      //      std::vector<ub_BeamData>::iterator it=std::unique(vdt.begin(), vdt.end(), beamRun::compareBeamData);

      pair<ub_BeamHeader, std::vector<ub_BeamData> > p(hdr,vdt);
      data_map.insert(p);

      vdt.clear();
      number_of_bytes_in_record=0;
      number_of_devices=0;
    } 

    if (vdt.size() == 0) {
      last_secs=secs;
      last_msecs=msecs;
    }

    vdt.insert(vdt.end(),itr->second.begin(), itr->second.end());
    //    number_of_bytes_in_record += ??;
    number_of_devices+=itr->second.size();

    is_first=false;
  }

}

void beamRun::WriteData(std::string beamline, map<ub_BeamHeader, std::vector<ub_BeamData> > &data_map) 
{
   /****************************************************************
   * Function writes data to file (file is identified by beamline))
   * Only new data written to file even if data_map contains already 
   * stored data.
   ***************************************************************/

  if (data_map.size()==0) return;

  stringstream ss;
  ss<<fOutputDirData<<"/beam_"<<beamline<<"_"
    <<setfill('0')<<setw(7)<<fRunHeader.fRun<<"_"
    <<setfill('0')<<setw(4)<<fRunHeader.fSubRun
    <<".dat.lock";
  
  fopen(ss.str().c_str(),"w");

  int n=0;
  std::map<ub_BeamHeader, std::vector<ub_BeamData> >::iterator it = data_map.begin();

  int ibeamline=0;
  while (beamline!=fBeamLine[ibeamline]) ibeamline++;

  while (it != data_map.end()) {
    if ( !(it->first <= fLastBeamHeader[beamline]) &&
	 (ToPtime(it->first) <= fRunHeader.fRunEnd)) {
      fRunHeader.fCounter[beamline]=fRunHeader.fCounter[beamline]+1;
      //    (*fOA[beamline]) << it->first;
      boost::archive::binary_oarchive oa(fOut[ibeamline]);
      oa<<it->first;
      for (unsigned int i=0;i<it->second.size();i++)  oa<<it->second[i];
	// (*fOA[beamline])<<it->second[i];
      n++;
    } 
    it++;
  }
  
  remove(ss.str().c_str());

  it=data_map.end();
  it--;

  fLastBeamHeader[beamline] = it->first < fLastBeamHeader[beamline] ? fLastBeamHeader[beamline] : it->first;
}

ptime beamRun::ToPtime(ub_BeamHeader bh) 
{ 
  return from_time_t(bh.getSeconds()) + milliseconds(bh.getMilliSeconds()) + fZoneOffset;
}
