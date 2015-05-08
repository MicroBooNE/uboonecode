/**
 * \file UBooneIOVTimeStamp.h
 * 
 * \brief Class def header for a class UBooneIOVTimeStamp
 *
 * @author eberly@fnal.gov
 */

#ifndef UBOONEIOVTIMESTAMP_H
#define UBOONEIOVTIMESTAMP_H

#include "art/Framework/Principal/Event.h"
#include "CalibrationDBI/IOVData/IOVTimeStamp.h"
#include "CalibrationDBI/IOVData/IOVDataConstants.h"
#include "CalibrationDBI/IOVData/IOVDataError.h"
#include <string>

namespace lariov {
  class UBooneIOVTimeStamp : public IOVTimeStamp {
    public: 
      //new constructor
      UBooneIOVTimeStamp(const art::Event& evt);
  };

  //implementation
  UBooneIOVTimeStamp::UBooneIOVTimeStamp(const art::Event& evt) :
    IOVTimeStamp(0,0) {

    std::string time = std::to_string(evt.time().value());

    //microboone stores timestamp as ns from epoch, so there should be 19 digits.
    if (time.length() == 19) {
      //make timestamp conform to database precision
      time = time.substr(0, 10+kMAX_SUBSTAMP_LENGTH);

      //insert decimal point
      time.insert(10,".");

      //finish construction
      IOVTimeStamp tmp = GetFromString(time);
      fStamp    = tmp.Stamp();
      fSubStamp = tmp.SubStamp();
      fDBStamp  = tmp.DBStamp();
    }
    else {
      std::string msg = "UBooneIOVTimeStamp: I do not know how to convert this timestamp: " + time;
      throw IOVDataError(msg);
    }

  }
}

#endif
