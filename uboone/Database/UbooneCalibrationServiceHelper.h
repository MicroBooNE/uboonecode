#ifndef UBOONECALIBRATIONSERVICEHELPER_H
#define UBOONECALIBRATIONSERVICEHELPER_H

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "uboone/DataOverlay/DataOverlayProducts/EventMixingSummary.h"

namespace lariov {

  class UbooneCalibrationServiceHelper {
  
    public:
    
      UbooneCalibrationServiceHelper(fhicl::ParameterSet const& pset);
      ~UbooneCalibrationServiceHelper() {}
      
      DBTimeStamp_t GetTimeStamp(const art::Event& evt, std::string dbname="") const;
      
    private:
      
      std::string fMixingModuleLabel;
        
  };
   
  UbooneCalibrationServiceHelper::UbooneCalibrationServiceHelper(fhicl::ParameterSet const& pset)
  : fMixingModuleLabel(pset.get<std::string>("EventMixingModuleLabel"))
  {}


  DBTimeStamp_t UbooneCalibrationServiceHelper::GetTimeStamp(const art::Event& evt, std::string dbname /*=""*/ ) const {

    art::Handle< std::vector<mix::EventMixingSummary> > eventMixingSummary;
    evt.getByLabel(fMixingModuleLabel, eventMixingSummary);
    if (eventMixingSummary.isValid() && eventMixingSummary->size()>0) {
      if (eventMixingSummary->size() > 1) {
	std::cout<<"  INFO: "<<eventMixingSummary->size()<<" EventMixingSummary objects"<<std::endl;
      }
      art::Timestamp time_stamp = eventMixingSummary->front().Timestamp();
      DBTimeStamp_t  the_time = time_stamp.timeHigh()*1.0e9 + time_stamp.timeLow();

      std::cout<<"Using EventMixingSummary timestamp to query "<<dbname<<" database: "<<the_time<<std::endl;
      return the_time;
    }
    else {    
      art::Timestamp time_stamp = evt.time();
      DBTimeStamp_t the_time;
      if (evt.isRealData()) the_time = time_stamp.timeHigh()*1.0e9 + time_stamp.timeLow();
      else                  the_time = time_stamp.value();

      std::cout<<"Using art::Event timestamp to query "<<dbname<<" database: "<<the_time<<std::endl;
      return the_time;
    }
  }
  
} //end namespace lariov

#endif
