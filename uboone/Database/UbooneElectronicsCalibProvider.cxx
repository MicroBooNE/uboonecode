#ifndef UBOONEELECTRONICSCALIBPROVIDER_CXX
#define UBOONEELECTRONICSCALIBPROVIDER_CXX

#include "UbooneElectronicsCalibProvider.h"
#include "larevt/CalibrationDBI/Providers/WebError.h"

// art/LArSoft libraries
#include "cetlib/exception.h"
#include "larcore/Geometry/Geometry.h"


#include <fstream>

namespace lariov {

  //constructor      
  UbooneElectronicsCalibProvider::UbooneElectronicsCalibProvider(fhicl::ParameterSet const& p) :
    DatabaseRetrievalAlg(p.get<fhicl::ParameterSet>("DatabaseRetrievalAlg")) {	
    
    this->Reconfigure(p);
  }
      
  void UbooneElectronicsCalibProvider::Reconfigure(fhicl::ParameterSet const& p) {
    
    this->DatabaseRetrievalAlg::Reconfigure(p.get<fhicl::ParameterSet>("DatabaseRetrievalAlg"));
    fData.Clear();
    IOVTimeStamp tmp = IOVTimeStamp::MaxTimeStamp();
    tmp.SetStamp(tmp.Stamp()-1, tmp.SubStamp());
    fData.SetIoV(tmp, IOVTimeStamp::MaxTimeStamp());

    bool UseDB      = p.get<bool>("UseDB", false);
    bool UseFile    = p.get<bool>("UseFile", false);
    std::string fileName = p.get<std::string>("FileName", "");

    //priority:  (1) use db, (2) use table, (3) use defaults
    //If none are specified, use defaults
    if ( UseDB )      fDataSource = DataSource::Database;
    else if (UseFile) fDataSource = DataSource::File;
    else              fDataSource = DataSource::Default;

    if (fDataSource == DataSource::Default) {
      float default_gain     = p.get<float>("DefaultGain");
      float default_gain_err = p.get<float>("DefaultGainErr");
      float default_st       = p.get<float>("DefaultShapingTime");
      float default_st_err   = p.get<float>("DefaultShapingTimeErr");

      ElectronicsCalib defaultCalib(0);
      CalibrationExtraInfo extra_info("ElectronicsCalib");
      extra_info.AddOrReplaceBoolData("is_misconfigured", false);

      defaultCalib.SetGain(default_gain);
      defaultCalib.SetGainErr(default_gain_err);
      defaultCalib.SetShapingTime(default_st);
      defaultCalib.SetShapingTimeErr(default_st_err);     
      defaultCalib.SetExtraInfo(extra_info);
      
      art::ServiceHandle<geo::Geometry> geo;
      geo::wire_id_iterator itW = geo->begin_wire_id();
      for (; itW != geo->end_wire_id(); ++itW) {
	DBChannelID_t ch = geo->PlaneWireToChannel(*itW);
	defaultCalib.SetChannel(ch);
	fData.AddOrReplaceRow(defaultCalib);
      }
      
    }
    else if (fDataSource == DataSource::File) {
      std::cout << "Using electronics calibrations from local file: "<<fileName<<"\n";
      std::ifstream file(fileName);
      if (!file) {
        throw cet::exception("UbooneElectronicsCalibProvider")
          << "File "<<fileName<<" is not found.";
      }
      
      std::string line;
      ElectronicsCalib dp(0);
      while (std::getline(file, line)) {
        size_t current_comma = line.find(',');
        DBChannelID_t ch = (DBChannelID_t)std::stoi(line.substr(0, current_comma));     
        float gain = std::stof(line.substr(current_comma+1, line.find(',',current_comma+1)));
        
        current_comma = line.find(',',current_comma+1);
        float gain_err = std::stof(line.substr(current_comma+1, line.find(',',current_comma+1)));
	
	current_comma = line.find(',',current_comma+1);
        float shaping_time = std::stof(line.substr(current_comma+1, line.find(',',current_comma+1)));
	
	current_comma = line.find(',',current_comma+1);
        float shaping_time_err = std::stof(line.substr(current_comma+1, line.find(',',current_comma+1)));
        
	current_comma = line.find(',',current_comma+1);
        int is_misconfigured = std::stoi(line.substr(current_comma+1, line.find(',',current_comma+1)));
        
	CalibrationExtraInfo extra_info("ElectronicsCalib");
	extra_info.AddOrReplaceBoolData("is_misconfigured", (bool)is_misconfigured);

        dp.SetChannel(ch);
        dp.SetGain(gain);
        dp.SetGainErr(gain_err);
	dp.SetShapingTime(shaping_time);
        dp.SetShapingTimeErr(shaping_time_err);
	dp.SetExtraInfo(extra_info);
        
        fData.AddOrReplaceRow(dp);
      }
    }
    else {
      std::cout << "Using electronics calibrations from conditions database"<<std::endl;
    }
  }

  bool UbooneElectronicsCalibProvider::Update(DBTimeStamp_t ts) {
    
    if (fDataSource != DataSource::Database) return false;
      
    if (!this->UpdateFolder(ts)) return false;

    //DBFolder was updated, so now update the Snapshot
    fData.Clear();
    fData.SetIoV(this->Begin(), this->End());

    std::vector<DBChannelID_t> channels;
    fFolder->GetChannelList(channels);
    for (auto it = channels.begin(); it != channels.end(); ++it) {

      double gain, gain_err, shaping_time, shaping_time_err;
      bool is_misconfigured;
      fFolder->GetNamedChannelData(*it, "gain",     gain);
      fFolder->GetNamedChannelData(*it, "gain_err", gain_err); 
      fFolder->GetNamedChannelData(*it, "shaping_time",     shaping_time);
      fFolder->GetNamedChannelData(*it, "shaping_time_err", shaping_time_err); 
      fFolder->GetNamedChannelData(*it, "is_misconfigured", is_misconfigured);
      
      ElectronicsCalib pg(*it);
      CalibrationExtraInfo extra_info("ElectronicsCalib");
      extra_info.AddOrReplaceBoolData("is_misconfigured", is_misconfigured);
      
      pg.SetGain( (float)gain );
      pg.SetGainErr( (float)gain_err );
      pg.SetShapingTime( (float)shaping_time );
      pg.SetShapingTimeErr( (float)shaping_time_err );
      pg.SetExtraInfo(extra_info);

      fData.AddOrReplaceRow(pg);
    }

    return true;
  }
  
  const ElectronicsCalib& UbooneElectronicsCalibProvider::ElectronicsCalibObject(DBChannelID_t ch) const { 
    return fData.GetRow(ch);
  }
      
  float UbooneElectronicsCalibProvider::Gain(DBChannelID_t ch) const {
    return this->ElectronicsCalibObject(ch).Gain();
  }
  
  float UbooneElectronicsCalibProvider::GainErr(DBChannelID_t ch) const {
    return this->ElectronicsCalibObject(ch).GainErr();
  }
  
  float UbooneElectronicsCalibProvider::ShapingTime(DBChannelID_t ch) const {
    return this->ElectronicsCalibObject(ch).ShapingTime();
  }
  
  float UbooneElectronicsCalibProvider::ShapingTimeErr(DBChannelID_t ch) const {
    return this->ElectronicsCalibObject(ch).ShapingTimeErr();
  }
  
  CalibrationExtraInfo const& UbooneElectronicsCalibProvider::ExtraInfo(DBChannelID_t ch) const {
    return this->ElectronicsCalibObject(ch).ExtraInfo();
  }


}//end namespace lariov
	
#endif
        
