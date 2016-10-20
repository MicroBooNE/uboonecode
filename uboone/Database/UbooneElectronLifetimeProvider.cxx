#ifndef UBOONEELECTRONLIFETIMEPROVIDER_CXX
#define UBOONEELECTRONLIFETIMEPROVIDER_CXX

#include "UbooneElectronLifetimeProvider.h"
//#include "larevt/CalibrationDBI/Providers/WebError.h"

// art/LArSoft libraries
#include "cetlib/exception.h"

//C/C++
#include <fstream>

namespace lariov {

  //constructors
  UbooneElectronLifetimeProvider::UbooneElectronLifetimeProvider(const std::string& foldername, 
      			      			   const std::string& url, 
			      			   const std::string& tag /*=""*/) : 
    DatabaseRetrievalAlg(foldername, url, tag),
    fDataSource(DataSource::Database) {
    
    fData.Clear();
    IOVTimeStamp tmp = IOVTimeStamp::MaxTimeStamp();
    tmp.SetStamp(tmp.Stamp()-1, tmp.SubStamp());
    fData.SetIoV(tmp, IOVTimeStamp::MaxTimeStamp());
  }
	
      
  UbooneElectronLifetimeProvider::UbooneElectronLifetimeProvider(fhicl::ParameterSet const& p) :
    DatabaseRetrievalAlg(p.get<fhicl::ParameterSet>("DatabaseRetrievalAlg")) {	
    
    this->Reconfigure(p);
  }
      
  void UbooneElectronLifetimeProvider::Reconfigure(fhicl::ParameterSet const& p) {
    
    this->DatabaseRetrievalAlg::Reconfigure(p.get<fhicl::ParameterSet>("DatabaseRetrievalAlg"));
    fData.Clear();
    IOVTimeStamp tmp = IOVTimeStamp::MaxTimeStamp();
    tmp.SetStamp(tmp.Stamp()-1, tmp.SubStamp());
    fData.SetIoV(tmp, IOVTimeStamp::MaxTimeStamp());

    bool UseDB      = p.get<bool>("UseDB", false);
    bool UseFile   = p.get<bool>("UseFile", false);
    std::string fileName = p.get<std::string>("FileName", "");

    //priority:  (1) use db, (2) use table, (3) use defaults
    //If none are specified, use defaults
    if ( UseDB )      fDataSource = DataSource::Database;
    else if (UseFile) fDataSource = DataSource::File;
    else              fDataSource = DataSource::Default;

    if (fDataSource == DataSource::Default) {
      std::cout << "Using default lifetime fit values\n";
      float default_expoffset      = p.get<float>("DefaultExpOffset");
      float default_timeconstant   = p.get<float>("DefaultTimeConstant");
      float default_expoffseterr      = p.get<float>("DefaultExpOffsetErr");
      float default_timeconstanterr   = p.get<float>("DefaultTimeConstantErr");
      
      ElectronLifetimeContainer default_pars(fLifetimeChannel);
      
      default_pars.SetExpOffset(default_expoffset);
      default_pars.SetTimeConstant(default_timeconstant);
      default_pars.SetExpOffsetErr(default_expoffseterr);
      default_pars.SetTimeConstantErr(default_timeconstanterr);
      fData.AddOrReplaceRow(default_pars);
 
    }
    else if (fDataSource == DataSource::File) {
      std::cout << "Using pedestals from local file: "<<fileName<<"\n";
      std::ifstream file(fileName);
      if (!file) {
        throw cet::exception("UbooneElectronLifetimeProvider")
	  << "File "<<fileName<<" is not found.";
      }
      
      std::string line;
      ElectronLifetimeContainer dp(fLifetimeChannel);
      while (std::getline(file, line)) {
        size_t current_comma = line.find(',');	
	float exp_offset = std::stof(line.substr(current_comma+1, line.find(',',current_comma+1)));
	
	current_comma = line.find(',',current_comma+1);
	float exp_offset_err = std::stof(line.substr(current_comma+1, line.find(',',current_comma+1)));
	
	current_comma = line.find(',',current_comma+1);
	float time_constant = std::stof(line.substr(current_comma+1, line.find(',',current_comma+1)));
	
	current_comma = line.find(',',current_comma+1);
	float time_constant_err = std::stof(line.substr(current_comma+1, line.find(',',current_comma+1)));

	dp.SetChannel(fLifetimeChannel);
	dp.SetExpOffset(exp_offset);
        dp.SetTimeConstant(time_constant);
	dp.SetExpOffsetErr(exp_offset_err);
        dp.SetTimeConstantErr(time_constant_err);
	fData.AddOrReplaceRow(dp);
	
	break; //only should have one line
      }
    } // if source from file
    else { 
      std::cout << "Using electron lifetime from conditions database\n";
    }
  }


  bool UbooneElectronLifetimeProvider::Update(DBTimeStamp_t ts) {
    
    if (fDataSource != DataSource::Database) return false;
      
    if (!this->UpdateFolder(ts)) return false;

    //DBFolder was updated, so now update the Snapshot
    fData.Clear();
    fData.SetIoV(this->Begin(), this->End());

    double exp_offset, exp_offset_err, time_constant, time_constant_err;
    fFolder->GetNamedChannelData(fLifetimeChannel, "exponential_offset",     exp_offset);
    fFolder->GetNamedChannelData(fLifetimeChannel, "time_constant",          time_constant);
    fFolder->GetNamedChannelData(fLifetimeChannel, "err_exponential_offset", exp_offset_err);
    fFolder->GetNamedChannelData(fLifetimeChannel, "err_time_constant",      time_constant_err);

    ElectronLifetimeContainer pd(fLifetimeChannel);
    pd.SetExpOffset( (float)exp_offset );
    pd.SetTimeConstant( (float)time_constant );
    pd.SetExpOffsetErr( (float)exp_offset_err );
    pd.SetTimeConstantErr( (float)time_constant_err );

    fData.AddOrReplaceRow(pd);

    return true;

  }
  
  const ElectronLifetimeContainer& UbooneElectronLifetimeProvider::LifetimeContainer() const {     
    return fData.GetRow(fLifetimeChannel);
  }
  
  float UbooneElectronLifetimeProvider::Lifetime(float t) const {
    float offset = this->ExpOffset();
    float c = this->TimeConstant();
    return (1.0+exp(offset))/(1.0+exp(offset + c*t));
  }
  
  float UbooneElectronLifetimeProvider::LifetimeErr(float t) const {
    float offset = this->ExpOffset();
    float c = this->TimeConstant();
    float offset_err = this->ExpOffsetErr();
    float c_err = this->TimeConstantErr();
    
    float vb = pow(exp(offset)*(1.0-exp(c*t))*offset_err,2.0);
    float vc = pow(t*(1.0+exp(offset))*exp(offset+c*t)*c_err,2.0);
    return sqrt(vb+vc)/pow(1.0+exp(offset+c*t),2.0);
  }
   
  float UbooneElectronLifetimeProvider::Purity() const {
    float offset = this->ExpOffset();
    float c = this->TimeConstant();
    return (1.0 + exp(offset + 2200.0*c))/(1.0 + exp(offset));
  }
  
  float UbooneElectronLifetimeProvider::PurityErr() const {
    float offset = this->ExpOffset();
    float c = this->TimeConstant();
    float offset_err = this->ExpOffsetErr();
    float c_err = this->TimeConstantErr();
    
    float vb = pow(offset_err*exp(offset)*(exp(2200.0*c)-1.0)/pow(1.0+exp(offset),2.0),2.0);
    float vc = pow(c_err*2200.0*exp(offset+2200.0*c)/(1.0+exp(offset)),2.0);
    return sqrt(vb+vc);
  }   
      
  float UbooneElectronLifetimeProvider::ExpOffset() const {
    return this->LifetimeContainer().ExpOffset();
  }
  
  float UbooneElectronLifetimeProvider::TimeConstant() const {
    return this->LifetimeContainer().TimeConstant();
  }

  float UbooneElectronLifetimeProvider::ExpOffsetErr() const {
    return this->LifetimeContainer().ExpOffsetErr();
  }
  
  float UbooneElectronLifetimeProvider::TimeConstantErr() const {
    return this->LifetimeContainer().TimeConstantErr();
  }

}//end namespace lariov
	
#endif
        
