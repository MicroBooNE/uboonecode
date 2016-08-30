#ifndef UBOONEELECTRONLIFETIMEPROVIDER_CXX
#define UBOONEELECTRONLIFETIMEPROVIDER_CXX

#include "UbooneElectronLifetimeProvider.h"
#include "larevt/CalibrationDBI/Providers/WebError.h"
#include "larevt/CalibrationDBI/IOVData/IOVDataConstants.h"

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
      std::cout << "Using default pedestal values\n";
      float default_exppar     = p.get<float>("DefaultExpPar");
      float default_constpar   = p.get<float>("DefaultConstPar");
      
      ElectronLifetime default_pars(0);
      
      default_pars.SetExpPar(default_exppar);
      default_pars.SetConstPar(default_constpar);
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
      ElectronLifetime dp(0);
      while (std::getline(file, line)) {
        size_t current_comma = line.find(',');
        DBChannelID_t ch = 0;	
	float exp_par = std::stof(line.substr(current_comma+1, line.find(',',current_comma+1)));
	
	current_comma = line.find(',',current_comma+1);
	float const_par = std::stof(line.substr(current_comma+1, line.find(',',current_comma+1)));

	dp.SetChannel(ch);
	dp.SetExpPar(exp_par);
        dp.SetConstPar(const_par);
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

    double exp_par, const_par;
    fFolder->GetNamedChannelData(0, "exponential",     exp_par);
    fFolder->GetNamedChannelData(0, "constant", const_par);

    ElectronLifetime pd(0);
    pd.SetExpPar( (float)exp_par );
    pd.SetConstPar( (float)const_par );

    fData.AddOrReplaceRow(pd);

    return true;

  }
  
  const ElectronLifetime& UbooneElectronLifetimeProvider::Lifetime() const {     
    return fData.GetRow(0);
  }
      
  float UbooneElectronLifetimeProvider::ExpPar() const {
    return this->Lifetime().ExpPar();
  }
  
  float UbooneElectronLifetimeProvider::ConstPar() const {
    return this->Lifetime().ConstPar();
  }



}//end namespace lariov
	
#endif
        
