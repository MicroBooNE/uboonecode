
// Sanford-Wang xsec fit used for pi+, pi-, K0_s

#include "../WeightCalcCreator.h"
#include "../WeightCalc.h"
#include "uboone/EventWeight/IFDHFileTransfer.h"


#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"

#include <vector>
#include "TH1.h"
#include "TMatrixD.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include <TROOT.h>
#include <TChain.h>

using namespace std;

namespace evwgh {
  class PrimaryHadronSWWeightCalc : public WeightCalc
  {
    evwgh::IFDHFileTransfer IFDH;
  public:
    PrimaryHadronSWWeightCalc();
    void Configure(fhicl::ParameterSet const& p);
    virtual std::vector<std::vector<double> > GetWeight(art::Event & e);
    std::vector<double> xSecFitParameters;
    std::vector< std::vector<double> > fitCovarianceMatrix;
    void ExternalData(std::vector<double>& xSecFitParameters, std::vector< std::vector<double> >& fitCovarianceMatrix);

  private:
    double xsec_sw(TLorentzVector had, TLorentzVector incidentP, std::vector<double> xSecFitParameterstmp, std::vector< std::vector<double> > fitCovarianceMatrixtmp); 
    std::vector<double> ConvertToVector(TArrayD const* array);
    std::string fGenieModuleLabel;
    std::vector<std::string> fParameter_list;
    float fParameter_sigma;
    int primaryHad;
    int fNmultisims;
    TMatrixD crossSection;
    TMatrixD covarianceMatrix;
    CLHEP::RandGaussQ *GaussRandom;
    std::vector<std::vector<double> > fakeData;
    std::string ExternalDataInput;

    DECLARE_WEIGHTCALC(PrimaryHadronSWWeightCalc)
  };
  PrimaryHadronSWWeightCalc::PrimaryHadronSWWeightCalc()
  {
  }

  std::vector<double> PrimaryHadronSWWeightCalc::ConvertToVector(TArrayD const* array) {
    std::vector<double> v(array->GetSize());
    std::copy(array->GetArray(), array->GetArray() + array->GetSize(),
	      v.begin());
    return v;
  } // ConvertToVector()

  void PrimaryHadronSWWeightCalc::ExternalData(std::vector<double>& xSecFitParameters, std::vector< std::vector<double> >& fitCovarianceMatrix)
  {

    TFile file(ExternalDataInput.c_str(),"READ");

    std::string pname;
    pname = "SW/PiPlus/SWPiPlusFitVal";
    TArrayD* SWPiPlusFitValArray = (TArrayD*) file.Get(pname.c_str());
    if (!SWPiPlusFitValArray) {
      throw art::Exception(art::errors::NotFound)
        << "Could not find parameter '" << pname << "' in file '"
        << file.GetPath() << "'\n";
    }

    std::vector<double> SWPiPlusFitVal = PrimaryHadronSWWeightCalc::ConvertToVector(SWPiPlusFitValArray);
    delete SWPiPlusFitValArray;
    pname = "SW/PiPlus/SWPiPlusFitCov";
    TMatrixD* SWPiPlusFitCov = (TMatrixD*) file.Get(pname.c_str());
    if (!SWPiPlusFitCov) {
      throw art::Exception(art::errors::NotFound)
        << "Could not find parameter '" << pname << "' in file '"
        << file.GetPath() << "'\n";
    }

    pname = "SW/PiMinus/SWPiMinusFitVal";
    TArrayD* SWPiMinusFitValArray = (TArrayD*) file.Get(pname.c_str());
    if (!SWPiMinusFitValArray) {
      throw art::Exception(art::errors::NotFound)
        << "Could not find parameter '" << pname << "' in file '"
        << file.GetPath() << "'\n";
    }
    std::vector<double> SWPiMinusFitVal = PrimaryHadronSWWeightCalc::ConvertToVector(SWPiMinusFitValArray);
    delete SWPiMinusFitValArray;

    pname = "SW/PiMinus/SWPiMinusFitCov";
    TMatrixD* SWPiMinusFitCov = (TMatrixD*) file.Get(pname.c_str());

    pname = "SW/K0s/SWK0sFitVal";
    TArrayD* SWK0sFitValArray = (TArrayD*) file.Get(pname.c_str());
    if (!SWK0sFitValArray) {
      throw art::Exception(art::errors::NotFound)
        << "Could not find parameter '" << pname << "' in file '"
        << file.GetPath() << "'\n";
    }
    std::vector<double> SWK0sFitVal = PrimaryHadronSWWeightCalc::ConvertToVector(SWK0sFitValArray);
    delete SWK0sFitValArray;

    pname = "SW/K0s/SWK0sFitCov";
    TMatrixD* SWK0sFitCov = (TMatrixD*) file.Get(pname.c_str());

    if(primaryHad == 211){
	  std::vector< std::vector<double> > v2d;
	  xSecFitParameters = SWPiPlusFitVal;
	  for (int entry=0; entry<9; entry++){
	    std::vector<double> tmp(9);
	    for(Int_t ixsCov =0; ixsCov <9; ++ixsCov) tmp.at(ixsCov) = (*SWPiPlusFitCov)(entry,ixsCov);
	    v2d.push_back(tmp);

	  }

	  fitCovarianceMatrix = v2d;
    }


    if(primaryHad == -211){
	  std::vector< std::vector<double> > v2d;
	  xSecFitParameters = SWPiMinusFitVal;
          for (int entry=0; entry<9 ; ++entry){
	    std::vector<double> tmp(9);
	    for(Int_t ixsCov =0; ixsCov <9; ++ixsCov) tmp.at(ixsCov) = (*SWPiMinusFitCov)(entry,ixsCov);
	    v2d.push_back(tmp);
	  }

	  fitCovarianceMatrix = v2d;

    }

    if(primaryHad == 310){
      std::vector< std::vector<double> > v2d;
      xSecFitParameters = SWK0sFitVal;
      for (int entry=0; entry<9 ; ++entry){
	std::vector<double> tmp(9);
	for(Int_t ixsCov =0; ixsCov <9; ++ixsCov) tmp.at(ixsCov) = (*SWK0sFitCov)(entry,ixsCov);
	v2d.push_back(tmp);
      }

      fitCovarianceMatrix = v2d;
    }

  }



  //Sanford Wang cross section fit from MiniBooNE's Flux paper: Phys. Rev. D 79, 072002
 
  double PrimaryHadronSWWeightCalc::xsec_sw(TLorentzVector had, TLorentzVector incidentP, std::vector<double> xSecFitParameterstmp,std::vector< std::vector<double> > fitCovarianceMatrixtmp)
  {

    
    if(xSecFitParameterstmp.size() != 9)//check number of xSec fit parameters
      {
	throw art::Exception(art::errors::StdException)
	  << "Failed to run, PrimaryHadronSW requires 9 fit parameters";
	return 0;
      }
    
    double a = xSecFitParameterstmp.at(0);
    double b = xSecFitParameterstmp.at(1);
    double c = xSecFitParameterstmp.at(2);
    double d = xSecFitParameterstmp.at(3);
    double e = xSecFitParameterstmp.at(4);
    double f = xSecFitParameterstmp.at(5);
    double g = xSecFitParameterstmp.at(6);
    double h = xSecFitParameterstmp.at(7);
    double i = xSecFitParameterstmp.at(8);

    //    std::cout<<" primary hadron momentum   "<<had.P()<< "incident proton momentum  "<< incidentP.P()<<std::endl;
    //    std::cout<<" primary hadron theta wrt proton  "<<had.Theta()<<std::endl;

    
    double xsec = a*pow(had.P(),b)*(1. - ( had.P()/(incidentP.P()-i) ))
          *exp(-1.*(c*pow(had.P(),d))/pow(incidentP.P(),e)
    		 - f*had.Theta()*(had.P() - g*incidentP.P()*pow(cos(had.Theta()),h) ));
    
    if(xsec < 0) return 0.;

    return xsec;
  }//xsec

  void PrimaryHadronSWWeightCalc::Configure(fhicl::ParameterSet const& p)
  {
    fGenieModuleLabel= p.get< std::string > ("genie_module_label");
    //get configuration for this function
    fhicl::ParameterSet const &pset=p.get<fhicl::ParameterSet> (GetName());
    std::cout << pset.to_string() << std::endl;

    fParameter_list			=   pset.get<std::vector<std::string> >("parameter_list");
    fParameter_sigma		=   pset.get<float>("parameter_sigma");
    primaryHad					=   pset.get<int>("PrimaryHadronGeantCode");
    fNmultisims      = 			pset.get<int>("number_of_multisims");
    std::string dataInput       =   pset.get< std::string >("ExternalData");
    ExternalDataInput = IFDH.fetch(dataInput);
   
    //Prepare random generator
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::RandGaussQ GaussRandom(rng->getEngine(GetName()));
 
    PrimaryHadronSWWeightCalc::ExternalData( xSecFitParameters, fitCovarianceMatrix);

   fakeData = WeightCalc::MultiGaussianSmearing(xSecFitParameters, fitCovarianceMatrix, fNmultisims, GaussRandom);

  }

   std::vector<std::vector<double> >  PrimaryHadronSWWeightCalc::GetWeight(art::Event & e)

  {

    //calculate weight(s) here 
    std::vector<std::vector<double> > weight;

    //get the MC generator information out of the event       
    //these are all handles to mc information.
    art::Handle< std::vector<simb::MCFlux> > mcFluxHandle;

    //actually go and get the stuff
    e.getByLabel(fGenieModuleLabel,mcFluxHandle);
   
    std::vector<simb::MCFlux> const& fluxlist = *mcFluxHandle;
    std::vector<double> weightsForOneNeutrino;

    for ( unsigned int inu = 0; inu < fluxlist.size(); ++inu )//each neutrino
      {

	if (fluxlist[inu].ftptype != primaryHad){ ///continue;
	  for(int it =0; it<fNmultisims; it++) weightsForOneNeutrino.push_back(1.);
	  weight.push_back(weightsForOneNeutrino);
	  //continue;
	}

	if (fluxlist[inu].ftptype != primaryHad) continue;

    	//determine mass
    	double m = 0; //mass
    	switch (abs(primaryHad))
	  {
	  case 211:
	    m = 0.13957010; //pi+ and pi- mass [GeV/c^2]
	    break;
	  case 310:
	    m = 0.4976140;//k0_s mass [GeV/c^2]
	    break;
	  case 321:
	    m = 0.4936770; //kaon+ and kaon- mass [GeV/c^2]
	    break;
	  default:
	    //	    throw art::Exception(art::errors::StdException)
	    {
	      std::cout<< "Failed, please add the mass of the requested Primary Hadron."<< primaryHad<<"   ******************"<<std::endl;
	      continue;
	    }

	  }
	/*
	double px = fluxlist[inu].fpdpx;
	double py = fluxlist[inu].fpdpy;
	double pz = fluxlist[inu].fpdpz;
	std::cout<<"particle type  "<<fluxlist[inu].fptype<<std::endl;
	std::cout<<"particle type accoding default "<<fluxlist[inu].ftptype<<std::endl;
	*/

	
	double px = fluxlist[inu].ftpx;
	double py = fluxlist[inu].ftpy;
	double pz = fluxlist[inu].ftpz;
	double e  = sqrt(px*px + py*py + pz*pz + m*m);
	TLorentzVector hadronVec;
	hadronVec.SetPxPyPzE(px,py,pz,e);

	//find hadron's incident proton info
	TLorentzVector incidentPVec;
	//	double px_proton = fluxlist[inu].fbeampx;
	//	double py_proton = fluxlist[inu].fbeampy;
	//	double pz_proton = fluxlist[inu].fbeampz;
	double px_proton = 0.;
	double py_proton = 0.;
	double pz_proton = 8.89;/// nominal peak of proton momentum in GeV, assuming along z-direction
	double m_proton  = 0.9382720;
	double e_proton  = sqrt(px_proton*px_proton + py_proton*py_proton + pz_proton*pz_proton + m_proton*m_proton);  
	incidentPVec.SetPxPyPzE(px_proton,py_proton,pz_proton,e_proton);

	PrimaryHadronSWWeightCalc::ExternalData( xSecFitParameters, fitCovarianceMatrix);

	double crossSection = xsec_sw(hadronVec, incidentPVec, xSecFitParameters,fitCovarianceMatrix);

	for(unsigned int fitparam = 0;fitparam < xSecFitParameters.size(); fitparam++){
	  for(vector<vector<double> >::iterator itf = fakeData.begin(); itf != fakeData.end(); ++itf)
	  {
	    std::vector<double> tmp = *itf;
	    //	    std::cout<<" fake fit parameter: "<<tmp[fitparam]<<" at parameter : "<< fitparam<<std::endl;

	  }
	  
	}

	//	if(primaryHad ==211) h4->Fill(crossSection);
      
	//	std::cout<<"crossSection SW  "<<crossSection<<"  hadron type "<< fluxlist[inu].ftptype<<std::endl;


	int kint = 0;
	for(vector<vector<double> >::iterator it = fakeData.begin(); it != fakeData.end(); ++it)
	  {

	    kint++;
	    std::vector<double> tmp = *it;
	    double fakeCrossSection = xsec_sw(hadronVec, incidentPVec, tmp,fitCovarianceMatrix);

	    //  for(Int_t ixsCov =0; ixsCov <9; ++ixsCov) std::cout<<"fakedata fitparam(entry,ixsCov) "<<"  "<<tmp[ixsCov]<<std::endl;

	    if(fakeCrossSection == 0) std::cout<<" fakeCrossSection SW =0, for   "<<primaryHad<<std::endl;
	    if(fakeCrossSection/crossSection > 2 ) std::cout<<"weight bigger than 2 for "<<primaryHad<<" weight " <<fakeCrossSection/crossSection<<std::endl;

	    if(primaryHad != 211) std::cout<<"no pi+ weight SW  "<< fakeCrossSection/crossSection << "   primary hadron  "<<  primaryHad <<"  mulstissim   "<< kint <<std::endl;
	    //if(primaryHad == 211) std::cout<<"pi+ fakeCrossSection SW  "<<fakeCrossSection<<std::endl;
	    //if(primaryHad == 211) std::cout<<"pi+ weight SW  "<< fakeCrossSection/crossSection << "   primary hadron  "<<  primaryHad <<"  mulstissim   "<< kint <<std::endl;


	    weightsForOneNeutrino.push_back(fakeCrossSection/crossSection);
	  }///over number of sims

	kint =0;
       	weight.push_back(weightsForOneNeutrino);


      } //loop over each neutrino

	///////////
    /*
   for ( unsigned int inu = 0; inu < fluxlist.size(); ++inu )//each neutrino
      {

	std::vector<double> weightsForOneNeutrino;
	if (fluxlist[inu].ftptype != primaryHad){
	  for(int it =0; it<fNmultisims; it++)
	    {
	      weightsForOneNeutrino.push_back(1.);
	    }
	  weight.push_back(weightsForOneNeutrino);
	}
      } //loop over each neutrino
    */

    return weight;
  }

  REGISTER_WEIGHTCALC(PrimaryHadronSWWeightCalc)
}
