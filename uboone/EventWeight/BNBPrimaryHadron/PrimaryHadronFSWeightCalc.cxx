//Feynmann scaling xsec fit for kaon+

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
  class PrimaryHadronFSWeightCalc : public WeightCalc
  {
    evwgh::IFDHFileTransfer IFDH;
  public:
    PrimaryHadronFSWeightCalc();
    void Configure(fhicl::ParameterSet const& p);
    virtual std::vector<std::vector<double> > GetWeight(art::Event & e);
    std::vector<double> xSecFitParameters;
    std::vector< std::vector<double> > fitCovarianceMatrix;
    void ExternalData(std::vector<double>& xSecFitParameters, std::vector< std::vector<double> >& fitCovarianceMatrix);

  private:

    CLHEP::RandGaussQ *GaussRandom;
    double xsec_fs(TLorentzVector had, TLorentzVector incidentP, std::vector<double> xSecFitParameterstmp,std::vector< std::vector<double> > fitCovarianceMatrixtmp);
    std::vector<double> ConvertToVector(TArrayD const* array);
    std::vector< std::vector<double> > fakeData;				 		
    std::string fGenieModuleLabel;
    std::vector<std::string> fParameter_list;
    float fParameter_sigma;
    int fNmultisims;
    int primaryHad;
    TMatrixD crossSection;
    TMatrixD covarianceMatrix;
    std::string ExternalDataInput;
    
     DECLARE_WEIGHTCALC(PrimaryHadronFSWeightCalc)
  };
  PrimaryHadronFSWeightCalc::PrimaryHadronFSWeightCalc()
  {
  }

  std::vector<double> PrimaryHadronFSWeightCalc::ConvertToVector(TArrayD const* array) {
    std::vector<double> v(array->GetSize());
    std::copy(array->GetArray(), array->GetArray() + array->GetSize(),
	      v.begin());
    return v;
  } // ConvertToVector()

  void PrimaryHadronFSWeightCalc::ExternalData(std::vector<double>& xSecFitParameters, std::vector< std::vector<double> >& fitCovarianceMatrix)
  {

    TFile file(ExternalDataInput.c_str(),"READ");

    std::string pname;

    pname = "FS/KPlus/FSKPlusFitVal";
    TArrayD* FSKPlusFitValArray = (TArrayD*) file.Get(pname.c_str());
    if (!FSKPlusFitValArray) {
      throw art::Exception(art::errors::NotFound)
        << "Could not find parameter '" << pname << "' in file '"
        << file.GetPath() << "'\n";
    }
    std::vector<double> FSKPlusFitVal = PrimaryHadronFSWeightCalc::ConvertToVector(FSKPlusFitValArray);
    delete FSKPlusFitValArray;

    pname = "FS/KPlus/FSKPlusFitCov";
    TMatrixD* FSKPlusFitCov = (TMatrixD*) file.Get(pname.c_str());


    //// Get fit values and covariances for specific particle type and method ///////

    if(primaryHad == 321){
      std::vector< std::vector<double> > v2d;
      xSecFitParameters = FSKPlusFitVal;
      for (int entry=0; entry<7 ; ++entry){
	std::vector<double> tmp(7);
	for(Int_t ixsCov =0; ixsCov <7; ++ixsCov) tmp.at(ixsCov) = (*FSKPlusFitCov)(entry,ixsCov);
	v2d.push_back(tmp);
      }

      fitCovarianceMatrix = v2d;
    }

    else std::cout<<"primary hadron code"<<  primaryHad <<std::endl;

  }


  //cross section fit from MiniBooNE's Flux paper: Phys. Rev. D 79, 072002

  double PrimaryHadronFSWeightCalc::xsec_fs(TLorentzVector had, TLorentzVector incidentP, std::vector<double> xSecFitParameterstmp,std::vector< std::vector<double> > fitCovarianceMatrixtmp)
  {

    //   PrimaryHadronFSWeightCalc::ExternalData( xSecFitParameters, fitCovarianceMatrix);
    if(xSecFitParameterstmp.size() != 7)//check number of xSec fit parameters
      {
	throw art::Exception(art::errors::StdException)
	  << "Failed to run, PrimaryHadronFS requires 7 fit parameters";
	return 0;
      }

    //Compute the Feynman scaling variable x_f = p_cm_parallel_hadron / p_cm_parallel_hadron_MAX

    double E_cm = sqrt(2.*incidentP.M()*incidentP.M() + 2.*incidentP.M()*incidentP.E());///corrected
    double gamma = (incidentP.E()+incidentP.M())/E_cm;
    double gammaBeta = incidentP.P()/E_cm;
		
    //To produce this particular hadron, what is it's maximum
    //CM energy? First figure out the max mass energy, it varies by
    //the hadron produced...

    double maxMassEnergy = 0.;
    if(abs(primaryHad)==311 || abs(primaryHad)==310 || abs(primaryHad)==130)//K0s//only data for k+
      maxMassEnergy = 2.130;//Sigma+ + p
    else if(primaryHad==321)//K+
      maxMassEnergy = 2.053;// lambda + p
    else // assume a pion
      maxMassEnergy = 2.*incidentP.M(); // p + p		

    //maximum energy available for hadron in CM frame
    double E_max_cm = (E_cm*E_cm + had.M()*had.M() - maxMassEnergy*maxMassEnergy)/(2.*E_cm);
    //maximum momentum available for hadron in CM frame
    double p_max_cm = sqrt(E_max_cm*E_max_cm - had.M()*had.M());
			
    double pT = had.P()*sin(had.Theta());
    double pParallel = had.P()*cos(had.Theta());
    double p_cm_parallel = gamma*pParallel - gammaBeta*had.E();
		
    double xF = p_cm_parallel / p_max_cm;
 
    if(fabs(xF) > 1.) return 0.;

  double a = xSecFitParameterstmp.at(0);//c1
  double b = xSecFitParameterstmp.at(1);//c2
  double c = xSecFitParameterstmp.at(2);//c3
  double d = xSecFitParameterstmp.at(3);//c4
  double e = xSecFitParameterstmp.at(4);//c5
  double f = xSecFitParameterstmp.at(5);//c6
  double g = xSecFitParameterstmp.at(6);//c7


  //  std::cout<<" primary hadron momentum   "<<had.P()<< "incident proton momentum  "<< incidentP.P()<<std::endl;
  //  std::cout<<" primary hadron theta wrt proton  "<<had.Theta()<<std::endl;


    double xsec = (had.P()*had.P()/had.E()) * a* exp( -1.*b*pT - c*pow(fabs(xF),d) - e*pT*pT)
      * exp( -1.*g*pow(fabs(pT*xF),f) );

    if(xsec < 0) return 0.;

    return xsec;
  }///xsec

  void PrimaryHadronFSWeightCalc::Configure(fhicl::ParameterSet const& p)
  {
    fGenieModuleLabel= p.get< std::string > ("genie_module_label");
    //get configuration for this function
    fhicl::ParameterSet const &pset=p.get<fhicl::ParameterSet> (GetName());
    std::cout << pset.to_string() << std::endl;

    fParameter_list		=   pset.get<std::vector<std::string> >("parameter_list");
    fParameter_sigma		=   pset.get<float>("parameter_sigma");
    fNmultisims			=   pset.get<int>("number_of_multisims");
    primaryHad			=   pset.get<int>("PrimaryHadronGeantCode");//debiera ser un vector
    std::string dataInput       =   pset.get< std::string >("ExternalData");
    ExternalDataInput = IFDH.fetch(dataInput);
 
    //Prepare random generator
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::RandGaussQ GaussRandom(rng->getEngine(GetName()));
 
    PrimaryHadronFSWeightCalc::ExternalData( xSecFitParameters, fitCovarianceMatrix);
 
    //  for(Int_t entry =0; entry <7; ++entry){
 
     //    std::cout<<"xSecFitParameters "<< xSecFitParameters.at(entry)<< " entry   "<<entry<<std::endl;

     //   std::cout<<" size cov   "<<fitCovarianceMatrix.size()<<std::endl;
     //   std::cout<<" size cov at 0   "<<fitCovarianceMatrix.at(0).size()<<std::endl;

    //  for(Int_t ixsCov =0; ixsCov <7; ++ixsCov) std::cout<<"fitCov(entry,ixsCov) "<< " row   "<< entry<<"  column  "<< ixsCov<<"  "<<fitCovarianceMatrix.at(entry)[ixsCov]<<"   "<< primaryHad<<std::endl;
    //  }

    fakeData = WeightCalc::MultiGaussianSmearing(xSecFitParameters, fitCovarianceMatrix, fNmultisims, GaussRandom);		
  
  }

  std::vector<std::vector<double> > PrimaryHadronFSWeightCalc::GetWeight(art::Event & e)
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
	if (fluxlist[inu].ftptype != primaryHad){/// continue;
	  for(int it =0; it<fNmultisims; it++) weightsForOneNeutrino.push_back(1.);
	  weight.push_back(weightsForOneNeutrino);
	  continue;
	}
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
	    m = 0.4936770; //kaon+ and kaon- mass [GeV/c^2]//check where it is the k+/- data
	    break;
	  default:
	    //	    throw art::Exception(art::errors::StdException)
	    {
	      std::cout<< "Failed, please add the mass of the requested Primary Hadron."<< primaryHad<<"   ******************"<<std::endl;
	      continue;
	    }
	  }

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

	PrimaryHadronFSWeightCalc::ExternalData( xSecFitParameters, fitCovarianceMatrix);

	double crossSection = xsec_fs(hadronVec, incidentPVec, xSecFitParameters,fitCovarianceMatrix);

	//	std::cout<<"crossSection FS  "<<crossSection<<"  hadron type "<< fluxlist[inu].ftptype<<std::endl;

	//	weightsForOneNeutrino;////***poner por defecto a 1!!!

	for(vector<vector<double> >::iterator it = fakeData.begin(); it != fakeData.end(); ++it)
	  {

	    std::vector<double> tmp = *it;
	    double fakeCrossSection = xsec_fs(hadronVec, incidentPVec, tmp,fitCovarianceMatrix);

	    //  std::cout<<"fakeCrossSection FS  "<<fakeCrossSection<<std::endl;

	    weightsForOneNeutrino.push_back(fakeCrossSection/crossSection);
	  }///over number of sims
    	weight.push_back(weightsForOneNeutrino);
      }//loop over each neutrino
    /*
    for ( unsigned int inu = 0; inu < fluxlist.size(); ++inu )//each neutrino
      {

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
  REGISTER_WEIGHTCALC(PrimaryHadronFSWeightCalc)
}
