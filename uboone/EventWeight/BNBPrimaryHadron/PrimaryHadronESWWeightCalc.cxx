
// extended sanford-wang (Athula) xsec fit used for pi+, pi-, K0_s
//// parameters form pi+ available. No idea from where Athula gets pi- data sets
/// no k0_s found!!

#include "../WeightCalcCreator.h"
#include "../WeightCalc.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandomEngine.h"
#include "nusimdata/SimulationBase/MCFlux.h"

#include "TLorentzVector.h"

#include <iostream>
using namespace std;

namespace evwgh {
  class PrimaryHadronESWWeightCalc : public WeightCalc
  {
     public:
       PrimaryHadronESWWeightCalc();
       void Configure(fhicl::ParameterSet const& p);
       std::vector<std::vector<double> > GetWeight(art::Event & e);
     private:
    //       CLHEP::RandGaussQ *fGaussRandom; ***---> no tiene que estar aqui

    double xsec_esw(TLorentzVector had, TLorentzVector incidentP, std::vector<double> myFitParameters);
    std::vector< std::vector<double> > fakeData;	
    std::string fGenieModuleLabel;
    std::vector<std::string> fParameter_list;
    float fParameter_sigma;
    int fNmultisims;
    int primaryHad;
    std::vector<double> xSecFitParameters;
    std::vector< std::vector<double> > fitCovarianceMatrix;

    DECLARE_WEIGHTCALC(PrimaryHadronESWWeightCalc)
  };
  PrimaryHadronESWWeightCalc::PrimaryHadronESWWeightCalc()
  {
  }

  // Extended Sanford-Wang cross section fit from HARP paper: Eur. Phys. J. C 52, 29-53 (2007)
  double PrimaryHadronESWWeightCalc::xsec_esw(TLorentzVector had, TLorentzVector incidentP, std::vector<double> myFitParameters)
  {
    if(xSecFitParameters.size() != 9)//check number of xSec fit parameters
      {
	throw art::Exception(art::errors::StdException)
	  << "Failed to run, PrimaryHadronESW requires 9 fit parameters";
	return 0;
      }

    double a = myFitParameters[0];//c1
    double b = myFitParameters[1];//c2
    double c = myFitParameters[2];//c3
    double d = myFitParameters[3];//c4
    double e = myFitParameters[4];//c5
    double f = myFitParameters[5];//c6
    double g = myFitParameters[6];//c7
    double h = myFitParameters[7];//c8
    double i = myFitParameters[8];//c9

    double factor1 = had.Theta()*(had.P() - g*incidentP.P()*pow(cos( had.Theta() ),h));
    double factor2 = a - c*pow(had.P(),d)/pow(incidentP.P(),e) - f*factor1;

    double xsec = exp(factor2)*pow(had.P(),b)*(1. - ( had.P()/incidentP.P() ))*pow( (1. + ( had.P()/incidentP.P() )), (i*factor1));

    return xsec;
  }//xsec

  void PrimaryHadronESWWeightCalc::Configure(fhicl::ParameterSet const& p)////esto se introduce desde las fcl files!!!
  {

    fGenieModuleLabel= p.get< std::string > ("genie_module_label");
    //get configuration for this function
    fhicl::ParameterSet const &pset = p.get<fhicl::ParameterSet> (GetName());
    fParameter_list	       	=   pset.get<std::vector<std::string> >("parameter_list");
    fParameter_sigma		=   pset.get<float>("parameter_sigma");
    fNmultisims			=   pset.get<int>("number_of_multisims");
    primaryHad			=   pset.get<int>("PrimaryHadronGeantCode");
    xSecFitParameters		=   pset.get<std::vector<double> >("CrossSecFitParameters");///maybe different name?
    fitCovarianceMatrix	        =   pset.get<std::vector< std::vector<double> > >("CrossSecFitCovarianceMatrix");///maybe different name?
   
    //Prepare random generator
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::RandGaussQ GaussRandom(rng->getEngine(GetName()));
		
    fakeData = WeightCalc::MultiGaussianSmearing(xSecFitParameters, fitCovarianceMatrix, fNmultisims, GaussRandom);

  }

  std::vector<std::vector<double> > PrimaryHadronESWWeightCalc::GetWeight(art::Event & e)
  {
    //calculate weight(s) here 
    std::vector<std::vector<double> > weight;

    //get the MC generator information out of the event       
    //these are all handles to mc information.
    art::Handle< std::vector<simb::MCFlux> > mcFluxHandle;

    //actually go and get the stuff
    e.getByLabel(fGenieModuleLabel,mcFluxHandle);

    std::vector<simb::MCFlux> const& fluxlist = *mcFluxHandle;

    for ( unsigned int inu = 0; inu < fluxlist.size(); ++inu )//each neutrino
      {
	if (fluxlist[inu].ftptype != primaryHad) continue;

    	//determine mass
    	double m = 0; //mass
    	switch (abs(primaryHad))
	  {
	  case 211:
	  m = 0.13957010; //pi+ and pi- mass [GeV/c^2]/// check if we really have the pi- data
	  break;
	  case 310:
	    m = 0.4976140;//k0_s mass [GeV/c^2]/// check where it is the k0 data
	    break;
	  case 321:
	    m = 0.4936770; //kaon+ and kaon- mass [GeV/c^2]//check where it is the k+/- data
	  default:
	    throw art::Exception(art::errors::StdException)
	      << "Failed, please add the mass of the requested Primary Hadron.";
	  }

	double px = fluxlist[inu].ftpx;
	double py = fluxlist[inu].ftpy;
	double pz = fluxlist[inu].ftpz;
	double e  = sqrt(px*px + py*py + pz*pz + m*m);
	TLorentzVector hadronVec;
	hadronVec.SetPxPyPzE(px,py,pz,e);

	//find hadron's incident proton info
	TLorentzVector incidentPVec;
	double px_proton = fluxlist[inu].fbeampx;
	double py_proton = fluxlist[inu].fbeampy;
	double pz_proton = fluxlist[inu].fbeampz;
	double m_proton  = 0.9382720;
	double e_proton  = sqrt(px_proton*px_proton + py_proton*py_proton + pz_proton*pz_proton + m_proton*m_proton);  
	incidentPVec.SetPxPyPzE(px_proton,py_proton,pz_proton,e_proton);

	double crossSection = xsec_esw(hadronVec, incidentPVec, xSecFitParameters);
	std::vector<double> weightsForOneNeutrino;
    	//loop over each multisim
    	for(unsigned int i = 0; i < fakeData.size(); ++i)
	  {
	    double fakeCrossSection = xsec_esw(hadronVec, incidentPVec, fakeData[i]);
	    weightsForOneNeutrino.push_back(fakeCrossSection/crossSection);
	  }//over number of sims
    	weight.push_back(weightsForOneNeutrino);
      }//loop over each neutrino

    return weight;
  }
  REGISTER_WEIGHTCALC(PrimaryHadronESWWeightCalc)
}
