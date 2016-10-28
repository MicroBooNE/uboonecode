// Splines method used for pi+

#include "../WeightCalcCreator.h"
#include "../WeightCalc.h"
#include "uboone/EventWeight/IFDHFileTransfer.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "canvas/Utilities/Exception.h"

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
#include "TSpline.h"
#include "TDecompChol.h"

using namespace std;

namespace evwgh {
  class PrimaryHadronSplinesWeightCalc : public WeightCalc
  {
    evwgh::IFDHFileTransfer IFDH;
  public:
    PrimaryHadronSplinesWeightCalc();
    void Configure(fhicl::ParameterSet const& p);
    virtual std::vector<std::vector<double> > GetWeight(art::Event & e);
    TMatrixD myXSecSurfacetmp;
    std::vector<double> momentumBounds;
    std::vector<double> thetaBounds;
    int momentumSplineSampling;
    int thetaSplineSampling;
    TMatrixD covarianceMatrix;
    TMatrixD crossSection;

    void ExternalData(TMatrixD& crossSection, std::vector<double>& momentumBounds,std::vector<double>& thetaBounds,TMatrixD& covarianceMatrix);

  private:
 
    double xsec_splines(TLorentzVector had, TLorentzVector incidentP, TMatrixD myXSecSurfacetmp, std::vector<double> momentumBoundstmp,std::vector<double> thetaBoundstmp);
    std::vector<TMatrixD> GetFakeData();
    std::vector<double> ConvertToVector(TArrayD const* array);
    std::string fGenieModuleLabel;
    std::vector<std::string> fParameter_list;
    float fParameter_sigma;
    int primaryHad;
    int fNmultisims;
    CLHEP::RandGaussQ *GaussRandom;
    std::vector<TMatrixD> fakeData;
    std::string ExternalDataInput;

    DECLARE_WEIGHTCALC(PrimaryHadronSplinesWeightCalc)
  };
  PrimaryHadronSplinesWeightCalc::PrimaryHadronSplinesWeightCalc()
  {
  }

  std::vector<double> PrimaryHadronSplinesWeightCalc::ConvertToVector(TArrayD const* array) {
    std::vector<double> v(array->GetSize());
    std::copy(array->GetArray(), array->GetArray() + array->GetSize(),
	      v.begin());
    return v;
  } // ConvertToVector()

  void PrimaryHadronSplinesWeightCalc::ExternalData(TMatrixD& crossSection, std::vector<double>& momentumBounds,std::vector<double>& thetaBounds,TMatrixD& covarianceMatrix)
  {
    TFile file(ExternalDataInput.c_str(),"READ");

    std::string pname;

    ///reading Pi+ information from external file

    pname = "HARPData/PiPlus/PPCrossSection";

    TMatrixD* PiPluscrossSection = (TMatrixD*) file.Get(pname.c_str());
    if (!PiPluscrossSection) {
      throw art::Exception(art::errors::NotFound)
        << "Could not find parameter '" << pname << "' in file '"
        << file.GetPath() << "'\n";
    }

    pname = "HARPData/PiPlus/PPcovarianceMatrix";

    TMatrixD* PiPluscovarianceMatrix = (TMatrixD*) file.Get(pname.c_str());
    if (!PiPluscovarianceMatrix) {
      throw art::Exception(art::errors::NotFound)
        << "Could not find parameter '" << pname << "' in file '"
        << file.GetPath() << "'\n";
    }

    pname = "HARPData/PiPlus/PPmomentumBoundsArray";

    TArrayD* PiPlusmomentumBoundsArray = (TArrayD*) file.Get(pname.c_str());
    if (!PiPlusmomentumBoundsArray) {
      throw art::Exception(art::errors::NotFound)
        << "Could not find parameter '" << pname << "' in file '"
        << file.GetPath() << "'\n";
    }

    std::vector<double> PiPlusmomentumBounds = PrimaryHadronSplinesWeightCalc::ConvertToVector(PiPlusmomentumBoundsArray);
    delete PiPlusmomentumBoundsArray;

    pname = "HARPData/PiPlus/PPthetaBoundsArray";

    TArrayD* PiPlusthetaBoundsArray = (TArrayD*) file.Get(pname.c_str());
    if (!PiPlusthetaBoundsArray) {
      throw art::Exception(art::errors::NotFound)
        << "Could not find parameter '" << pname << "' in file '"
        << file.GetPath() << "'\n";
    }

    std::vector<double> PiPlusthetaBounds = PrimaryHadronSplinesWeightCalc::ConvertToVector(PiPlusthetaBoundsArray);
    delete PiPlusthetaBoundsArray;

    if(primaryHad == 211){
 
      crossSection.ResizeTo(13,6);///*** do it more generic//read #rows #columns
      covarianceMatrix.ResizeTo(78,78);

      crossSection = *PiPluscrossSection;
      momentumBounds = PiPlusmomentumBounds;
      thetaBounds = PiPlusthetaBounds;
      covarianceMatrix = *PiPluscovarianceMatrix;
    }
    
  }/////end ExternalData function

 
  double PrimaryHadronSplinesWeightCalc::xsec_splines(TLorentzVector had, TLorentzVector incidentP, TMatrixD myXSecSurfacetmp, std::vector<double> momentumBoundstmp,std::vector<double> thetaBoundstmp)
  {

    double requestedP = had.P();
    double requestedT = had.Theta();

    std::vector<double> thetaBinCenter;
    std::vector<double> momentumBinCenter;

    for(unsigned int t = 0; t < thetaBoundstmp.size()-1; ++t){
      thetaBinCenter.push_back(0.5*(thetaBoundstmp[t]+thetaBoundstmp[t+1]));
    }      

    for(unsigned int p = 0; p < momentumBoundstmp.size()-1; ++p){
      momentumBinCenter.push_back(0.5*(momentumBoundstmp[p]+momentumBoundstmp[p+1]));
    }

    //step 1: spline the xsec arg across constant theta bins
    std::vector<TSpline3> thetaSplines;

    for(unsigned int t = 0; t < thetaBinCenter.size(); ++t)
      {
	TH1F primaryHistConstT("primaryHistConstT","primaryHistConstT",momentumBinCenter.size()-1,momentumBinCenter.data());

	for(unsigned int p = 0; p < momentumBinCenter.size(); ++p)
	  {
	    primaryHistConstT.SetBinContent(p+1,myXSecSurfacetmp[p][t]);
	  }
	TSpline3 thetaSpline(&primaryHistConstT);
	thetaSplines.push_back(thetaSpline);
      }//loop over theta bins

    //step 2: spline the constant theta splines of step 1 across momentum

    TH1F histConstantP("histConstantP","histConstantP",thetaBinCenter.size()-1,thetaBinCenter.data());
  
    for(unsigned int t = 0; t < thetaSplines.size(); ++t)
      {
	histConstantP.SetBinContent(t+1,thetaSplines[t].Eval(requestedP));

      }
    TSpline3 splineConstantP(&histConstantP);

    return splineConstantP.Eval(requestedT);
  }//xsec_spline function definition


  void PrimaryHadronSplinesWeightCalc::Configure(fhicl::ParameterSet const& p)
  {
    fGenieModuleLabel= p.get< std::string > ("genie_module_label");
    //get configuration for this function
    fhicl::ParameterSet const &pset=p.get<fhicl::ParameterSet> (GetName());
    std::cout << pset.to_string() << std::endl;

    fParameter_list		=   pset.get<std::vector<std::string> >("parameter_list");
    fParameter_sigma		=   pset.get<float>("parameter_sigma");
    primaryHad			=   pset.get<int>("PrimaryHadronGeantCode");
    fNmultisims                 =   pset.get<int>("number_of_multisims");
    std::string dataInput       =   pset.get< std::string >("ExternalData");
    ExternalDataInput = IFDH.fetch(dataInput);
    
    PrimaryHadronSplinesWeightCalc::ExternalData(crossSection, momentumBounds, thetaBounds, covarianceMatrix);
 

  }

  //Cannot use WeightCalc because need to pass vector of TMatrixD/// could be re-write (in WeightCalc but be careful with other methods)
 
  std::vector<TMatrixD> PrimaryHadronSplinesWeightCalc::GetFakeData()
  {
    std::vector<TMatrixD> myFakeData;
    std::vector<double> thetaBinCenter;
    std::vector<double> momentumBinCenter;
    PrimaryHadronSplinesWeightCalc::ExternalData(crossSection, momentumBounds, thetaBounds, covarianceMatrix);

    for(unsigned int t = 0; t < thetaBounds.size(); ++t){
      thetaBinCenter.push_back(0.5*(thetaBounds[t]+thetaBounds[t+1]));
    }      

    for(unsigned int p = 0; p < momentumBounds.size(); ++p){
      momentumBinCenter.push_back(0.5*(momentumBounds[p]+momentumBounds[p+1]));
    }

    //Transform crossSection vector of vectors into a flat TArrayD
	
    int dim = int((thetaBinCenter.size()-1)*(momentumBinCenter.size()-1));
    double cvArray[dim];
    crossSection.GetMatrix2Array(cvArray);

    //test
    /*
    std::cout << "cvArray = "<<std::endl;
    for (int element = 0; element < dim; ++element) std::cout << element <<std::endl;
    */
	
    TDecompChol dc = TDecompChol(covarianceMatrix);
    if(!dc.Decompose())
      {
	throw art::Exception(art::errors::StdException)
	  << "Cannot decompose matrix to generate fake data.";
	return myFakeData;
      }
	
    TMatrixD U = dc.GetU();

    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::RandGaussQ GaussRandom(rng->getEngine(GetName()));

    for(int n = 0; n < fNmultisims; ++n)
      {
	///get some random numbers      
	std::vector<double> rands(dim);
	GaussRandom.fireArray(dim, rands.data());
		
	//make some constrained fake data
	double fakeDataArray[dim];
	for(int col = 0; col < dim; ++col)
	  {
	    double weightFromU = 0.;
	    for(int row = 0; row < col+1; ++row)
	      {
		weightFromU += U(row,col)*rands[row];
	      }
	    fakeDataArray[col] = weightFromU + cvArray[col];
	  }
	
	//Convert the fake data array into a TMatrixD
	TMatrixD oneFakeXSec(momentumBinCenter.size(),thetaBinCenter.size());
		
	int entry = 0;
	for(unsigned int p = 0; p < momentumBinCenter.size(); ++p)
	  {
	    std::vector<double> thetasForThisMomentum;
	    for(unsigned int t = 0; t < thetaBinCenter.size(); ++t)
	      {
		oneFakeXSec[p][t] = fakeDataArray[entry];
		entry++;
	      }
	  }

	myFakeData.push_back(oneFakeXSec);
      }//loop over multisims
	
    return myFakeData;
  }//getFakeData function


  std::vector<std::vector<double> >  PrimaryHadronSplinesWeightCalc::GetWeight(art::Event & e)
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
	if (fluxlist[inu].ftptype != primaryHad){ 
	  for(int it =0; it<fNmultisims; it++) weightsForOneNeutrino.push_back(1.);
	  weight.push_back(weightsForOneNeutrino);
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
	    {
	      std::cout<< "Failed, please add the mass of the requested Primary Hadron."<< primaryHad<<"   ***"<<std::endl;
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

	PrimaryHadronSplinesWeightCalc::ExternalData(crossSection, momentumBounds, thetaBounds, covarianceMatrix);

	double crossSection_at_pt = xsec_splines(hadronVec, incidentPVec, crossSection,momentumBounds, thetaBounds);

      	fakeData = GetFakeData();

	//loop over each multisim
	for(Int_t i = 0; i < fNmultisims; ++i)
	  {
	    double fakeCrossSection_at_pt = xsec_splines(hadronVec, incidentPVec, fakeData[i],momentumBounds, thetaBounds);

	    weightsForOneNeutrino.push_back(fakeCrossSection_at_pt/crossSection_at_pt);
	  }//over number of sims
	
	weight.push_back(weightsForOneNeutrino);

      } //loop over each neutrino

    return weight;
  }

  REGISTER_WEIGHTCALC(PrimaryHadronSplinesWeightCalc)
}
