////////////////////////////////////////////////////////////////////////
/// \file  MuCSRecoToMC_module.cc
/// \brief Generate MCTruth muons with trajectory matching MuCSRecoData
///
///
/// \author  matthew.bass@physics.ox.ac.uk
////////////////////////////////////////////////////////////////////////

// ROOT includes
#include "TRandom3.h"
#include "TH1.h"
#include "TH2.h"
#include "TDatabasePDG.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// larsoft includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/EventGeneratorBase/evgenbase.h"

#include "MuCSRecoData.h"

/// A module to check the results from the Monte Carlo generator
class MuCSRecoToMC : public art::EDProducer {
public:
  explicit MuCSRecoToMC(fhicl::ParameterSet const& pset);
  virtual ~MuCSRecoToMC();                        

  void produce(art::Event& evt);  
  //void beginJob();
  //void beginRun(art::Run& run);
  void reconfigure(fhicl::ParameterSet const& p);

private:

  double fSampleTime;//<t0 for the particle (at the top box position) [ns]
  double fMuonKE; //<Produce muons with this energy [GeV]
};


//____________________________________________________________________________
MuCSRecoToMC::MuCSRecoToMC(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);
  
  produces< std::vector<simb::MCTruth> >();
}

MuCSRecoToMC::~MuCSRecoToMC(){
}

  
void MuCSRecoToMC::reconfigure(fhicl::ParameterSet const& p){
  fSampleTime=p.get<double>("SampleTime",0.);
  fMuonKE=p.get<double>("MuonKE",10.);
  
  return;
}

void MuCSRecoToMC::produce(art::Event& evt)
{
  std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);
  simb::MCTruth mctruth;

  //get MuCSRecoData object
  std::vector< art::Handle< std::vector<MuCS::MuCSRecoData> > > colMuCSRecoData;
  evt.getManyByType( colMuCSRecoData );
  
  //if (recodata->at(0).y()!=0){
  if(colMuCSRecoData.size()==0){
    mf::LogInfo("MuCRecoToMC") << "No MuCSRecoData found; skipping event";
  }else if(colMuCSRecoData[0]->at(0).y()!=0){ //if the MuCS Reco Data isn't just empty, proceed
    art::Handle< std::vector<MuCS::MuCSRecoData> > recodata = colMuCSRecoData[0];
    //make a new MCTruth particle and push it onto mctruth
    double m = 0.; // in GeV
    const int pdg=13; //define in a variable in case we want to vary it later
    static TDatabasePDG*  pdgt = TDatabasePDG::Instance();
    TParticlePDG* pdgp = pdgt->GetParticle(pdg);
    if (pdgp) m = pdgp->Mass();
    
    double etot = fMuonKE + m;
    double ptot = etot*etot-m*m;
    if (ptot>0.0) ptot=sqrt(ptot);
    else ptot=0.0;
    
    
    //particle start momentum
    double theta=recodata->at(0).theta();
    double phi=recodata->at(0).phi();
    double px = ptot*(cos(theta))*cos(phi);
    double py = ptot*sin(theta);
    double pz = ptot*(cos(theta))*sin(phi);
    
    //particle start position
    double vx = recodata->at(0).x();
    double vy = recodata->at(0).y();
    double vz = recodata->at(0).z();
    double t  = fSampleTime; // seconds
    
    int istatus    =  1;
    int imother1   = evgb::kCosmicRayGenerator;
    std::string primary("primary");
    simb::MCParticle p(0,pdg,primary,imother1,m,istatus);
    TLorentzVector pos(vx,vy,vz,t);// time needs to be in ns to match GENIE, etc
    TLorentzVector mom(px,py,pz,etot);
    p.AddTrajectoryPoint(pos,mom);
    
    mctruth.Add(p);
  }
    
  mctruth.SetOrigin(simb::kCosmicRay);

  truthcol->push_back(mctruth);
  evt.put(std::move(truthcol));

  return;
}// end produce


DEFINE_ART_MODULE(MuCSRecoToMC)

