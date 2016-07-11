////////////////////////////////////////////////////////////////////////
/// \file  MuCSMCGen_module.cc
/// \brief Generator for cosmic-ray muons passing through MuCS
///
/// Module to produce cosmic ray MC for the MuCS
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
#include "larsim/RandomUtils/LArSeedService.h"

// larsoft includes
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"
#include "EventGeneratorBase/evgenbase.h"
#include "larcore/Geometry/geo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/SummaryData/RunData.h"


// CRY related include files
#include "CRYSetup.h"
#include "CRYParticle.h"
#include "CRYGenerator.h"
#include "CLHEP/Random/RandEngine.h"
#include "EventGeneratorBase/CRY/CRYHelper.h"


/// A module to check the results from the Monte Carlo generator
class MuCSMCGen : public art::EDProducer {
public:
  explicit MuCSMCGen(fhicl::ParameterSet const& pset);
  virtual ~MuCSMCGen();                        

  void produce(art::Event& evt);  
  //void beginJob();
  void beginRun(art::Run& run);
  void reconfigure(fhicl::ParameterSet const& p);

private:

  CRYSetup*      fSetup; 
  CRYGenerator*  fGen;
  std::vector<double> fTopLayerDims; //<Dimensions of each layer (linearized multi-dim array layers vs dims (x1,x2,y1,y2,z1,z2))
  std::vector<double> fBottomLayerDims; //<Dimensions of each layer (linearized multi-dim array layers vs dims (x1,x2,y1,y2,z1,z2))
  double fSampleTime;//<t0 for the particle (at the top box position)
  double fEnergyThreshold; //<Only keep particles above this energy
  std::string fCRYConfigStr; //<Configuration string that gets passed to CRY
  int fParticlesPerEvent; //<How many particles to keep per event
};


//____________________________________________________________________________
MuCSMCGen::MuCSMCGen(fhicl::ParameterSet const& pset)
{
  // create a default random engine; obtain the random seed from LArSeedService,
  // unless overridden in configuration with key "Seed"
  art::ServiceHandle<sim::LArSeedService>()->createEngine(*this, pset, "Seed");
  
  this->reconfigure(pset);
  
  produces< std::vector<simb::MCTruth> >();
  produces< sumdata::RunData, art::InRun >();
  
  // Find the pointer to the CRY data tables
  std::string crydatadir;
  const char* datapath = getenv("CRYDATAPATH");
  if( datapath != 0) crydatadir = datapath;
  else{
    mf::LogError("CRYHelper") << "no variable CRYDATAPATH set for cry data location, bail";
    exit(0);
  }
    
  // Construct the event generator object
  fSetup = new CRYSetup(fCRYConfigStr, crydatadir);
  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine& engine = rng->getEngine();
  evgb::RNGWrapper<CLHEP::HepRandomEngine>::set(&engine, &CLHEP::HepRandomEngine::flat);
  fSetup->setRandomFunction(evgb::RNGWrapper<CLHEP::HepRandomEngine>::rng);
  fGen = new CRYGenerator(fSetup);  
}

MuCSMCGen::~MuCSMCGen(){
  delete fGen;
  delete fSetup;
}

void MuCSMCGen::beginRun(art::Run& run)
{
  // grab the geometry object to see what geometry we are using
  art::ServiceHandle<geo::Geometry> geo;
  std::unique_ptr<sumdata::RunData> runcol(new sumdata::RunData(geo->DetectorName()));
  run.put(std::move(runcol));
  return;
}
  
//____________________________________________________________________________
void MuCSMCGen::reconfigure(fhicl::ParameterSet const& p){

  art::ServiceHandle<geo::Geometry> geo;
  fTopLayerDims = p.get< std::vector<double> >( "TopLayerDims" );
  fBottomLayerDims = p.get< std::vector<double> >( "BottomLayerDims" );
  fSampleTime=p.get<double>("SampleTime",0);
  fEnergyThreshold=p.get<double>("EnergyThreshold",0.05);
  fCRYConfigStr=p.get<std::string>("CRYConfigStr");
  fParticlesPerEvent=p.get<int>("ParticlesPerEvent",1);
  
  return;
}


//____________________________________________________________________________
void MuCSMCGen::produce(art::Event& evt)
{
  std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);

  simb::MCTruth mctruth;

   // Generator time at start of sample
  int idctr = 0;

  while (idctr<fParticlesPerEvent) {
    std::vector<CRYParticle*> parts;
    fGen->genEvent(&parts);
    for (unsigned int i=0; i<parts.size(); ++i) {
      // Take ownership of the particle from the vector
      std::unique_ptr<CRYParticle> cryp(parts[i]);
      
      // Get the energies of the particles
      double ke = cryp->ke()*1.0E-3; // MeV to GeV conversion
      if (ke<fEnergyThreshold) continue;
      
      double m = 0.; // in GeV
      static TDatabasePDG*  pdgt = TDatabasePDG::Instance();
      TParticlePDG* pdgp = pdgt->GetParticle(cryp->PDGid());
      if (pdgp) m = pdgp->Mass();
      
      double etot = ke + m;
      double ptot = etot*etot-m*m;
      if (ptot>0.0) ptot = sqrt(ptot);
      else          ptot = 0.0;
      
      // Sort out the momentum components. Remember that the NOvA
      // frame has y up and z along the beam. So uvw -> zxy
      double px = ptot * cryp->v();
      double py = ptot * cryp->w();
      double pz = ptot * cryp->u();
      
      // Particle start position. CRY distributes uniformly in x-y
      // plane at fixed z, where z is the vertical direction. This
      // requires some offsets and rotations to put the particles at
      // the surface in the geometry as well as some rotations
      // since the coordinate frame has y up and z along the
      // beam.
      double vx = cryp->y()*100.0 + 0.5*(fTopLayerDims[1]+fTopLayerDims[0]);
      double vy = fTopLayerDims[3];
      double vz = cryp->x()*100.0 + 0.5*(fTopLayerDims[5]+fTopLayerDims[4]);
      double t  = fSampleTime; // seconds

      //project to bottom boxes y position and keep if new position is in bottom box
      if(py==0.) continue; //ignore horizontal particles
      double dt=(fBottomLayerDims[3]-vy)/py;
      double nv[]={vx + dt*px, vy + dt*py, vz + dt*pz};
      
      if( nv[0]<fBottomLayerDims[0] || nv[0]>fBottomLayerDims[1]
          || nv[2]<fBottomLayerDims[4] || nv[2]>fBottomLayerDims[5])
          continue;
      
      // Boiler plate...
      int istatus    =  1;
      int imother1   = evgb::kCosmicRayGenerator;
      std::string primary("primary");
      simb::MCParticle p(idctr,cryp->PDGid(),primary,imother1,m,istatus);
      TLorentzVector pos(vx,vy,vz,t*1e9);// time needs to be in ns to match GENIE, etc
      TLorentzVector mom(px,py,pz,etot);
      p.AddTrajectoryPoint(pos,mom);
      
      mctruth.Add(p);
      ++idctr;
      break;
    } // Loop on particles in event
  } // Loop on events simulated
    
  mctruth.SetOrigin(simb::kCosmicRay);

  truthcol->push_back(mctruth);
  evt.put(std::move(truthcol));

  return;
}// end produce


DEFINE_ART_MODULE(MuCSMCGen)

