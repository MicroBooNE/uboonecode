////////////////////////////////////////////////////////////////////////
// Class:       ToyOneShowerGen
// Module Type: producer
// File:        ToyOneShowerGen_module.cc
//
// Generated at Mon Aug 11 14:14:59 2014 by Kazuhiro Terao using artmod
// from cetpkgsupport v1_05_04.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include <memory>

#include "CLHEP/Random/RandFlat.h"

#include "TRandom.h"
#include "TLorentzVector.h"
#include "TF1.h"

#include "Geometry/Geometry.h"
#include "SummaryData/RunData.h"
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"
#include "Simulation/sim.h"

class ToyOneShowerGen;

class ToyOneShowerGen : public art::EDProducer {
public:
  explicit ToyOneShowerGen(fhicl::ParameterSet const & p);
  virtual ~ToyOneShowerGen();

  void produce(art::Event & e) override;
  void beginRun(art::Run & run) override;
  std::vector<double> GetXYZDirection(double uz=0);
  std::vector<double> GetXYZPosition();

private:

  TF1 *fShapeEnergy;
  TF1 *fShapeTheta;
  CLHEP::RandFlat *fFlatRandom;
  size_t fNEvents;
  double fMass;
  double fTime;
  int    fPDGCode;
};


ToyOneShowerGen::ToyOneShowerGen(fhicl::ParameterSet const & p)
  : fShapeEnergy(nullptr), fShapeTheta(nullptr), fFlatRandom(nullptr)
{
  gRandom->SetSeed(0);
  TRandom3 r(0);
  produces< std::vector<simb::MCTruth>   >();
  produces< sumdata::RunData, art::InRun >();

  // Mass
  fMass = p.get<double>("Mass");
  
  // Time
  fTime = p.get<double>("Time");

  // PDGCode
  fPDGCode = p.get<int>("PDGCode");

  //
  // Energy distribution
  //
  
  // Read-in formula
  fShapeEnergy = new TF1("fShapeEnergy",
			 (p.get<std::string>("EnergyShapeFormula")).c_str(),
			 p.get<double>("EnergyLowerBound"),
			 p.get<double>("EnergyUpperBound"));
  // Check formula
  if(!fShapeEnergy)

    throw cet::exception(__PRETTY_FUNCTION__) 
      << "Failed making energy spectrum shape from provided formula!" << std::endl;

  // Read-in parameter values
  std::vector<double> shape_par_values = p.get<std::vector<double> >("EnergyShapeParameters");

  // Check parameter values
  if((int)(shape_par_values.size()) != fShapeEnergy->GetNpar())

    throw cet::exception(__PRETTY_FUNCTION__)
      << "Number of parameters provided to EnergyShapeFormula does not match with the formula!" << std::endl;

  // Set parameter values
  for(size_t i=0; i<shape_par_values.size(); ++i)

    fShapeEnergy->SetParameter(i,shape_par_values.at(i));

  //
  // Angular distribution
  //
  
  // Read-in formula
  fShapeTheta = new TF1("fShapeTheta",
			(p.get<std::string>("ThetaShapeFormula")).c_str(),
			p.get<double>("ThetaLowerBound"),
			p.get<double>("ThetaUpperBound"));
  // Check formula
  if(!fShapeTheta)

    throw cet::exception(__PRETTY_FUNCTION__) 
      << "Failed making energy spectrum shape from provided formula!" << std::endl;

  // Read-in parameter values
  shape_par_values = p.get<std::vector<double> >("ThetaShapeParameters");

  // Check parameter values
  if((int)(shape_par_values.size()) != fShapeTheta->GetNpar())

    throw cet::exception(__PRETTY_FUNCTION__)
      << "Number of parameters provided to ThetaShapeFormula does not match with the formula!" << std::endl;

  // Set parameter values
  for(size_t i=0; i<shape_par_values.size(); ++i)

    fShapeTheta->SetParameter(i,shape_par_values.at(i));

  //
  // Random engine initialization
  //
  createEngine(sim::GetRandomNumberSeed());
  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine &engine = rng->getEngine();
  fFlatRandom = new CLHEP::RandFlat(engine);
}

//------------------------------------------------------------------------------
void ToyOneShowerGen::beginRun(art::Run& run)
{
  // grab the geometry object to see what geometry we are using
  art::ServiceHandle<geo::Geometry> geo;

  std::unique_ptr<sumdata::RunData> runData(new sumdata::RunData(geo->DetectorName()));

  run.put(std::move(runData));

  return;
}

ToyOneShowerGen::~ToyOneShowerGen()
{
  delete fShapeTheta;
  delete fShapeEnergy;
  delete fFlatRandom;
}

std::vector<double> ToyOneShowerGen::GetXYZPosition() {

  std::vector<double> pos(3,0);

  //pos.at(0) = fFlatRandom->fire(170.,2390.);
  //pos.at(1) = fFlatRandom->fire(-995.,995.);
  //pos.at(2) = fFlatRandom->fire(170.,10190.);

  pos.at(0) = fFlatRandom->fire(17.,239.);
  pos.at(1) = fFlatRandom->fire(-99.5,99.5);
  pos.at(2) = fFlatRandom->fire(17.,1019.);

  return pos;
}

std::vector<double> ToyOneShowerGen::GetXYZDirection(double uz) {

  std::vector<double> udir(3,0);

  udir.at(2) = uz;

  double leftover = 1. - pow(udir.at(2),2);
  double frac = fFlatRandom->fire(0,leftover);
  
  udir.at(0) = sqrt(frac);
  udir.at(1) = sqrt(leftover-frac);

  double ux_sign = fFlatRandom->fire(0,1);
  double uy_sign = fFlatRandom->fire(0,1);
  if(ux_sign<0.5) udir.at(0) *= -1;
  if(uy_sign<0.5) udir.at(1) *= -1;

  return udir;
}

void ToyOneShowerGen::produce(art::Event & e)
{
  std::unique_ptr< std::vector<simb::MCTruth> > mctArray(new std::vector<simb::MCTruth>);
  double Evis = fShapeEnergy->GetRandom();
  double Uz   = TMath::Cos(fShapeTheta->GetRandom());

  std::vector<double> pos = GetXYZPosition();

  std::vector<double> dir = GetXYZDirection(Uz);

  simb::MCTruth truth;

  TLorentzVector pos_lorentz(pos.at(0), pos.at(1), pos.at(2), fTime);
  TLorentzVector mom_lorentz( dir.at(0) * Evis,
			      dir.at(1) * Evis,
			      dir.at(2) * Evis,
			      sqrt(pow(Evis,2)+pow(fMass,2)));

  simb::MCParticle part(0, fPDGCode, "primary", 0, fMass, 1);

  part.AddTrajectoryPoint(pos_lorentz, mom_lorentz);

  truth.Add(part);

  mctArray->push_back(truth);
  
  e.put(std::move(mctArray));

}

DEFINE_ART_MODULE(ToyOneShowerGen)
