////////////////////////////////////////////////////////////////////////
// Class:       CheckGenerator
// Module Type: analyzer
// File:        CheckGenerator_module.cc
//
// Generated at Tue Jun 17 06:24:01 2014 by Kazuhiro Terao using artmod
// from cetpkgsupport v1_05_04.
////////////////////////////////////////////////////////////////////////


// Framework
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

// NuTools
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"

// LArSoft
#include "Simulation/SimChannel.h"
#include "Geometry/Geometry.h"

// STL
#include <sstream>
#include <iostream>

class CheckGenerator;

class CheckGenerator : public art::EDAnalyzer {
public:
  explicit CheckGenerator(fhicl::ParameterSet const & p);
  virtual ~CheckGenerator();

  void analyze(art::Event const & e) override;


private:

  /// Producer name for particle generator module
  std::string fGenLabel;

};


CheckGenerator::CheckGenerator(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)
{
  // Obtain Module label to access data product
  fGenLabel = p.get<std::string>("GenLabel");
}

CheckGenerator::~CheckGenerator()
{
}

void CheckGenerator::analyze(art::Event const & e)
{

  //
  // Access Generator information
  //
  art::Handle<std::vector<simb::MCTruth> > mctHandle;
  e.getByLabel(fGenLabel,mctHandle);

  std::cout
    << "(Run, Sub-Run, Event-ID) = "
    << Form("( %d, %d, %d )",e.id().run(), e.id().subRun(), e.id().event())
    << std::endl;

  // Check if data product is found or not
  if(!mctHandle.isValid())

    std::cout 
      << std::endl
      << Form("\033[93mDid not find MCTruth from label %s\033[00m",fGenLabel.c_str())
      << std::endl
      << std::endl;

  else {
    //
    // Found data product. Access each MCTruth information
    //

    // Create string to store event-wise cout information
    std::ostringstream msg;

    // Report # of MCTruth data products
    msg << Form("Found \033[93m%zu MCTruth\033[00m data products...", mctHandle->size()) << std::endl;

    // Loop over std::vector<simb::MCTruth> vector data product
    for(auto const& mct : *mctHandle) {

      // Report generator type
      switch(mct.Origin()) {
      case ::simb::kUnknown:
	msg << "  One MCTruth from \033[93munknown origin\033[00m..." << std::endl; break;
      case ::simb::kBeamNeutrino:
	msg << "  One MCTruth from \033[93mNu generator\033[00m..." << std::endl; break;
      case ::simb::kCosmicRay:
	msg << "  One MCTruth from \033[93mCRY\033[00m..." << std::endl; break;
      case ::simb::kSuperNovaNeutrino:
	msg << "  One MCTruth from \033[93mSN generator\033[00m..." << std::endl; break;
      case ::simb::kSingleParticle:
	msg << "  One MCTruth from \033[93msingle particle gun\033[00m..." << std::endl; break;
      }

      // Report # of MCParticles stored in this particular MCTruth
      msg << Form("  This MCTruth holds \033[93m%d MC partciles\033[00m ingected into G4 world.",mct.NParticles()) << std::endl
	  << "  Listing MCParticle summary below:" << std::endl
	  << std::endl;

      // Loop over stored MCParticle and print information. This is info used to inject particle into G4 world
      for(size_t i=0; i<(size_t)(mct.NParticles()); ++i) {

	// Obtain MCParticle through reference copy
	const simb::MCParticle& mcp = mct.GetParticle((int)i);

	// For each particle, report TrackID, PDGID, Creation Process name, vtx and momentum
	msg 
	  << "    " << "Track ID   : " << mcp.TrackId() << std::endl
	  << "    " << "PDG ID     : " << mcp.PdgCode() << std::endl
	  << "    " << "Process    : " << mcp.Process().c_str() << std::endl
	  << "    " << "4-vector   : " << Form("(%g, %g, %g,%g)", 
					       mcp.Vx(),
					       mcp.Vy(),
					       mcp.Vz(),
					       mcp.T()) << std::endl
	  << "    " << "4-momentum : " << Form("(%g, %g, %g,%g)", 
					       mcp.Px(),
					       mcp.Py(),
					       mcp.Pz(),
					       mcp.E()) << std::endl
	  << std::endl;
	  
      }

    }

    //
    // Hack! to look at simchannel
    //
    art::Handle<std::vector<sim::SimChannel> > simchArray;
    e.getByLabel("largeant",simchArray);
    if(simchArray.isValid()) {
      art::ServiceHandle<geo::Geometry> geo;
      std::vector<double> sum(geo->Nplanes(),0);
      for(auto const& simch : *simchArray) {

	auto plane = geo->ChannelToWire(simch.Channel())[0].Plane;
	for(auto const& info : simch.TDCIDEMap()) {

	  for(auto const& edep : info.second)

	    sum.at(plane) += edep.numElectrons;
	  
	}
      }
      msg << "Found SimChannel electrons..." << std::endl;
      for(size_t i=0; i<sum.size(); ++i)
	msg << "Plane " << i << " : " << sum.at(i) << std::endl;
      msg<<std::endl;
    }


    // End of accessing data products
    msg << "\033[95mEvent ends...\033[00m" << std::endl << std::endl;

    // Spit out cout
    std::cout << msg.str() << std::endl;
  }

}

DEFINE_ART_MODULE(CheckGenerator)
