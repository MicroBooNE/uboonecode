////////////////////////////////////////////////////////////////////////
// Class:       CheckLArG4
// Module Type: analyzer
// File:        CheckLArG4_module.cc
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
#include "Simulation/SimChannel.h"
// STL
#include <sstream>
#include <iostream>

class CheckLArG4;

class CheckLArG4 : public art::EDAnalyzer {
public:
  explicit CheckLArG4(fhicl::ParameterSet const & p);
  virtual ~CheckLArG4();

  void analyze(art::Event const & e) override;


private:

  /// Producer name for G4 module
  std::string fG4Label;

};


CheckLArG4::CheckLArG4(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)
{
  // Obtain producer module label for LArG4
  fG4Label = p.get<std::string>("G4Label");
}

CheckLArG4::~CheckLArG4()
{
}

void CheckLArG4::analyze(art::Event const & e)
{

  //
  // Access MCParticle information from LArG4
  //
  art::Handle<std::vector<simb::MCParticle> > mcpHandle;
  e.getByLabel(fG4Label,mcpHandle);

  // Check if data product is found or not
  if(!mcpHandle.isValid())

    std::cout 
      << std::endl
      << Form("\033[93mDid not find MCParticle from label %s\033[00m",fG4Label.c_str())
      << std::endl
      << std::endl;

  else {
    //
    // Found data product. Access each MCParticle information
    //

    // Create string to store event-wise cout information
    std::ostringstream msg;

    // Report # of MCParticles stored by LArG4
    msg << Form("Found \033[93m%zu MCParticle\033[00m data products...", mcpHandle->size()) << std::endl
	<< "Listing MCParticle summary below:" << std::endl
	<< std::endl;

    // Loop over stored MCParticle and print information. This is a list of MCParticle stored by LArG4
    for(auto const& mcp : *mcpHandle) {

      if(mcp.Mother()>0) continue;

      // For each particle, report TrackID, PDGID, Creation Process name, vtx and momentum      
      msg 
	<< "    " << "Track ID   : " << mcp.TrackId() << std::endl
	<< "    " << "PDG ID     : " << mcp.PdgCode() << std::endl
	<< "    " << "Process    : " << mcp.Process().c_str() << std::endl
	<< "    " << "Mom Track  : " << mcp.Mother() << std::endl
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
    
    msg << "\033[95mEvent ends...\033[00m" << std::endl << std::endl;
      
    std::cout << msg.str() << std::endl;

  }

}

DEFINE_ART_MODULE(CheckLArG4)
