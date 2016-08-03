////////////////////////////////////////////////////////////////////////
//
// @file:  TPCNeutrinoID_module.cc
//
// @brief: Producer module to scan for neutrino candidates in a given event
//
// Class:       TPCNeutrinoID
// Module Type: producer
// File:        TPCNeutrinoID_module.cc
//
//              This is a producer module for scanning the end of reconstruction
//              for neutrino candidates. This serves primarily as a shell for
//              individual algorithms to perform the work, however it does output
//              a standard list of collections based on the output of the algorithms
//
// Configuration parameters:
//
// NeutrinoIDAlgName  - name of algorithm to instantiate for scanning
//
//
// @authors  usher@slac.stanford.edu
// Collating done by Tracy Usher (usher@slac.stanford.edu) based on work
// done by the Neutrino ID task force
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"

// Includes for the interface to our algorithms and their creator
#include "uboone/TPCNeutrinoIDFilter/Algorithms/NeutrinoIDAlgFactory.h"

class TPCNeutrinoID : public art::EDProducer
{
public:

    // Copnstructors, destructor.
    explicit TPCNeutrinoID(fhicl::ParameterSet const & pset);
    virtual ~TPCNeutrinoID();

    // Overrides.
    virtual void reconfigure(fhicl::ParameterSet const & pset);
    virtual void produce(art::Event & e);
    virtual void beginJob();
    virtual void endJob();

private:

    // Fcl parameters.
    std::string                                      fNeutrinoIDAlg;        ///< Algorithm used to do the work

    // Statistics.
    int                                              fNumEvent;             ///< Number of events seen.
    
    // Pointer to the algorithm to do the work
    std::unique_ptr< neutrinoid::NeutrinoIDAlgBase > fNeutrinoIDPtr;        ///< Algorithm to to the work
    
};

DEFINE_ART_MODULE(TPCNeutrinoID)

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
TPCNeutrinoID::TPCNeutrinoID(fhicl::ParameterSet const & pset) :
                      fNumEvent(0)
{
    // Instantiate the algorithm we'll use to do the work
    fNeutrinoIDPtr = neutrinoid::NeutrinoIDAlgFactory().MakeNeutrinoIDAlg(pset);
    
    // "Set" fhicl parameters this module uses
    reconfigure(pset);
    
    // We are a producer, say so here
    //produces<std::vector<recob::Track> >();
    fNeutrinoIDPtr->produces(this);

    // Report.
    mf::LogInfo("TPCNeutrinoID") << "TPCNeutrinoID instantiated\n";
}

//----------------------------------------------------------------------------
/// Destructor.
TPCNeutrinoID::~TPCNeutrinoID()
{
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void TPCNeutrinoID::reconfigure(fhicl::ParameterSet const & pset)
{
    fNeutrinoIDAlg = pset.get<std::string>("NeutrinoIDAlgName", "TrackPairPlusVertexAlg");
}

//----------------------------------------------------------------------------
/// Begin job method.
void TPCNeutrinoID::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs;

    fNeutrinoIDPtr->beginJob(tfs);
}

//----------------------------------------------------------------------------
/// Produce method.
///
/// Arguments:
///
/// evt - Art event.
///
/// This is the primary method.
///
void TPCNeutrinoID::produce(art::Event & event)
{
    ++fNumEvent;
    
    // Call the algorithm to do the work
    fNeutrinoIDPtr->findNeutrinoCandidates(event);
    
    return;
}

//----------------------------------------------------------------------------
/// End job method.
void TPCNeutrinoID::endJob()
{
    mf::LogInfo("TPCNeutrinoID") << "Looked at " << fNumEvent << " events" << std::endl;
}
