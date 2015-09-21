/**
 * @file    TPCNeutrinoIDFilter_module.cc
 * @brief   Module to filter neutrino candidate events based on TPC topology
 * @authors aschu@fnal.gov
 * 
 ******************************************************************************/

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Core/FindManyP.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Utilities/LArProperties.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"
#include "SimpleTypesAndConstants/geo_types.h"

// LArSoft data definitions
#include "RecoBase/Track.h"
#include "RecoBase/Shower.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "RecoBase/EndPoint2D.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/OpFlash.h"
#include "AnalysisBase/CosmicTag.h"
#include "AnalysisBase/FlashMatch.h"

// root
#include "TH1D.h"


namespace TPCNeutrinoIDFilter
{


class TPCNeutrinoIDFilter : public art::EDFilter {

public:

    explicit TPCNeutrinoIDFilter(fhicl::ParameterSet const&);
    virtual ~TPCNeutrinoIDFilter();
    
    // Called when job begins for definitions of histograms/tuples/etc.
    void beginJob();
    
    // Called when job completes to deal with output of stuff from beginJob
    void endJob();
    
    // Recover information from the start of a run (if processing across runs)
    void beginRun(const art::Run&);
    
    // Allow for fhicl parameters to possibly change during processing...
    void reconfigure(fhicl::ParameterSet const&);
    
    // The actual method for filtering the events
    bool filter(art::Event&);

private:

    std::string fVertexModuleLabel;
    std::string fVtxTrackAssnsModuleLabel;

    TH1D*       fTotNumTracks;
};


TPCNeutrinoIDFilter::TPCNeutrinoIDFilter(fhicl::ParameterSet const& pset) 
{
    this->reconfigure(pset);
}
    
TPCNeutrinoIDFilter::~TPCNeutrinoIDFilter()
{
}
    
void TPCNeutrinoIDFilter::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs;
    
    fTotNumTracks = tfs->make<TH1D>("TotNumTracks", ";# Tracks", 50, 0., 50.);
}

void TPCNeutrinoIDFilter::endJob()
{
    
}

void TPCNeutrinoIDFilter::beginRun(const art::Run& run)
{
    
}

void TPCNeutrinoIDFilter::reconfigure(fhicl::ParameterSet const& pset)
{
    fVertexModuleLabel        = pset.get< std::string >("VertexModuleLabel",       "pandoraNu");
    fVtxTrackAssnsModuleLabel = pset.get< std::string >("VtxTrackAssnModuleLabel", "vertextrackpair");
    
    return;
}

bool TPCNeutrinoIDFilter::filter(art::Event& event)
{
    bool pass = false;

    // Recover a handle to the collection of vertices
    art::Handle< std::vector<recob::Vertex> > vertexVecHandle;
    event.getByLabel(fVertexModuleLabel, vertexVecHandle);
    
    if (vertexVecHandle.isValid())
    {
        // Recover associations relating cosmic tags and track
        art::FindManyP<recob::Track> vertexTrackAssns(vertexVecHandle, event, fVtxTrackAssnsModuleLabel);
        
        // First check that we have something
        if (vertexTrackAssns.isValid() && vertexTrackAssns.size() > 0)
        {
            // Actually, the fact that there is a non-empty association vector is good enough for now
            pass = true;
        }
    }

    return pass;

} // microboone::TPCNeutrinoIDFilter::filter()


DEFINE_ART_MODULE(TPCNeutrinoIDFilter)

}
