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

    // Need vectors here because we have have several instantiations
    // of the module running depending on matching
    std::vector<std::string> fVertexModuleLabelVec;
    std::vector<std::string> fVtxTrackAssnsModuleLabelVec;
    
    // For Cluster2D we only have one version
    std::string              fCosmicProducerLabel;
    std::string              fCosmicClusterAssnsLabel;

    TH1D*                    fTotNumTracks;
    TH1D*                    fTotNumClusters;
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
    fTotNumClusters = tfs->make<TH1D>("TotNumClusters", ";# Clusters", 50, 0., 1000);
}

void TPCNeutrinoIDFilter::endJob()
{
    
}

void TPCNeutrinoIDFilter::beginRun(const art::Run& run)
{
    
}

void TPCNeutrinoIDFilter::reconfigure(fhicl::ParameterSet const& pset)
{
    fVertexModuleLabelVec        = pset.get< std::vector<std::string> >("VertexModuleLabelVec",       std::vector<std::string>() ={"pandoraNu"});
    fVtxTrackAssnsModuleLabelVec = pset.get< std::vector<std::string> >("VtxTrackAssnModuleLabelVec", std::vector<std::string>() ={"neutrinoID"});
    
    if (fVertexModuleLabelVec.size() != fVtxTrackAssnsModuleLabelVec.size())
    {
        mf::LogError("TPCNeutrinoIDFilter") << "Mismatch between string vector lengths input from fhicl!" << std::endl;
    }
    
    // For cluster 2D approach
    fCosmicProducerLabel         = pset.get< std::string > ("Cluster2DCosmicProducerLabel", "ccclustertag");
    fCosmicClusterAssnsLabel     = pset.get< std::string > ("Cluster2DCosmicClusterAssns",  "cluster2D");
    
    return;
}

bool TPCNeutrinoIDFilter::filter(art::Event& event)
{
    bool pass = false;

    // In principle we can have several producers running over various configurations of vertices and tracks.
    // The output associations we want to check are then encapsuated in the input vectors of strings
    // So the outer loop is over the indices
    for(size_t assnIdx = 0; assnIdx < fVertexModuleLabelVec.size(); assnIdx++)
    {
        // Recover a handle to the collection of vertices
        art::Handle< std::vector<recob::Vertex> > vertexVecHandle;
        event.getByLabel(fVertexModuleLabelVec[assnIdx], vertexVecHandle);
    
        if (vertexVecHandle.isValid())
        {
            // Recover associations relating vertices and tracks
            art::FindManyP<recob::Track> vertexTrackAssns(vertexVecHandle, event, fVtxTrackAssnsModuleLabelVec[assnIdx]);
        
            // First check that we have something
            if (vertexTrackAssns.isValid() && vertexTrackAssns.size() > 0)
            {
                // Look for the first valid association
                for (size_t vtxIdx = 0; vtxIdx < vertexVecHandle->size(); vtxIdx++)
                {
                    if (vertexTrackAssns.at(vtxIdx).size() > 0)
                    {
                        // Actually, the fact that there is a non-empty association vector is good enough for now
                        pass = true;
                        break;
                    }
                }
            }
        }
    }
    
    // To check for the 2D cluster ID results we need to get cosmic to cluster associations
    // For that we need to start with the overall cosmic tag producer module
    // Fortunately, we only do this if the above failed...
   // if (!pass)
    //{
        art::Handle<std::vector<anab::CosmicTag>> cosmicVecHandle;
        event.getByLabel(fCosmicProducerLabel, cosmicVecHandle);
    
        if (cosmicVecHandle.isValid())
        {
            art::FindManyP<recob::Cluster> cosmicClusterAssns(cosmicVecHandle, event, fCosmicClusterAssnsLabel);
        
            if (cosmicClusterAssns.isValid() && cosmicClusterAssns.size() > 0)
            {
                for(size_t tagIdx = 0; tagIdx < cosmicVecHandle->size(); tagIdx++)
                {
                    if (cosmicClusterAssns.at(tagIdx).size() > 0)
                    {
                        pass = true;
                        break;
                    }
                }
            }
        }
    //}

    return pass;

} // microboone::TPCNeutrinoIDFilter::filter()


DEFINE_ART_MODULE(TPCNeutrinoIDFilter)

}
