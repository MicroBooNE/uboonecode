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

    std::string fTrackModuleLabel;
    std::string fVertexModuleLabel;
    std::string fCosmicTaggerAssocLabel;

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
    fTrackModuleLabel       = pset.get< std::string >("TrackModuleLabel",       "trackkalmanhit");
    fVertexModuleLabel      = pset.get< std::string >("VertexModuleLabel",      "pandoraNu");
    fCosmicTaggerAssocLabel = pset.get< std::string >("CosmicTaggerAssocLabel", "trackkalmanhittag");
    
    return;
}

bool TPCNeutrinoIDFilter::filter(art::Event& evt)
{
    bool pass = false;

    // Recover a handle to the collection of vertices
    art::Handle< std::vector<recob::Vertex> > vertexListHandle;
  
    // Recover a handle to the collection of tracks
    art::Handle< std::vector<recob::Track> > trackListHandle;
    
    // Make sure there are tracks and vertices before doing more work
    if (vertexListHandle.isValid() && vertexListHandle->size() > 0 && trackListHandle.isValid() && trackListHandle->size() > 0)
    {
        double trkstartx = 0;
        double trkstarty = 0;
        double trkstartz = 0;
        double trkendx   = 0;
        double trkendy   = 0;
        double trkendz   = 0;
    
        // Recover the associations between the tracks above and cosmic tags
        art::FindManyP<anab::CosmicTag> cosmicAssns(trackListHandle, evt, fCosmicTaggerAssocLabel);

        size_t NTracks = trackListHandle->size();
    
        fTotNumTracks->Fill(NTracks);
    
        for(size_t iTrk=0; iTrk < trackListHandle->size(); ++iTrk) {//loop over tracks
      
            // Recover track pointer from handle
            art::Ptr<recob::Track> ptrack(trackListHandle, iTrk);
            const recob::Track& track = *ptrack;
            
            //Cosmic Tagger information
            float cosmicScore = 0.;
            
            if (cosmicAssns.isValid()) cosmicScore = cosmicAssns.at(track.ID()).front()->CosmicScore();

            TVector3 pos, end;

            int ntraj = track.NumberTrajectoryPoints();
      
            if(ntraj > 0) {
                pos = track.Vertex();
                end = track.End();
            }

            std::cout << "trajectory points " << ntraj << std::endl;
            if(ntraj > 0) {
                trkstartx = pos.X();
                trkstarty = pos.Y();
                trkstartz = pos.Z();
                trkendx = end.X();
                trkendy = end.Y();
                trkendz = end.Z();
            }

            std::cout << "Track cosmic score: " << cosmicScore << std::endl;
            std::cout << trkstartx << "\t" << trkstarty << "\t" << trkstartz << std::endl;
            std::cout << trkendx << "\t" << trkendy << "\t" << trkendz << std::endl;
        
            // Just put this here for now to enable testing of the fhicl
            pass = true;
        }
  }

  return pass;

} // microboone::TPCNeutrinoIDFilter::filter()


DEFINE_ART_MODULE(TPCNeutrinoIDFilter)

}
