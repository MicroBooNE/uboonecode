/**
 *  @file   TrackPairPlusVertexAlg.cxx
 * 
 *  @brief  Implementation of the Track/Vertex Neutrino ID alg
 *          This module outputs associations between vertices
 *          and tracks that are found to be within the cut value
 *
 *  Original implementation September 20, 2015 by usher@slac.stanford.edu
 *  This is based on work by Anne Schukraft (aschu@fnal.gov) and her
 *  Neutrino ID task force
 */

// The main include
#include "uboone/TPCNeutrinoIDFilter/Algorithms/TrackPairPlusVertexAlg.h"

// Framework Includes
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"

//------------------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace neutrinoid {

TrackPairPlusVertexAlg::TrackPairPlusVertexAlg(fhicl::ParameterSet const &pset) : fMyProducerModule(0)
{
    this->reconfigure(pset);
    
    m_geometry = lar::providerFrom<geo::Geometry>();
    m_detector = lar::providerFrom<detinfo::DetectorPropertiesService>();
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackPairPlusVertexAlg::~TrackPairPlusVertexAlg()
{
}
    
void TrackPairPlusVertexAlg::reconfigure(fhicl::ParameterSet const &pset)
{
    // Assume we could be called externally with the top level module's complete parameter set
    const fhicl::ParameterSet& myPset = pset.get<fhicl::ParameterSet>("TPCTracksPlusVertexAlg");
    
    fTrackModuleLabel        = myPset.get<std::string> ("TrackModuleLabel",  "trackkalmanhit");
    fVertexModuleLabel       = myPset.get<std::string> ("VertexModuleLabel", "pandoraNu");
    fCosmicModuleLabel       = myPset.get<std::string> ("CosmicModuleLabel", "trackKalmanHitTag");
    fCosmicScoreCut          = myPset.get<double>      ("CosmicScoreCut",    0.4);
    fNeutrinoVtxTrackDistCut = myPset.get<double>      ("NuVtxTrackDistCut", 4.5);
    fDoHists                 = myPset.get<bool>        ("FillHistograms",    true);
}
    
void TrackPairPlusVertexAlg::beginJob(art::ServiceHandle<art::TFileService>& tfs)
{
    if (fDoHists)
    {
        // Define the histograms. Putting semi-colons around the title
        // causes it to be displayed as the x-axis label if the histogram
        // is drawn.
        fMaxDistHists     = tfs->make<TH1D>("TriangleMaxDist", "Max distance for each triangle found",            2000, 0, 1000);
        fBestMaxDistHists = tfs->make<TH1D>("TriBestMaxDist",  "Max distance for the best triangle in the event", 2000, 0, 1000);
    }
    
    return;
}
    
void TrackPairPlusVertexAlg::produces(art::EDProducer* owner)
{
    fMyProducerModule = owner;
    fMyProducerModule->produces< art::Assns<recob::Vertex, recob::Track> >();
}

    
bool TrackPairPlusVertexAlg::findNeutrinoCandidates(art::Event & event) const
{
    // Agreed convention is to ALWAYS output to the event store so get a pointer to our collection
    std::unique_ptr<std::vector<recob::Track> > neutrinoTracks(new std::vector<recob::Track>);
    std::unique_ptr<art::Assns<recob::Vertex, recob::Track> > vertexTrackAssociations(new art::Assns<recob::Vertex, recob::Track>);
    
    // Recover the hanles to the vertex and track collections we want to analyze.
    art::Handle< std::vector<recob::Vertex> > vertexVecHandle;
    art::Handle< std::vector<recob::Track> >  trackVecHandle;
    
    event.getByLabel(fVertexModuleLabel, vertexVecHandle);
    event.getByLabel(fTrackModuleLabel,  trackVecHandle);
    
    // Require valid handles, otherwise nothing to do
    if (vertexVecHandle.isValid() && vertexVecHandle->size() > 0 && trackVecHandle.isValid() && trackVecHandle->size() > 0)
    {
        // Recover associations relating cosmic tags and track
        art::FindManyP<anab::CosmicTag> cosmicAssns(trackVecHandle, event, fCosmicModuleLabel);
        
        // We need to keep track of the best combination
        // Can we assign art ptrs? I don't think so...
        std::vector<art::Ptr<recob::Vertex>> bestVertexVec;
        std::vector<art::Ptr<recob::Track> > bestTrackVec;
        double bestDistance(100000.);
        
        // Outer loop is over the vertices in the collection
        for(size_t vertexIdx = 0; vertexIdx < vertexVecHandle->size(); vertexIdx++)
        {
            // Recover art ptr to vertex
            art::Ptr<recob::Vertex> vertex(vertexVecHandle, vertexIdx);
            
            // Get the position of the vertex
            // Ultimately we really want the vertex position in a TVector3 object...
            double vertexXYZ[3];
            
            vertex->XYZ(vertexXYZ);
            
            TVector3 vertexPos(vertexXYZ[0],vertexXYZ[1],vertexXYZ[2]);
            
            // For each vertex we loop over all tracks looking for matching pairs
            // The outer loop here, then is over one less than all tracks
            for(size_t track1Idx = 0; track1Idx < trackVecHandle->size() - 1; track1Idx++)
            {
                // Work with an art Ptr here
                art::Ptr<recob::Track> track1(trackVecHandle,track1Idx);
                
                // Is there a CosmicTag associated with this track?
                // There are other/better ways to handle this but I think this covers worst case scenario
                if (cosmicAssns.isValid() && cosmicAssns.size() > 0)
                {
                    std::vector<art::Ptr<anab::CosmicTag> > cosmicVec = cosmicAssns.at(track1.key());
                    
                    if (!cosmicVec.empty())
                    {
                        art::Ptr<anab::CosmicTag>& cosmicTag(cosmicVec.front());
                        
                        if (cosmicTag->CosmicScore() > fCosmicScoreCut) continue;
                    }
                }
                
                // Currently we have the problem that tracks can be fit in the "wrong" direction
                // so we need to get the track direction sorted out.
                TVector3 track1Pos = track1->Vertex();
                TVector3 track1End = track1->End();
                
                // Take the closer end
                double track1ToVertexDist = (track1Pos - vertexPos).Mag();
                
                if ((track1End - vertexPos).Mag() < track1ToVertexDist)
                {
                    track1Pos          = track1->End();
                    track1End          = track1->Vertex();
                    track1ToVertexDist = (track1Pos - vertexPos).Mag();
                }
                
                // Is there a cut at this point?
                
                // Now loop over rest of tracks looking for best match
                for(size_t track2Idx = track1Idx+1; track2Idx < trackVecHandle->size(); track2Idx++)
                {
                    // Still working the art ptrs
                    art::Ptr<recob::Track> track2(trackVecHandle,track2Idx);
                    
                    // Check cosmic ray tags again
                    if (cosmicAssns.isValid() && cosmicAssns.size() > 0)
                    {
                        std::vector<art::Ptr<anab::CosmicTag> > cosmicVec = cosmicAssns.at(track2.key());
                        
                        if (!cosmicVec.empty())
                        {
                            art::Ptr<anab::CosmicTag>& cosmicTag(cosmicVec.front());
                            
                            if (cosmicTag->CosmicScore() > fCosmicScoreCut) continue;
                        }
                    }
                    
                    // Same dance for closest position
                    TVector3 track2Pos = track2->Vertex();
                    TVector3 track2End = track2->End();
                    
                    // Take the closer end
                    double track2ToVertexDist = (track2Pos - vertexPos).Mag();
                    
                    if ((track2End - vertexPos).Mag() < track2ToVertexDist)
                    {
                        track2Pos          = track2->End();
                        track2End          = track2->Vertex();
                        track2ToVertexDist = (track2Pos - vertexPos).Mag();
                    }
                    
                    // Which distance is larger?
                    double maxDist = std::max(track1ToVertexDist,track2ToVertexDist);
                    
                    // Now also get the distance between the start of the two tracks
                    double track1ToTrack2Dist = (track1Pos - track2Pos).Mag();
                    
                    // is it larger?
                    maxDist = std::max(maxDist,track1ToTrack2Dist);
                    
                    // fill hist if asked
                    if (fDoHists) fMaxDistHists->Fill(maxDist, 1.);
                    
                    // Is this the best?
                    if (maxDist < bestDistance)
                    {
                        double length1 = (track1End - track1Pos).Mag();
                        double length2 = (track2End - track2Pos).Mag();
                        double cos = 0;

                        if(length1 > 10 && length1 > length2) {
                           double length1z = fabs(track1End.z() - track1Pos.z());
                           cos = length1z/length1;
                        }
                        else if(length2 > 10) {
                           double length2z = fabs(track2End.z() - track2Pos.z());
                           cos = length2z/length2;
                        }

                        if(cos > 0.85) {

			   // Clear out the old results
                           bestVertexVec.clear();
                           bestTrackVec.clear();
                        
                           //Now store away
                           bestVertexVec.push_back(vertex);
                           bestTrackVec.push_back(track1);
                           bestTrackVec.push_back(track2);
                           bestDistance = maxDist;
                       }
                    }
                }
            }
        }
        
        // Fill hists if asked
        if (fDoHists) fBestMaxDistHists->Fill(bestDistance, 1.);
        
        // Check to see if we think we have a candidate
        if (bestDistance < fNeutrinoVtxTrackDistCut)
        {
            // Make an association between the best vertex and the matching tracks
            util::CreateAssn(*fMyProducerModule, event, bestVertexVec[0], bestTrackVec, *vertexTrackAssociations);
        }
    }
    
    // Add tracks and associations to event.
    event.put(std::move(vertexTrackAssociations));
    
    return true;
}

} // namespace
