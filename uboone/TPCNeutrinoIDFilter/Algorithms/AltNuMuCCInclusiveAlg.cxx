/**
 *  @file   AltNuMuCCInclusiveAlg.cxx
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
#include "uboone/TPCNeutrinoIDFilter/Algorithms/AltNuMuCCInclusiveAlg.h"

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Track.h"
#include "lardata/RecoBase/Vertex.h"
#include "lardata/RecoBase/PFParticle.h"
#include "lardata/RecoBase/OpFlash.h"
#include "lardata/AnalysisBase/CosmicTag.h"

#include <tuple>

//------------------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace neutrinoid {

AltNuMuCCInclusiveAlg::AltNuMuCCInclusiveAlg(fhicl::ParameterSet const &pset) :
    fMyProducerModule(0),
    fGeometry(lar::providerFrom<geo::Geometry>()),
    fDetector(lar::providerFrom<detinfo::DetectorPropertiesService>())
{
    this->reconfigure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

AltNuMuCCInclusiveAlg::~AltNuMuCCInclusiveAlg()
{
}
    
void AltNuMuCCInclusiveAlg::reconfigure(fhicl::ParameterSet const &inputPset)
{
    // Assume we could be called externally with the top level module's complete parameter set
    const fhicl::ParameterSet& pset = inputPset.get<fhicl::ParameterSet>("AltNuMuCCInclusiveAlg");
    
    fPFParticleModuleLabel   = pset.get<std::string> ("PFParticleModuleLabel",         "pandoraNu");
    fTrackModuleLabel        = pset.get<std::string> ("TrackModuleLabel",              "pandoraNu");
    fVertexModuleLabel       = pset.get<std::string> ("VertexModuleLabel",             "pandoraNu");
    fOpFlashModuleLabel      = pset.get<std::string> ("OpFlashModuleLabel",           "opFlashSat");
    
    fDistToEdgeX             = fGeometry->DetHalfWidth()   - pset.get<double>("DistToEdgeX",    0.);
    fDistToEdgeY             = fGeometry->DetHalfHeight()  - pset.get<double>("DistToEdgeY",   20.);
    fDistToEdgeZ             = fGeometry->DetLength() / 2. - pset.get<double>("DistToEdgeZ",   10.);
    
    fFlashWidth              = pset.get<double>      ("FlashWidth",                            10.);
    fBeamMin                 = pset.get<double>      ("BeamMin",                              3.55);
    fBeamMax                 = pset.get<double>      ("BeamMax",                              5.15);
    fPEThresh                = pset.get<double>      ("PEThresh",                              50.);
    fMinTrk2VtxDist          = pset.get<double>      ("MinTrk2VtxDist",                         5.);
    fMinTrackLen             = pset.get<double>      ("MinTrackLen",                           75.);
    
    fDoHists                 = pset.get<bool>        ("FillHistograms",                      false);
}
    
void AltNuMuCCInclusiveAlg::beginJob(art::ServiceHandle<art::TFileService>& tfs)
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
    
void AltNuMuCCInclusiveAlg::produces(art::EDProducer* owner)
{
    fMyProducerModule = owner;
    fMyProducerModule->produces< art::Assns<recob::Vertex, recob::Track> >();
}

    
bool AltNuMuCCInclusiveAlg::findNeutrinoCandidates(art::Event & event) const
{
    // Agreed convention is to ALWAYS output to the event store so get a pointer to our collection
    std::unique_ptr<art::Assns<recob::Vertex, recob::Track> > vertexTrackAssociations(new art::Assns<recob::Vertex, recob::Track>);
    
    // Recover the hanles to the vertex and track collections we want to analyze.
    art::Handle<std::vector<recob::Vertex>>  vertexVecHandle;
    art::Handle<std::vector<recob::Track>>   trackVecHandle;
    art::Handle<std::vector<recob::OpFlash>> flashListHandle;
    art::Handle<std::vector<recob::PFParticle>> pfParticleHandle;
    
    event.getByLabel(fVertexModuleLabel,     vertexVecHandle);
    event.getByLabel(fTrackModuleLabel,      trackVecHandle);
    event.getByLabel(fPFParticleModuleLabel, pfParticleHandle);
    
    //----------------------------------------------------
    std::vector<art::Ptr<recob::OpFlash> > flashlist;
    
    if (event.getByLabel(fOpFlashModuleLabel,flashListHandle))
        art::fill_ptr_vector(flashlist, flashListHandle);
    
    // Require valid handles, otherwise nothing to do
    if (vertexVecHandle.isValid() && vertexVecHandle->size() > 0 && trackVecHandle.isValid() && trackVecHandle->size() > 0 && pfParticleHandle.isValid() && pfParticleHandle->size() > 0)
    {
        //----loop over all the flashes and check if there are flashes within the beam
        //    window and above the PE threshold
        std::vector<art::Ptr<recob::OpFlash>> flashVec;
        
        for(const auto& opFlash : flashlist)
        {
            if (opFlash->Time() > fBeamMin && opFlash->Time() < fBeamMax && opFlash->TotalPE() > fPEThresh)
                flashVec.push_back(opFlash);
        }  //end of loop over all the flashes
        
        // Do we have any "in time" flashes?
        if(!flashVec.empty())
        {
            // Sort the flashes by largest PE down
            std::sort(flashVec.begin(), flashVec.end(), [](const auto& left, const auto& right){return left->TotalPE() > right->TotalPE();});
            
            // Now we go through the PFParticle collection and examine in detail the PFParticle hierarchies input to us
            // To facilitate this we need to recover a few associations
            art::FindManyP<recob::Vertex>     pfPartToVertexAssns(pfParticleHandle, event, fPFParticleModuleLabel);
            art::FindManyP<recob::Track>      pfPartToTrackAssns(pfParticleHandle,  event, fTrackModuleLabel);
            art::FindManyP<recob::PFParticle> vertexToPFPartAssns(vertexVecHandle,  event, fVertexModuleLabel);
            
            // Define a data structure to keep track of the winners
            using FlashMatch      = std::tuple<size_t, double, size_t, size_t>;
            using FlashMatchTuple = std::vector<FlashMatch>;
            
            FlashMatchTuple flashMatchTuple;
            
            // Loop through the entire PFParticle collection...
            for(size_t pfPartIdx = 0; pfPartIdx < pfParticleHandle->size(); pfPartIdx++)
            {
                art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle,pfPartIdx);
                
                // But only examine if this is the primary in the hierarchy
                if (!pfParticle->IsPrimary()) continue;
                
                // Recover vertices associated to this PFParticle
                std::vector<art::Ptr<recob::Vertex>> primaryVertexVec = pfPartToVertexAssns.at(pfParticle.key());

                if (!primaryVertexVec.empty())
                {
                    // Start by looking at tracks associated to the Primary PFParticle
                    std::vector<art::Ptr<recob::Track>>  pfPartTrackVec  = pfPartToTrackAssns.at(pfParticle.key());
                    
                    // Primary track, if there is one, gets special handling
                    if (!pfPartTrackVec.empty())
                    {
                        const art::Ptr<recob::Track>& track = pfPartTrackVec.at(0);
                        size_t                        bestFlashIdx(0);
                        double                        bestDist(1000.);
                        
                        getBestFlashTrackDist(flashVec, track->Vertex().Z(), track->End().Z(), bestFlashIdx, bestDist);
                        
                        flashMatchTuple.emplace_back(FlashMatch(bestFlashIdx, bestDist, pfParticle.key(), track.key()));
                    }
                    
                    // Now loop over the Primary's daughters
                    for(auto& daughterIdx : pfParticle->Daughters())
                    {
                        art::Ptr<recob::PFParticle> daughter(pfParticleHandle, daughterIdx);
                        
                        // recover tracks associated to this daughter
                        pfPartTrackVec = pfPartToTrackAssns.at(daughter.key());
                    
                        // In all likelihood there is but one track associated to each daughter but loop just in case
                        for(const auto& track : pfPartTrackVec)
                        {
                            size_t bestFlashIdx(0);
                            double bestDist(0.);
                        
                            getBestFlashTrackDist(flashVec, track->Vertex().Z(), track->End().Z(), bestFlashIdx, bestDist);
                            
                            flashMatchTuple.emplace_back(FlashMatch(bestFlashIdx, bestDist, pfParticle.key(), track.key()));
                        }
                    }
                }
            }
            
            // Process possible matches
            if (!flashMatchTuple.empty())
            {
                std::sort(flashMatchTuple.begin(),flashMatchTuple.end(),[](const auto& left, const auto& right){return std::get<1>(left) < std::get<1>(right);});
                
                for(const auto& tupleVal : flashMatchTuple)
                {
                    // If the match distance is large then we're done
                    if (std::get<1>(tupleVal) > fFlashWidth) break;

                    art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle,std::get<2>(tupleVal));
                    art::Ptr<recob::Track>      track(trackVecHandle,std::get<3>(tupleVal));
                    
                    // Recover vertices associated to this PFParticle
                    std::vector<art::Ptr<recob::Vertex>> primaryVertexVec = pfPartToVertexAssns.at(pfParticle.key());
                    
                    // Recover the position of the vertex associated to the primary particle
                    double primaryVertexXYZ[3];
                    
                    primaryVertexVec[0]->XYZ(primaryVertexXYZ);
                    
                    TVector3 primaryVertex(primaryVertexXYZ[0],primaryVertexXYZ[1],primaryVertexXYZ[2]);
                    
                    // so we need to get the track direction sorted out.
                    TVector3 trackPos = track->Vertex();
                    TVector3 trackEnd = track->End();
                    
                    // Take the closer end---------------------------------
                    double trackToVertexDist = (trackPos - primaryVertex).Mag();
                    
                    if ((trackEnd - primaryVertex).Mag() < trackToVertexDist)
                    {
                        trackPos          = track->End();
                        trackEnd          = track->Vertex();
                        trackToVertexDist = (trackPos - primaryVertex).Mag();
                    }
                    
                    bool inFidVol = inFV(trackPos) && inFV(trackEnd);
                    
                    if (trackToVertexDist < fMinTrk2VtxDist && inFidVol)
                    {
                        util::CreateAssn(*fMyProducerModule, event, track, primaryVertexVec[0], *vertexTrackAssociations);
                    }
                }
            }
        }  //end of if flag
    }
    
    // Add associations to event.
    event.put(std::move(vertexTrackAssociations));
    
    return true;
}
    
bool AltNuMuCCInclusiveAlg::inFV(const TVector3& pos) const
{
    double distInX = pos.X() - fGeometry->DetHalfWidth();
    double distInY = pos.Y();
    double distInZ = pos.Z() - 0.5 * fGeometry->DetLength();
    
    if (fabs(distInX) < fDistToEdgeX && fabs(distInY) < fDistToEdgeY && fabs(distInZ) < fDistToEdgeZ) return true;
    
    return false;
}
    
int AltNuMuCCInclusiveAlg::traversePFParticleHierarchy(art::Handle<std::vector<recob::PFParticle>>& pfParticleHandle,
                                                       size_t                                       pfParticleIdx,
                                                       const art::FindManyP<recob::Track>&          trackAssns,
                                                       const art::FindManyP<recob::Vertex>&         vertexAssns,
                                                       std::vector<art::Ptr<recob::PFParticle>>&    pfParticleVec,
                                                       std::vector<art::Ptr<recob::Track>>&         trackVec,
                                                       std::vector<art::Ptr<recob::Vertex>>&        vertexVec) const
{
    // So far no daughters...
    int nDaughters(0);
    
    // Get pointer to PFParticle
    art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle,pfParticleIdx);
    
    // Recover tracks/vertices associated to this PFParticle
    std::vector<art::Ptr<recob::Track>>  pfPartTrackVec  = trackAssns.at(pfParticle.key());
    std::vector<art::Ptr<recob::Vertex>> pfPartVertexVec = vertexAssns.at(pfParticle.key());
    
    if (pfPartTrackVec.size() > 0 && pfPartVertexVec.size() > 0)
    {
        pfParticleVec.push_back(pfParticle);
        trackVec.push_back(pfPartTrackVec.at(0));
        vertexVec.push_back(pfPartVertexVec.at(0));
        nDaughters++;
    }
    
    for(auto& daughterIdx : pfParticle->Daughters())
    {
        nDaughters += traversePFParticleHierarchy(pfParticleHandle, daughterIdx, trackAssns, vertexAssns, pfParticleVec, trackVec, vertexVec);
    }
    
    return nDaughters;
}
    
void AltNuMuCCInclusiveAlg::getBestFlashTrackDist(const std::vector<art::Ptr<recob::OpFlash>>& flashVec,
                                                  double                                       trackStart,
                                                  double                                       trackEnd,
                                                  size_t&                                      flashIdx,
                                                  double&                                      bestDist) const
{
    bestDist = 10000.;
    
    for(size_t flashIdx = 0; flashIdx < flashVec.size(); flashIdx++)
    {
        const art::Ptr<recob::OpFlash>& flash = flashVec.at(flashIdx);
        
        double dist = FlashTrackDist(flash->ZCenter(), trackStart, trackEnd);
        
        if (dist < bestDist)
        {
            bestDist = dist;
            flashIdx = flashIdx;
        }
    }
    
    return;
}

//This function returns the distance between a flash and
//a track (in one dimension, here used only for z direction)
double AltNuMuCCInclusiveAlg::FlashTrackDist(double flash, double start, double end) const
{
    double flashStartDist = flash - start;
    double flashEndDist   = flash - end;
    double bestDist       = 0.;

    // if product of the two above is negative (or zero) then track overlaps flash z position
    // if product is positive then start and end are both on same size of flash and need to
    // compute the distance to the closest track end
    if (flashStartDist * flashEndDist > 0.)
        bestDist = std::min(fabs(flashStartDist),fabs(flashEndDist));
    
    return bestDist;
}

} // namespace
