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
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"

#include <tuple>

//------------------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace neutrinoid {
    
enum class TH1DLabels : size_t
{
    NFlashEvent,
    NFlashBeam,
    FlashPE,
    FlashPEgt10,
    FlashTime,
    FlashTimegt10,
    FlashPEBeam,
    FlashTimeBeam,
    FlashPEEmpty,
    FlashTimeEmpty,
    TrackFlash,
    TrackMatches,
    VertexContained,
    TrackMatchGood,
    DistToFlash,
    ProjLength,
    LengthInZ,
    DistToFlashA,
    ProjLengthA,
    LengthInZA,
    DistToVertex,
    DocaToVertex,
    ArcLenToVertex,
    Count
};
    
enum class TH2DLabels : size_t
{
    DocaVsArcLen,
    Count
};

AltNuMuCCInclusiveAlg::AltNuMuCCInclusiveAlg(fhicl::ParameterSet const &pset) :
    fMyProducerModule(0),
    fGeometry(lar::providerFrom<geo::Geometry>()),
    fDetector(lar::providerFrom<detinfo::DetectorPropertiesService>()),
    fClocks(lar::providerFrom<detinfo::DetectorClocksService>())
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
    
    fDistToEdgeX             = fGeometry->DetHalfWidth()   - pset.get<double>("DistToEdgeX",    2.);
    fDistToEdgeY             = fGeometry->DetHalfHeight()  - pset.get<double>("DistToEdgeY",   20.);
    fDistToEdgeZ             = fGeometry->DetLength() / 2. - pset.get<double>("DistToEdgeZ",   10.);
    
    fFlashWidth              = pset.get<double>      ("FlashWidth",                            10.);
    fBeamMin                 = pset.get<double>      ("BeamMin",                              3.55);
    fBeamMax                 = pset.get<double>      ("BeamMax",                              5.15);
    fPEThresh                = pset.get<double>      ("PEThresh",                              50.);
    fMinTrackLen             = pset.get<double>      ("MinTrackLen",                           75.);
    fMaxTrackDoca            = pset.get<double>      ("MinTrackDoca",                           4.);
    fMaxTrackArcLen          = pset.get<double>      ("MinTrackArcLen",                         4.);
    
    fDoHists                 = pset.get<bool>        ("FillHistograms",                      false);
}
    
void AltNuMuCCInclusiveAlg::beginJob(art::ServiceHandle<art::TFileService>& tfs)
{
    if (fDoHists)
    {
        // Define the histograms. Putting semi-colons around the title
        // causes it to be displayed as the x-axis label if the histogram
        // is drawn.
        fTH1DVec.resize(size_t(TH1DLabels::Count));
        
        fTH1DVec[size_t(TH1DLabels::NFlashEvent)]     = tfs->make<TH1D>("NFlashEvent",     ";Flash/Event",     200,   0.,   200.);
        fTH1DVec[size_t(TH1DLabels::NFlashBeam)]      = tfs->make<TH1D>("NFlashBeam",      ";Flash/Beam Win",   10,   0.,    10.);
        fTH1DVec[size_t(TH1DLabels::FlashPE)]         = tfs->make<TH1D>("FlashPE",         ";PE",              100,   0.,   100.);
        fTH1DVec[size_t(TH1DLabels::FlashPEgt10)]     = tfs->make<TH1D>("FlashPEgt10",     ";PE",              100,   0.,   100.);
        fTH1DVec[size_t(TH1DLabels::FlashTime)]       = tfs->make<TH1D>("FlashTime",       ";Flash Time(us)",  100, -10.,    30.);
        fTH1DVec[size_t(TH1DLabels::FlashTimegt10)]   = tfs->make<TH1D>("FlashTimegt10",   ";Flash Time(us)",  100, -10.,    30.);
        
        fTH1DVec[size_t(TH1DLabels::FlashPEBeam)]     = tfs->make<TH1D>("FlashPEBeam",     ";PE",              100,   0.,   100.);
        fTH1DVec[size_t(TH1DLabels::FlashTimeBeam)]   = tfs->make<TH1D>("FlashTimeBeam",   ";Flash Time(us)",  100, -10.,    30.);
        fTH1DVec[size_t(TH1DLabels::FlashPEEmpty)]    = tfs->make<TH1D>("FlashPEEmpty",    ";PE",              100,   0.,   100.);
        fTH1DVec[size_t(TH1DLabels::FlashTimeEmpty)]  = tfs->make<TH1D>("FlashTimeEmpty",  ";Flash Time(us)",  100, -10.,    30.);
        
        fTH1DVec[size_t(TH1DLabels::TrackFlash)]      = tfs->make<TH1D>("TracksFlash",     ";# tracks",         25,   0.,    25.);
        fTH1DVec[size_t(TH1DLabels::TrackMatches)]    = tfs->make<TH1D>("TracksMatched",   ";# tracks",         25,   0.,    25.);
        fTH1DVec[size_t(TH1DLabels::VertexContained)] = tfs->make<TH1D>("VertexContained", ";# vertices",       25,   0.,    25.);
        fTH1DVec[size_t(TH1DLabels::TrackMatchGood)]  = tfs->make<TH1D>("TracksMatGood",   ";# tracks",         25,   0.,    25.);
        fTH1DVec[size_t(TH1DLabels::DistToFlash)]     = tfs->make<TH1D>("DistToFlash",     ";distance",        100,   0.,   100.);
        fTH1DVec[size_t(TH1DLabels::ProjLength)]      = tfs->make<TH1D>("ProjLength",      ";Track Length",    250,   0.,  1000.);
        fTH1DVec[size_t(TH1DLabels::LengthInZ)]       = tfs->make<TH1D>("LengthInZ",       ";Length in Z",     250,   0.,   500.);
        fTH1DVec[size_t(TH1DLabels::DistToFlashA)]    = tfs->make<TH1D>("DistToFlashA",    ";distance",        100,   0.,   100.);
        fTH1DVec[size_t(TH1DLabels::ProjLengthA)]     = tfs->make<TH1D>("ProjLengthA",     ";Track Length",    250,   0.,  1000.);
        fTH1DVec[size_t(TH1DLabels::LengthInZA)]      = tfs->make<TH1D>("LengthInZA",      ";Length in Z",     250,   0.,   500.);
        
        fTH1DVec[size_t(TH1DLabels::DistToVertex)]    = tfs->make<TH1D>("DistToVertex",    ";Dist to Vtx",     100,   0.,    25.);
        fTH1DVec[size_t(TH1DLabels::DocaToVertex)]    = tfs->make<TH1D>("DocaToVertex",    ";Doca to Vtx",     100,   0.,    25.);
        fTH1DVec[size_t(TH1DLabels::ArcLenToVertex)]  = tfs->make<TH1D>("ArcLenToVertex",  ";ArcLen to Vtx",   100,   0.,    25.);
        
        fTH2DVec.resize(size_t(TH2DLabels::Count));
        fTH2DVec[size_t(TH2DLabels::DocaVsArcLen)]    = tfs->make<TH2D>("DocaVsArcLen",    ";Doca;ArcLen",     100,   0.,    25., 100, 0., 25.);
    }
    
    return;
}
    
void AltNuMuCCInclusiveAlg::produces(art::EDProducer* owner)
{
    fMyProducerModule = owner;
    fMyProducerModule->produces< art::Assns<recob::Vertex, recob::Track> >();
    fMyProducerModule->produces< art::Assns<recob::Vertex, recob::PFParticle> >();
}

    
bool AltNuMuCCInclusiveAlg::findNeutrinoCandidates(art::Event & event) const
{
    // Agreed convention is to ALWAYS output to the event store so get a pointer to our collection
    std::unique_ptr<art::Assns<recob::Vertex, recob::Track>>      vertexTrackAssociations(new art::Assns<recob::Vertex, recob::Track>);
    std::unique_ptr<art::Assns<recob::Vertex, recob::PFParticle>> vertexPFParticleAssociations(new art::Assns<recob::Vertex, recob::PFParticle>);
    
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
            
            if (fDoHists)
            {
                fTH1DVec[size_t(TH1DLabels::FlashPE)]->Fill(opFlash->TotalPE(), 1.);
                fTH1DVec[size_t(TH1DLabels::FlashTime)]->Fill(opFlash->Time(), 1.);
                
                if (opFlash->TotalPE() > 10.)
                {
                    fTH1DVec[size_t(TH1DLabels::FlashPEgt10)]->Fill(opFlash->TotalPE(), 1.);
                    fTH1DVec[size_t(TH1DLabels::FlashTimegt10)]->Fill(opFlash->Time(), 1.);
                }
            }
        }  //end of loop over all the flashes
        
        if (fDoHists)
        {
            fTH1DVec[size_t(TH1DLabels::NFlashEvent)]->Fill(flashlist.size(), 1.);
            fTH1DVec[size_t(TH1DLabels::NFlashBeam)]->Fill(flashVec.size(), 1.);
            
            for(const auto& opFlash : flashlist)
            {
                if (!flashVec.empty())
                {
                    fTH1DVec[size_t(TH1DLabels::FlashPEBeam)]->Fill(opFlash->TotalPE(), 1.);
                    fTH1DVec[size_t(TH1DLabels::FlashTimeBeam)]->Fill(opFlash->Time(), 1.);
                }
                else
                {
                    fTH1DVec[size_t(TH1DLabels::FlashPEEmpty)]->Fill(opFlash->TotalPE(), 1.);
                    fTH1DVec[size_t(TH1DLabels::FlashTimeEmpty)]->Fill(opFlash->Time(), 1.);
                }
            }
        }
        
        // install some counters to keep track of things...
        int nTrackFlash(0);
        int nTracksMatched(0);
        int nVerticesContained(0);
        int nTrackMatchGood(0);
        
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
            
            using PFParticleTrackMap = std::map<size_t,size_t>;
            
            PFParticleTrackMap pfParticleTrackMap;
            
            // Loop through the entire PFParticle collection and match any tracks associated to either the primary
            // or its immediate daughters (if a neutrino like PFParticle) to flashes
            for(size_t pfPartIdx = 0; pfPartIdx < pfParticleHandle->size(); pfPartIdx++)
            {
                art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle,pfPartIdx);
                
                // But only examine if this is the primary in the hierarchy
                if (!pfParticle->IsPrimary()) continue;
                
                // Recover vertices associated to this PFParticle
                std::vector<art::Ptr<recob::Vertex>> primaryVertexVec = pfPartToVertexAssns.at(pfParticle.key());
                
                if (primaryVertexVec.empty()) mf::LogDebug("AltNuMuCCInclusiveAlg") << "*****>> Primary PFParticle with no associated vertex! Key: " << pfParticle.key() << std::endl;

                if (!primaryVertexVec.empty())
                {
                    // Start by looking at tracks associated to the Primary PFParticle
                    std::vector<art::Ptr<recob::Track>>  pfPartTrackVec  = pfPartToTrackAssns.at(pfParticle.key());
                    
                    int nAssociatedTracks(0);
                    
                    pfParticleTrackMap[pfParticle.key()] = pfPartTrackVec.size();
                    
                    // Primary track, if there is one, gets special handling
                    if (!pfPartTrackVec.empty())
                    {
                        const art::Ptr<recob::Track>& track = pfPartTrackVec.at(0);
                        size_t                        bestFlashIdx(0);
                        double                        bestDist(1000.);
                        
                        getBestFlashTrackDist(flashVec, track->Vertex().Z(), track->End().Z(), bestFlashIdx, bestDist);
                        
                        nAssociatedTracks++;
                        
                        mf::LogDebug("AltNuMuCCInclusiveAlg") << "*****>> PFParticle primary with track, key: " << pfParticle.key() << ", pdg: " << pfParticle->PdgCode() << std::endl;
                        
                        if (bestDist > -1.)
                            flashMatchTuple.emplace_back(FlashMatch(bestFlashIdx, bestDist, pfParticle.key(), track.key()));
                    }

                    // Require particle to be "neutrino like"
                    if (pfParticle->PdgCode() == 14 || pfParticle->PdgCode() == 12)
                    {
                        // Now loop over the Primary's daughters
                        for(auto& daughterIdx : pfParticle->Daughters())
                        {
                            art::Ptr<recob::PFParticle> daughter(pfParticleHandle, daughterIdx);
                        
                            // recover tracks associated to this daughter
                            pfPartTrackVec = pfPartToTrackAssns.at(daughter.key());
                            
                            pfParticleTrackMap[pfParticle.key()] += pfPartTrackVec.size();
                    
                            // In all likelihood there is but one track associated to each daughter but loop just in case
                            for(const auto& track : pfPartTrackVec)
                            {
                                size_t bestFlashIdx(0);
                                double bestDist(0.);
                            
                                nAssociatedTracks++;
                        
                                getBestFlashTrackDist(flashVec, track->Vertex().Z(), track->End().Z(), bestFlashIdx, bestDist);
                            
                                if (bestDist > -1.)
                                    flashMatchTuple.emplace_back(FlashMatch(bestFlashIdx, bestDist, pfParticle.key(), track.key()));
                            }
                        }
                    }
                    else mf::LogDebug("AltNuMuCCInclusiveAlg") << "*****>> PFParticle not neutrino like: " << pfParticle.key() << ", pdg: " << pfParticle->PdgCode() << std::endl;
                }
            }
            
            // Process possible matches
            if (!flashMatchTuple.empty())
            {
                std::sort(flashMatchTuple.begin(),flashMatchTuple.end(),[](const auto& left, const auto& right){return std::get<1>(left) < std::get<1>(right);});
                
                nTrackFlash = flashMatchTuple.size();
                
                if (fDoHists)
                {
                    for(const auto& tupleVal : flashMatchTuple)
                    {
                        double distToFlash(std::get<1>(tupleVal));
                        
                        art::Ptr<recob::Track> track(trackVecHandle,std::get<3>(tupleVal));
                        
                        double projLength = projectedLength(track.get());
                        double lengthInZ  = fabs((track->Vertex() - track->End()).Z());
                        
                        fTH1DVec[size_t(TH1DLabels::DistToFlashA)]->Fill(std::min(fabs(distToFlash),99.5), 1.);
                        fTH1DVec[size_t(TH1DLabels::ProjLengthA)]->Fill(projLength, 1.);
                        fTH1DVec[size_t(TH1DLabels::LengthInZA)]->Fill(lengthInZ, 1.);
                    }
                }
                
                for(const auto& tupleVal : flashMatchTuple)
                {
                    double distToFlash(std::get<1>(tupleVal));
                    
                    // If the match distance is large then we're done
                    if (distToFlash > fFlashWidth)
                    {
                        mf::LogDebug("AltNuMuCCInclusiveAlg") << "*****>> distance to flash above cut: " << distToFlash << std::endl;
                        break;
                    }
                    
                    nTracksMatched++;

                    art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle,std::get<2>(tupleVal));
                    art::Ptr<recob::Track>      track(trackVecHandle,std::get<3>(tupleVal));
                    
                    // Recover vertices associated to this PFParticle
                    std::vector<art::Ptr<recob::Vertex>> primaryVertexVec = pfPartToVertexAssns.at(pfParticle.key());
                    
                    // Recover the position of the vertex associated to the primary particle
                    double primaryVertexXYZ[3];
                    
                    primaryVertexVec[0]->XYZ(primaryVertexXYZ);
                    
                    TVector3 primaryVertex(primaryVertexXYZ[0],primaryVertexXYZ[1],primaryVertexXYZ[2]);
                    
                    if (!inFV(primaryVertex))
                    {
                        mf::LogDebug("AltNuMuCCInclusiveAlg") << "*****>> Vertex outside fiducial volume" << std::endl;
                        continue;
                    }
                    
                    nVerticesContained++;
                    
                    // so we need to get the track direction sorted out.
                    TVector3 trackPos = track->Vertex();
                    TVector3 trackEnd = track->End();
                    TVector3 trackDir = track->VertexDirection();
                    
                    // Take the closer end---------------------------------
                    double trackToVertexDist = (trackPos - primaryVertex).Mag();
                    
                    if ((trackEnd - primaryVertex).Mag() < trackToVertexDist)
                    {
                        trackPos          = track->End();
                        trackEnd          = track->Vertex();
                        trackDir          = -track->EndDirection();
                        trackToVertexDist = (trackPos - primaryVertex).Mag();
                    }
                    
                    bool   inFidVolStart = inFV(trackPos);
                    bool   inFidVolEnd   = inFV(trackEnd);
                    bool   endOkInY      = endPointOK(trackEnd);
                    bool   trackEndCheck = (pfParticle->Daughters().size() > 1 && endOkInY) || (pfParticle->Daughters().size() == 1 && inFidVolEnd);
//                    bool   inFidVol      = inFidVolStart && inFidVolEnd;
                    double projLength    = projectedLength(track.get());
                    double trackVtxDoca(0.);
                    double trackVtxArcLen(0.);
                    
                    getTrackVertexDCA(primaryVertex, trackPos, trackDir, trackVtxDoca, trackVtxArcLen);
                    
                    mf::LogDebug("AltNuMuCCInclusiveAlg") << "**> PFParticle key: " << pfParticle.key() << ", track key: " << track.key() << ", len: " << projectedLength(track.get()) << ", start: " << inFidVolStart << ", end: " << inFidVolEnd << ", distToFlash: " << distToFlash << ", # dt: " << pfParticleTrackMap[pfParticle.key()] << std::endl;
                    mf::LogDebug("AltNuMuCCInclusiveAlg") << "    Track/Vertex doca: " << trackVtxDoca << ", arcLen: " << trackVtxArcLen << ", trackToVtx: " << trackToVertexDist << std::endl;
                    
                    // Do some histogramming
                    if (fDoHists)
                    {
                        art::Ptr<recob::Track> track(trackVecHandle,std::get<3>(tupleVal));
                        double                 distToFlash(std::get<1>(tupleVal));
                        double                 lengthInZ  = fabs((track->Vertex() - track->End()).Z());
                        
                        fTH1DVec[size_t(TH1DLabels::DistToFlash)]->Fill(std::min(fabs(distToFlash),99.5), 1.);
                        fTH1DVec[size_t(TH1DLabels::ProjLength)]->Fill(projLength, 1.);
                        fTH1DVec[size_t(TH1DLabels::LengthInZ)]->Fill(lengthInZ, 1.);
                        fTH1DVec[size_t(TH1DLabels::DistToVertex)]->Fill(std::min(trackToVertexDist,24.8), 1.);
                        fTH1DVec[size_t(TH1DLabels::DocaToVertex)]->Fill(std::min(trackVtxDoca,24.8), 1.);
                        fTH1DVec[size_t(TH1DLabels::ArcLenToVertex)]->Fill(std::min(trackVtxArcLen,24.8), 1.);
                        fTH2DVec[size_t(TH2DLabels::DocaVsArcLen)]->Fill(std::min(trackVtxDoca,24.8), std::min(trackVtxArcLen,24.8), 1.);
                    }
                    
                    // Require that the track starts and passes to close to the vertex, that it has a minimum length, starts in the fiducial volume
                    if (trackVtxDoca < fMaxTrackDoca && trackVtxArcLen < fMaxTrackArcLen && projLength > fMinTrackLen && inFidVolStart && trackEndCheck)
                    {
                        nTrackMatchGood++;
                        util::CreateAssn(*fMyProducerModule, event, track,      primaryVertexVec[0], *vertexTrackAssociations);
                        util::CreateAssn(*fMyProducerModule, event, pfParticle, primaryVertexVec[0], *vertexPFParticleAssociations);
                    }
                }
            }
        }  //end of if on flashes

        if (fDoHists)
        {
            fTH1DVec[size_t(TH1DLabels::TrackFlash)]->Fill(nTrackFlash, 1.);
            fTH1DVec[size_t(TH1DLabels::TrackMatches)]->Fill(nTracksMatched, 1.);
            fTH1DVec[size_t(TH1DLabels::VertexContained)]->Fill(nVerticesContained, 1.);
            fTH1DVec[size_t(TH1DLabels::TrackMatchGood)]->Fill(nTrackMatchGood, 1.);
        }
    }
    
    // Add associations to event.
    event.put(std::move(vertexTrackAssociations));
    event.put(std::move(vertexPFParticleAssociations));
    
    return true;
}
    
// Length of reconstructed track.
//----------------------------------------------------------------------------
double AltNuMuCCInclusiveAlg::projectedLength(const recob::Track* track) const
{
    double   result(0.);
    TVector3 lastPoint(track->LocationAtPoint(0));
    TVector3 lastDir(track->DirectionAtPoint(0));
    int      n(track->NumberTrajectoryPoints());
    
    for(int i = 1; i < n; ++i)
    {
        const TVector3& newPoint = track->LocationAtPoint(i);
        
        TVector3 lastToNewPoint = newPoint - lastPoint;
        double   arcLenToDoca   = lastDir.Dot(lastToNewPoint);
        
        result    += arcLenToDoca;
        lastPoint  = lastPoint + arcLenToDoca * lastDir;
        lastDir    = track->DirectionAtPoint(i);
    }
    
    return result;
}
    
bool AltNuMuCCInclusiveAlg::inFV(const TVector3& pos) const
{
    double distInX = pos.X() - fGeometry->DetHalfWidth();
    double distInY = pos.Y();
    double distInZ = pos.Z() - 0.5 * fGeometry->DetLength();
    
    if (fabs(distInX) < fDistToEdgeX && fabs(distInY) < fDistToEdgeY && fabs(distInZ) < fDistToEdgeZ) return true;
    
    return false;
}
    
bool AltNuMuCCInclusiveAlg::endPointOK(const TVector3& pos) const
{
    if (fabs(pos.Y()) < fDistToEdgeY) return true;
    
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
    static const double veryLarge(1000.);
    
    double localBestDist(veryLarge);
    
    for(size_t flashIdx = 0; flashIdx < flashVec.size(); flashIdx++)
    {
        const art::Ptr<recob::OpFlash>& flash = flashVec.at(flashIdx);
        
        double dist = FlashTrackDistInZ(flash->ZCenter(), trackStart, trackEnd);
        
        if (dist < localBestDist)
        {
            localBestDist = dist;
            flashIdx = flashIdx;
        }
    }
    
    if (localBestDist < veryLarge) bestDist = localBestDist;
    else                           bestDist = -1.;
    
    return;
}

//This function returns the distance between a flash and
//a track (in one dimension, here used only for z direction)
double AltNuMuCCInclusiveAlg::FlashTrackDistInZ(double flash, double start, double end) const
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
    
void AltNuMuCCInclusiveAlg::getTrackVertexDCA(const TVector3& vertex, const TVector3& trackStart, TVector3& trackDir, double& doca, double& arcLen) const
{
    TVector3 trackDirUnit  = trackDir;
    TVector3 trackToVertex = vertex - trackStart;
    
    trackDirUnit.SetMag(1.);
    
    arcLen = trackToVertex.Dot(trackDirUnit);
    
    TVector3 trackDocaPos = trackStart + arcLen * trackDirUnit;
    
    doca = (vertex - trackDocaPos).Mag();
    
    return;
}

} // namespace
