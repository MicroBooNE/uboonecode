/**
 *  @file   NuMuCCInclusiveAlg.cxx
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
#include "uboone/TPCNeutrinoIDFilter/Algorithms/ModNuMuCCInclusiveAlg.h"

// Framework Includes
#include "art/Framework/Core/FindManyP.h"

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Track.h"
#include "lardata/RecoBase/PFParticle.h"
#include "lardata/RecoBase/Vertex.h"
#include "lardata/RecoBase/OpFlash.h"
#include "lardata/AnalysisBase/CosmicTag.h"

//------------------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace neutrinoid {

ModNuMuCCInclusiveAlg::ModNuMuCCInclusiveAlg(fhicl::ParameterSet const &pset) :
        fMyProducerModule(0),
        fGeometry(lar::providerFrom<geo::Geometry>()),
        fDetector(lar::providerFrom<detinfo::DetectorPropertiesService>())
{
    this->reconfigure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ModNuMuCCInclusiveAlg::~ModNuMuCCInclusiveAlg()
{
}

void ModNuMuCCInclusiveAlg::reconfigure(fhicl::ParameterSet const &inputPset)
{
    // Assume we could be called externally with the top level module's complete parameter set
    const fhicl::ParameterSet& pset = inputPset.get<fhicl::ParameterSet>("NuMuCCInclusiveAlg");

    fTrackModuleLabel        = pset.get<std::string> ("TrackModuleLabel");
    fVertexModuleLabel       = pset.get<std::string> ("VertexModuleLabel");
    fOpFlashModuleLabel      = pset.get<std::string> ("OpFlashModuleLabel");

    fDistToEdgeX             = fGeometry->DetHalfWidth()   - pset.get<double>("DistToEdgeX",   10.);
    fDistToEdgeY             = fGeometry->DetHalfHeight()  - pset.get<double>("DistToEdgeY",   20.);
    fDistToEdgeZ             = fGeometry->DetLength() / 2. - pset.get<double>("DistToEdgeZ",   10.);

    fFlashWidth              = pset.get<double>      ("FlashWidth",                            80.);
    fBeamMin                 = pset.get<double>      ("BeamMin",                              3.55);
    fBeamMax                 = pset.get<double>      ("BeamMax",                              5.15);
    fPEThresh                = pset.get<double>      ("PEThresh",                              50.);
    fMinTrk2VtxDist          = pset.get<double>      ("MinTrk2VtxDist",                         5.);
    fMinTrackLen             = pset.get<double>      ("MinTrackLen",                           75.);

    fDoHists                 = pset.get<bool>        ("FillHistograms",                      false);
}

void ModNuMuCCInclusiveAlg::beginJob(art::ServiceHandle<art::TFileService>& tfs)
{
    if (fDoHists)
    {
        // Define the histograms. Putting semi-colons around the title
        // causes it to be displayed as the x-axis label if the histogram
        // is drawn.
        fNFlashPerEvent   = tfs->make<TH1D>("NFlashEvent", ";Flash/Event",     200,   0.,  200.);
        fFlashPE          = tfs->make<TH1D>("FlashPE",     ";PE",              100,   0.,  100.);
        fFlashTime        = tfs->make<TH1D>("FlashTime",   ";Flash Time(us)",  100, -10.,   30.);
    }

    return;
}

void ModNuMuCCInclusiveAlg::produces(art::EDProducer* owner)
{
    fMyProducerModule = owner;
    fMyProducerModule->produces< art::Assns<recob::Vertex, recob::Track> >();
    fMyProducerModule->produces< art::Assns<recob::Vertex, recob::PFParticle> >();
}


bool ModNuMuCCInclusiveAlg::findNeutrinoCandidates(art::Event & event) const
{
    // Agreed convention is to ALWAYS output to the event store so get a pointer to our collection
    std::unique_ptr<art::Assns<recob::Vertex, recob::Track>>      vertexTrackAssociations(new art::Assns<recob::Vertex, recob::Track>);
    std::unique_ptr<art::Assns<recob::Vertex, recob::PFParticle>> vertexPFParticleAssociations(new art::Assns<recob::Vertex, recob::PFParticle>);

    // Recover the hanles to the vertex and track collections we want to analyze.
    art::Handle<std::vector<recob::Vertex>>  vertexVecHandle;
    art::Handle<std::vector<recob::Track>>   trackVecHandle;
    art::Handle<std::vector<recob::OpFlash>> flashListHandle;

    event.getByLabel(fVertexModuleLabel,    vertexVecHandle);
    event.getByLabel(fTrackModuleLabel,     trackVecHandle);

    //----------------------------------------------------
    std::vector<art::Ptr<recob::OpFlash> > flashlist;

    if (event.getByLabel(fOpFlashModuleLabel,flashListHandle))
        art::fill_ptr_vector(flashlist, flashListHandle);

    // Require valid handles, otherwise nothing to do
    if (vertexVecHandle.isValid() && vertexVecHandle->size() > 0 && trackVecHandle.isValid() && trackVecHandle->size() > 0)
    {
        // Recover associations to PFParticles...
        art::FindManyP<recob::PFParticle> trackToPFPartAssns(trackVecHandle,  event, fTrackModuleLabel);

        //----loop over all the flashes and check if there are flashes within the beam
        //window and above the PE threshold
        const recob::OpFlash* flashPtr(0);
        double                flashmax(0);
        bool                  flashtag(false);

for (const auto& opFlash : flashlist)
        {
            if (opFlash->Time() > fBeamMin && opFlash->Time() < fBeamMax && opFlash->TotalPE() > fPEThresh)
            {
                flashtag = true;

                // Keep track of the largest flash
                if (opFlash->TotalPE() > flashmax)
                {
                    flashPtr = opFlash.get();
                    flashmax = opFlash->TotalPE();
                }
            }

            if (fDoHists)
            {
                fFlashPE->Fill(opFlash->TotalPE(), 1.);
                fFlashTime->Fill(opFlash->Time(), 1.);
            }
        }  //end of loop over all the flashes

        if (fDoHists) fNFlashPerEvent->Fill(flashlist.size(), 1.);

        if (flashtag)
        {
            // We need to keep track of the best combination
            int    VertexCandidate=-1;
            int    TrackCandidate=-1;
            double TrackCandLength = 0;
	    
	    // Get the position of the vertex
	    double vertexXYZ[3];

            // Initialize a vertex and associated track collection
            std::map< int,std::vector<int> > VertexTrackCollection;

            //-----------------------------------------------------------
            for (size_t vertexIdx = 0; vertexIdx < vertexVecHandle->size(); vertexIdx++)
            {
                // Recover art ptr to vertex
                art::Ptr<recob::Vertex> vertex(vertexVecHandle, vertexIdx);

                // Reset track at vertex count
                unsigned int TrackCountAtVertex = 0;

                vertex->XYZ(vertexXYZ);

                TVector3 vertexPos(vertexXYZ[0],vertexXYZ[1],vertexXYZ[2]);


                // For each vertex we loop over all tracks looking for matching pairs
                // The outer loop here, then is over one less than all tracks
                for (size_t trackIdx = 0; trackIdx < trackVecHandle->size(); trackIdx++)
                {
                    // Work with an art Ptr here
                    art::Ptr<recob::Track> track(trackVecHandle,trackIdx);

                    // so we need to get the track direction sorted out.
                    TVector3 trackPos = track->Vertex();
                    TVector3 trackEnd = track->End();

                    // Take the closer end---------------------------------
                    double trackToVertexDist = (trackPos - vertexPos).Mag();

                    if ((trackEnd - vertexPos).Mag() < trackToVertexDist)
                    {
                        trackPos          = track->End();
                        trackEnd          = track->Vertex();
                        trackToVertexDist = (trackPos - vertexPos).Mag();
                    }

                    //--------------------------------------------------------------------------
                    if (trackToVertexDist<fMinTrk2VtxDist)
                    {
                        if ((trackEnd-trackPos).Mag()>TrackCandLength)
                        {
                            // If we are looking at the first track which fulfills the distance to vertex cut
                            if (!TrackCountAtVertex)
                            {
                                // Fill vertex ID into the collection map
                                VertexTrackCollection.insert(std::pair< int,std::vector<int> >(vertexIdx,std::vector<int>()));
                            }

                            // Push back track ID for vertex v
                            VertexTrackCollection.at(vertexIdx).push_back(trackIdx);

                            // Increase track at vertex count
                            TrackCountAtVertex++;
                        }
                    } //end of if track distance is within 5cm
                }  //end of loop over the tracks
            } //end of loop over all the vertex

            // Vertex candidate properties
            VertexCandidate = -1;
            float VertexCosTheta = 0.0;

            // Loop over the collection of vertices
	    for (auto const& VtxTrack : VertexTrackCollection)
            {
                // Get vertexID
                int VertexID = VtxTrack.first;

                // Weighted cos theta average
                float WeightedCosTheta = 0.0;

                // Normalization factor
                float NormFactor = 0.0;

                // Loop over all associated track IDs of this vertex
		for (auto const& TrackID : VtxTrack.second)
                {
                    art::Ptr<recob::Track> track(trackVecHandle,TrackID);

                    // Add all track lengths of tracks close to vertex (Normalization factor of the weighted average)
                    NormFactor += track->Length();
                    // Add cos(theta) weighted by track length
                    WeightedCosTheta += track->Length()*cos(track->Theta());
                }// track ID loop

                // Make average
                WeightedCosTheta /= NormFactor;

                // Check for flatest angle (also backwards pointing)
                if (fabs(WeightedCosTheta) > VertexCosTheta)
                {
                    VertexCandidate = VertexID;
                    VertexCosTheta = fabs(WeightedCosTheta);
                }
            }// vertex collection loop
            
            // Create the candidate vertex
            art::Ptr<recob::Vertex> vertex(vertexVecHandle,VertexCandidate);
	    
	    // Fill vertex candidate coordinates
	    vertex->XYZ(vertexXYZ);
	    
            // Check if the wertex candidate is contained and pick the longest track
            if (inFV(vertexXYZ[0],vertexXYZ[1],vertexXYZ[2]))
            {
                // Looping over track IDs of tracks associated with the vertex candidate
		for (auto const& TrackID : VertexTrackCollection.find(VertexCandidate)->second)
                {
                    art::Ptr<recob::Track> track(trackVecHandle,TrackID);
		    
                    // Check for if track is longer
                    if (track->Length() > TrackCandLength)
                    {
			// Pick the numbers of the longest track
                        TrackCandidate = TrackID;
                        TrackCandLength = track->Length();
                    }
                }
            }  //end of if the vertex is contained

            if (TrackCandidate > -1)
            {
		// Create the candidate track
                art::Ptr<recob::Track>  track(trackVecHandle,TrackCandidate);
		
		// Check if track and flash are matched
                bool trackFlashFlag = FlashTrackDist(flashPtr->ZCenter(), track->Vertex().z(), track->End().z()) < fFlashWidth;
		
		// Check if track is contained in FV
		bool trackContainedFlag = inFV(track->Vertex().x(), track->Vertex().y(), track->Vertex().z());
		trackContainedFlag &= inFV(track->End().x(), track->End().y(), track->End().z());

                // Check to see if we think we have a candidate
                if (TrackCandLength>fMinTrackLen && trackFlashFlag && trackContainedFlag)
                {
                    // Make an association between the best vertex and the matching tracks
                    util::CreateAssn(*fMyProducerModule, event, track, vertex, *vertexTrackAssociations);

                    // Find the associated PFParticle
                    std::vector<art::Ptr<recob::PFParticle>> pfParticleVec = trackToPFPartAssns.at(track.key());

                    if (!pfParticleVec.empty())
                    {
                        util::CreateAssn(*fMyProducerModule, event, pfParticleVec[0], vertex, *vertexPFParticleAssociations);
                    }
                }
            }
        }  //end of if flag
    }

    // Add associations to event.
    event.put(std::move(vertexTrackAssociations));
    event.put(std::move(vertexPFParticleAssociations));

    return true;
}

bool ModNuMuCCInclusiveAlg::inFV(double x, double y, double z) const
{
    double distInX = x - fGeometry->DetHalfWidth();
    double distInY = y;
    double distInZ = z - 0.5 * fGeometry->DetLength();

    if (fabs(distInX) < fDistToEdgeX && fabs(distInY) < fDistToEdgeY && fabs(distInZ) < fDistToEdgeZ) return true;

    return false;
}

//This function returns the distance between a flash and
//a track (in one dimension, here used only for z direction)
double ModNuMuCCInclusiveAlg::FlashTrackDist(double flash, double start, double end) const
{
    if (end >= start) {
        if (flash < end && flash > start) return 0;
        else return std::min(fabs(flash-start), fabs(flash-end));
    }
    else {
        if (flash > end && flash < start) return 0;
        else return std::min(fabs(flash-start), fabs(flash-end));
    }
}

} // namespace
