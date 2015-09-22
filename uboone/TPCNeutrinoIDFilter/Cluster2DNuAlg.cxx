/**
 *  @file   Cluster2DNuAlg.cxx
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
#include "TPCNeutrinoIDFilter/Cluster2DNuAlg.h"

// Framework Includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Services/Optional/TFileService.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "Utilities/AssociationUtil.h"

#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "AnalysisBase/CosmicTag.h"

//------------------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace neutrinoid {

Cluster2DNuAlg::Cluster2DNuAlg(fhicl::ParameterSet const &pset) : fMyProducerModule(0)
{
    this->reconfigure(pset);
    
    art::ServiceHandle<geo::Geometry>            geometry;
    art::ServiceHandle<util::DetectorProperties> detectorProperties;
    
    m_geometry = &*geometry;
    m_detector = &*detectorProperties;
}

//------------------------------------------------------------------------------------------------------------------------------------------

Cluster2DNuAlg::~Cluster2DNuAlg()
{
}
    
void Cluster2DNuAlg::reconfigure(fhicl::ParameterSet const &pset)
{
    // Assume we could be called externally with the top level module's complete parameter set
    const fhicl::ParameterSet& myPset = pset.get<fhicl::ParameterSet>("TPCCluster2DNuAlg");
    
    fClusterModuleLabel      = myPset.get<std::string> ("ClusterModuleLabel",   "fuzzycluster");
    fCosmicModuleLabel       = myPset.get<std::string> ("CosmicModuleLabel",    "trackKalmanHitTag");
    fPlaneToCheck            = myPset.get<size_t>      ("PlaneToCheck",         2);
    fMinimumHits             = myPset.get<size_t>      ("MinimumHits",          10);
    fMaxCosmicScore          = myPset.get<float>       ("MaxCosmicScore",       0.4);
    fMaximumTick             = myPset.get<float>       ("MaximumTick",          6370);
    fMinimumTick             = myPset.get<float>       ("MinimumTick",          3210);
    fMaximumWire             = myPset.get<float>       ("MaximumWire",          3420);
    fMinimumWire             = myPset.get<float>       ("MinimumWire",          5);
    fMaximumAngle            = myPset.get<float>       ("MaximumAngle",         0.5);
    fMaximumLengthCut        = myPset.get<float>       ("MaximumLengthCut",     200.);
    fMaximumLength           = myPset.get<float>       ("MaximumLength",        500.);
    fMinimumDeltaTicks       = myPset.get<float>       ("MinimumDeltaTicks",    30.);
    fMinCandidateClusters    = myPset.get<size_t>      ("MinCandidateClusters", 2);
    fMaximumDistance         = myPset.get<float>       ("MaximumDistance",      10.);
    fMaximumTime             = myPset.get<float>       ("MaximumTime",          30.);
}
    
void Cluster2DNuAlg::beginJob(art::ServiceHandle<art::TFileService>& tfs) {}
    
void Cluster2DNuAlg::produces(art::EDProducer* owner)
{
    fMyProducerModule = owner;
    fMyProducerModule->produces< art::Assns<anab::CosmicTag, recob::Cluster> >();
}

    
bool Cluster2DNuAlg::findNeutrinoCandidates(art::Event & event) const
{
    // Agreed convention is to ALWAYS output to the event store so get a pointer to our collection
    std::unique_ptr<art::Assns<anab::CosmicTag, recob::Cluster> > cosmicClusterAssociations(new art::Assns<anab::CosmicTag, recob::Cluster>);
    
    // Recover the hanles to the cluster collection we want to analyze.
    art::Handle< std::vector<recob::Cluster> > clusterVecHandle;
    
    event.getByLabel(fClusterModuleLabel, clusterVecHandle);
    
    // Require valid handle, otherwise nothing to do
    if (clusterVecHandle.isValid())
    {
        // Recover associations relating cluster and hits
        art::FindManyP<recob::Hit> clusterHitAssns(clusterVecHandle, event, fClusterModuleLabel);
        
        // Recover associations relating cosmic tags and Cluster
        art::FindManyP<anab::CosmicTag> cosmicAssns(clusterVecHandle, event, fCosmicModuleLabel);
        
        // Make sure valid handles (this can't not happen?)
        if (clusterHitAssns.isValid() && cosmicAssns.isValid())
        {
            // Set up first loop over clusters to find the collection of "good" clusters
            // Keep the good ones in a vector of art ptrs
            std::vector<art::Ptr<recob::Cluster>> goodClusterVec;
        
            // Loop over input clusters
            for(size_t clusterIdx = 0; clusterIdx < clusterVecHandle->size(); clusterIdx++)
            {
                art::Ptr<recob::Cluster> cluster(clusterVecHandle,clusterIdx);
            
                // Make sure we have the right view/plane
                if (cluster->View() != fPlaneToCheck) continue;
            
                // Make sure we have enough hits
                std::vector<art::Ptr<recob::Hit>> clusterHitVec = clusterHitAssns.at(cluster->ID());
            
                if (clusterHitVec.size() < fMinimumHits) continue;
                
                // Check start/stop conditions
                if (cluster->StartTick() < fMinimumTick || cluster->StartTick() > fMaximumTick) continue;
                if (cluster->EndTick()   < fMinimumTick || cluster->EndTick()   > fMaximumTick) continue;
                if (cluster->StartWire() < fMinimumWire || cluster->StartWire() > fMaximumWire) continue;
                if (cluster->EndWire()   < fMinimumWire || cluster->EndWire()   > fMaximumWire) continue;
                
                // length angle conditions
                float deltaWire = fabs(cluster->StartWire() - cluster->EndWire());
                float deltaTick = fabs(cluster->StartTick() - cluster->EndTick());
                
                if (fabs(cluster->StartAngle()) > fMaximumAngle && deltaWire > fMaximumLengthCut) continue;
                if (deltaWire < fMaximumLength && deltaTick < fMinimumDeltaTicks)                 continue;
            
                // Finally! Check cosmic tag status
                std::vector<art::Ptr<anab::CosmicTag> > cosmicVec = cosmicAssns.at(cluster->ID());
                
                if (!cosmicVec.empty())
                {
                    art::Ptr<anab::CosmicTag>& cosmicTag(cosmicVec.front());
                    
                    if (cosmicTag->CosmicScore() >= fMaxCosmicScore) continue;
                }
                
                // If here then we must have a good candidate cluster, store in our interim vector
                goodClusterVec.push_back(cluster);
            }
            
            // Enought clusters to proceed?
            if (goodClusterVec.size() >= fMinCandidateClusters)
            {
                int numVtcs(0);
                
                // Loop over the good clusters
                // Shouldn't this go to one less than the full number of clusters?
                for(size_t clusterIdx = 0; clusterIdx < goodClusterVec.size()-1; clusterIdx++)
                {
                    art::Ptr<recob::Cluster>& cluster = goodClusterVec.at(clusterIdx);
                    
                    // Recover the associated CosmicTag
                    art::Ptr<anab::CosmicTag> cosmicTag(cosmicAssns.at(cluster->ID()).front());
                    
                    // Container to hold the clusters to associate
                    std::vector<art::Ptr<recob::Cluster>> clusterPtrVec;
                    
                    clusterPtrVec.push_back(cluster);
                    
                    float deltaWire = fabs(cluster->StartWire() - cluster->EndWire());
//                    float deltaTick = fabs(cluster->StartTick() - cluster->EndTick());
                    
                    if(deltaWire > fMaximumLength && cluster->StartWire() < cluster->EndWire())
                    {
                        // Starts with StartWire
                        // **** shouldn't this start at clusterIdx + 1 to avoid double counting?
                        for(size_t k2 = clusterIdx+1; k2 < goodClusterVec.size(); k2++)
                        {
                            if(k2 == clusterIdx) continue;
                            
                            art::Ptr<recob::Cluster>& cluster2 = goodClusterVec.at(k2);
                            
                            float deltaWire2 = fabs(cluster2->StartWire() - cluster2->EndWire());
//                            float deltaTick2 = fabs(cluster2->StartTick() - cluster2->EndTick());
                            
                            if(deltaWire2 > deltaWire) continue;
                            
                            if((fabs(cluster->StartWire() - cluster2->StartWire()) <= fMaximumDistance) &&
                               (fabs(cluster->StartTick() - cluster2->StartTick()) <= fMaximumTime)       )
                            {
                                if(cluster2->StartCharge() > cluster->StartCharge())
                                {
                                    float openangle = fabs(cluster->StartAngle() - cluster2->StartAngle());
                                    
                                    if(openangle > 0.1 && openangle < 1.35)
                                    {
                                        numVtcs++;
                                        clusterPtrVec.push_back(cluster2);
/*
                                        Int_t   maxzlen  = fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]);
                                        Float_t difflen  = fabs(fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]) - fabs(cluster_StartWire[cid2] - cluster_EndWire[cid2]));
                                        Float_t minangle = std::min(fabs(cluster_StartAngle[cid1]),fabs(cluster_StartAngle[cid2]));
                                        Float_t charge   = cluster_StartCharge[cid1] + cluster_StartCharge[cid2];
                                        Float_t totalZ   = fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]) + fabs(cluster_StartWire[cid2] - cluster_EndWire[cid2]);
                                        Float_t HitIntensity= cluster_StartCharge[cid2];
                                        Float_t EvtHitIntensity= cluster_Integral[cid2]+cluster_Integral[cid2];
                                        hOpenAngle_nu   -> Fill(openangle);
                                        hStartAngle_nu  -> Fill(minangle);
                                        hStartCharge_nu -> Fill(charge);
                                        hHitIntensity_nu-> Fill(HitIntensity);
                                        hHitEvtIntensity_nu-> Fill(EvtHitIntensity);
                                        hTotalZ_nu      -> Fill(totalZ);
                                        hMaxLen_nu      -> Fill(maxzlen);
                                        hDiffLen_nu     -> Fill(difflen);
                                        hTotalHits_nu   -> Fill(cluster_NHits[cid1] + cluster_NHits[cid2]);
                                        hZLenVAngle_nu  -> Fill(maxzlen,minangle);
                                        hZLenVCharge_nu -> Fill(maxzlen,charge);
                                        myFile << "nu," << ientry << "," << minangle << "," << openangle << "," << charge << "," << maxzlen << "," << totalZ << "," <<difflen <<","<<EvtHitIntensity<<","<<HitIntensity<<std::endl;
 */
                                    }
                                }
                            }
                            
/*
                            else if((fabs(cluster_StartWire[cid1] - cluster_EndWire[cid2]) <= maxDistance) && (fabs(cluster_StartTick[cid1] - cluster_EndTick[cid2]) <= maxTime))
                            {
                                if(cluster_EndCharge[cid2]>cluster_StartCharge[cid1])
                                {
                                    Float_t openangle=fabs(cluster_StartAngle[cid1] - cluster_EndAngle[cid2]);
                                    if(openangle>0.1 && openangle<1.35)
                                    {
                                        numVtcs++;
                                        Int_t   maxzlen  = fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]);
                                        Float_t difflen  = fabs(fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]) - fabs(cluster_StartWire[cid2] - cluster_EndWire[cid2]));
                                        Float_t minangle = std::min(fabs(cluster_StartAngle[cid1]),fabs(cluster_EndAngle[cid2]));
                                        Float_t charge   = cluster_StartCharge[cid1] + cluster_EndCharge[cid2];
                                        Float_t totalZ   = fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]) + fabs(cluster_StartWire[cid2] - cluster_EndWire[cid2]);
                                        Float_t HitIntensity= cluster_EndCharge[cid2];
                                        Float_t EvtHitIntensity= cluster_Integral[cid2]+cluster_Integral[cid2];
                                        hOpenAngle_nu   -> Fill(openangle);
                                        hStartAngle_nu  -> Fill(minangle);
                                        hStartCharge_nu -> Fill(charge);
                                        hHitIntensity_nu-> Fill(HitIntensity);
                                        hHitEvtIntensity_nu-> Fill(EvtHitIntensity);
                                        hTotalZ_nu      -> Fill(totalZ);
                                        hMaxLen_nu      -> Fill(maxzlen);
                                        hDiffLen_nu     -> Fill(difflen);
                                        hTotalHits_nu   -> Fill(cluster_NHits[cid1] + cluster_NHits[cid2]);
                                        hZLenVAngle_nu  -> Fill(maxzlen,minangle);
                                        hZLenVCharge_nu -> Fill(maxzlen,charge);
                                        myFile << "nu," << ientry << "," << minangle << "," << openangle << "," << charge << "," << maxzlen << "," << totalZ << "," <<difflen <<","<<EvtHitIntensity<<","<<HitIntensity<<std::endl;
                                    }
                                }
                            }
*/
                        }
                    } //end startwire
/*
                    else if(fabs(cluster_EndWire[cid1] - cluster_StartWire[cid1]) > maxLength && cluster_StartWire[cid1] > cluster_EndWire[cid1])
                    {
                        // Starts with EndWire
                        for(unsigned int k2 = 0; k2 < goodClusts.size(); ++k2)
                        {
                            if(k2 == k1) continue;
                            int cid2 = goodClusts.at(k2);
                            if(fabs(cluster_StartWire[cid2] - cluster_EndWire[cid2]) > fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1])) continue;
                            if((fabs(cluster_EndWire[cid1] - cluster_StartWire[cid2]) <= maxDistance) && (fabs(cluster_EndTick[cid1] - cluster_StartTick[cid2]) <= maxTime))
                            {
                                if(cluster_StartCharge[cid2]>cluster_EndCharge[cid1])
                                {
                                    Float_t openangle=fabs(cluster_EndAngle[cid1] - cluster_StartAngle[cid2]);
                                    if(openangle>0.1 && openangle<1.35)
                                    {
                                        numVtcs++;
                                        Int_t   maxzlen  = fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]);
                                        Float_t difflen  = fabs(fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]) - fabs(cluster_StartWire[cid2] - cluster_EndWire[cid2]));
                                        Float_t minangle = std::min(fabs(cluster_EndAngle[cid1]),fabs(cluster_StartAngle[cid2]));
                                        Float_t charge   = cluster_EndCharge[cid1] + cluster_StartCharge[cid2];
                                        Float_t totalZ   = fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]) + fabs(cluster_StartWire[cid2] - cluster_EndWire[cid2]);
                                        Float_t HitIntensity  =  cluster_StartCharge[cid2];
                                        Float_t EvtHitIntensity  =  cluster_Integral[cid2]+cluster_Integral[cid2];
                                        hOpenAngle_nu   -> Fill(openangle);
                                        hStartAngle_nu  -> Fill(minangle);
                                        hStartCharge_nu -> Fill(charge);
                                        hHitIntensity_nu-> Fill(HitIntensity);
                                        hHitEvtIntensity_nu-> Fill(EvtHitIntensity);
                                        hTotalZ_nu      -> Fill(totalZ);
                                        hMaxLen_nu      -> Fill(maxzlen);
                                        hDiffLen_nu     -> Fill(difflen);
                                        hTotalHits_nu   -> Fill(cluster_NHits[cid1] + cluster_NHits[cid2]);
                                        hZLenVAngle_nu  -> Fill(maxzlen,minangle);
                                        hZLenVCharge_nu -> Fill(maxzlen,charge);
                                        myFile << "nu," << ientry << "," << minangle << "," << openangle << "," << charge << "," << maxzlen << "," << totalZ << "," <<difflen <<","<<EvtHitIntensity<<","<<HitIntensity<<std::endl;
                                    }
                                }
                            }
                            else if((fabs(cluster_EndWire[cid1] - cluster_EndWire[cid2]) <= maxDistance) && (fabs(cluster_EndTick[cid1] - cluster_EndTick[cid2]) <= maxTime))
                            {
                                if(cluster_EndCharge[cid2]>cluster_EndCharge[cid1]) 
                                {
                                    Float_t openangle=fabs(cluster_EndAngle[cid1] - cluster_EndAngle[cid2]);
                                    if(openangle>0.1 && openangle<1.35)
                                    {
                                        numVtcs++;
                                        Int_t   maxzlen  = fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]);
                                        Float_t difflen  = fabs(fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]) - fabs(cluster_StartWire[cid2] - cluster_EndWire[cid2]));
                                        Float_t minangle = std::min(fabs(cluster_EndAngle[cid1]),fabs(cluster_EndAngle[cid2]));
                                        Float_t charge   = cluster_EndCharge[cid1] + cluster_EndCharge[cid2];
                                        Float_t totalZ   = fabs(cluster_StartWire[cid1] - cluster_EndWire[cid1]) + fabs(cluster_StartWire[cid2] - cluster_EndWire[cid2]);
                                        Float_t HitIntensity= cluster_EndCharge[cid2];
                                        Float_t EvtHitIntensity= cluster_Integral[cid2]+cluster_Integral[cid2];
                                        hOpenAngle_nu   -> Fill(openangle);
                                        hStartAngle_nu  -> Fill(minangle);
                                        hStartCharge_nu -> Fill(charge);
                                        hHitIntensity_nu-> Fill(HitIntensity);
                                        hHitEvtIntensity_nu-> Fill(EvtHitIntensity);
                                        hTotalZ_nu      -> Fill(totalZ);
                                        hMaxLen_nu      -> Fill(maxzlen);
                                        hDiffLen_nu     -> Fill(difflen);
                                        hTotalHits_nu   -> Fill(cluster_NHits[cid1] + cluster_NHits[cid2]);
                                        hZLenVAngle_nu  -> Fill(maxzlen,minangle);
                                        hZLenVCharge_nu -> Fill(maxzlen,charge);
                                        myFile << "nu," << ientry << "," << minangle << "," << openangle << "," << charge << "," << maxzlen << "," << totalZ <<","<<difflen << ","<<EvtHitIntensity<<","<<HitIntensity<<std::endl;
                                    }
                                }
                            }
                        }
                    } //end endwire
*/
                    
                    // Handle the associations as everything related to the original cosmic tag for now
                    // I think we will need a better way going forward...
                    util::CreateAssn(*fMyProducerModule, event, cosmicTag, clusterPtrVec, *cosmicClusterAssociations);
                }
            }
        }
    }
    
    // Add tracks and associations to event.
    event.put(std::move(cosmicClusterAssociations));
    
    return true;
}

} // namespace
