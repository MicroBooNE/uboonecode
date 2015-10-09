/**
 *  @file   Cluster2DNuAlg.cxx
 * 
 *  @brief  Implementation of the 2DCluster Neutrino ID alg
 *          This module outputs associations between cosmicTag 
 *          and clusters that are found after all cuts have been 
 *          implemented. 
 *
 *
 *  Original implementation September 20, 2015 by usher@slac.stanford.edu
 *  This is based on work by Anne Schukraft (aschu@fnal.gov) and her
 *  Neutrino ID task force
 *
 *  Updated by Jessica Esquivel (jeesquiv@syr.edu) and Katherine Woodruff (kwoodruf@nmsu.edu)
 */

// The main include
#include "uboone/TPCNeutrinoIDFilter/Cluster2DNuAlg.h"

// ROOT Includes
#include "TMath.h"

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
    fCosmicModuleLabel       = myPset.get<std::string> ("CosmicModuleLabel",    "fuzzyclusterTag");
    fPlaneToCheck            = myPset.get<size_t>      ("PlaneToCheck",         2);
    fMinimumHits             = myPset.get<size_t>      ("MinimumHits",          10);
    fMaxCosmicScore          = myPset.get<float>       ("MaxCosmicScore",       0.4);
    fMaximumAngle            = myPset.get<float>       ("MaximumAngle",         0.5);
    fMaximumLengthCut        = myPset.get<float>       ("MaximumLengthCut",     200.);
    fMaximumMatchedLengthCut = myPset.get<float>       ("MaximumMatchedLengthCut",     100.);
    fMaximumLength           = myPset.get<float>       ("MaximumLength",        500.);
    fMinimumDeltaTicks       = myPset.get<float>       ("MinimumDeltaTicks",    30.);
    fMinimumDeltaWires       = myPset.get<float>       ("MinimumDeltaWires",    30.);
    fMinCandidateClusters    = myPset.get<size_t>      ("MinCandidateClusters", 2);
    fMaximumDistance         = myPset.get<float>       ("MaximumDistance",      10.);
    fMaximumTime             = myPset.get<float>       ("MaximumTime",          30.);
}
    
void Cluster2DNuAlg::beginJob(art::ServiceHandle<art::TFileService>& tfs) {}
    
void Cluster2DNuAlg::produces(art::EDProducer* owner)
{
    fMyProducerModule = owner;
    fMyProducerModule->produces< std::vector<anab::CosmicTag> >();
    fMyProducerModule->produces< std::vector<recob::Cluster> >();
    fMyProducerModule->produces< art::Assns <anab::CosmicTag, recob::Cluster> >();
    //fMyProducerModule->produces< art::Assns <recob::Cluster, recob::Cluster> >();
}

    
bool Cluster2DNuAlg::findNeutrinoCandidates(art::Event & event) const
{
    // Agreed convention is to ALWAYS output to the event store so get a pointer to our collection
    //std::unique_ptr<std::vector<recob::Cluster> > clusterPtrVec(new std::vector<recob::Cluster>);
    //std::unique_ptr<std::vector<anab::CosmicTag> > cosmiccol(new std::vector<anab::CosmicTag>);
    //std::vector<anab::CosmicTag>  cosmicTagPtrVec(new art::Ptr<anab::CosmicTag>);
    //std::vector<recob::Cluster>  clusterPtrVec(new art::Ptr<recob::Cluster>);
    std::unique_ptr<art::Assns<anab::CosmicTag, recob::Cluster> > cosmicClusterAssociations(new art::Assns<anab::CosmicTag, recob::Cluster>);
//    std::unique_ptr<art::Assns<recob::Cluster, recob::Cluster> > cosmicClusterAssociations(new art::Assns<recob::Cluster, recob::Cluster>);
 
    // Recover the handles to the cluster collection we want to analyze.
    art::Handle< std::vector<recob::Cluster> > clusterVecHandle;
    
    std::vector<art::Ptr<recob::Hit>> clusterHitVec;
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
                std::vector<art::Ptr<recob::Hit>> clusterHitVec = clusterHitAssns.at(cluster.key());
            
                if (clusterHitVec.size() < fMinimumHits) continue;
               
                // length angle conditions
                float deltaWire = fabs(cluster->StartWire() - cluster->EndWire());
                float deltaTick = fabs(cluster->StartTick() - cluster->EndTick());
                
                if (fabs(cluster->StartAngle()) > fMaximumAngle && deltaWire > fMaximumLengthCut) continue;
                if (deltaWire < fMinimumDeltaWires && deltaTick < fMinimumDeltaTicks)             continue;
            
                // Finally! Check cosmic tag status
                std::vector<art::Ptr<anab::CosmicTag> > cosmicVec = cosmicAssns.at(cluster.key());
                
                if (!cosmicVec.empty())
                {
                    art::Ptr<anab::CosmicTag>& cosmicTag(cosmicVec.front());
                    //cosmiccol->push_back(cosmicTag->CosmicScore()); 
                    if (cosmicTag->CosmicScore() >= fMaxCosmicScore) continue;
                }
                
                // If here then we must have a good candidate cluster, store in our interim vector
                goodClusterVec.push_back(cluster);
            }
            
            // Enough clusters to proceed?
            if (goodClusterVec.size() >= fMinCandidateClusters)
            {
                // Loop over the good clusters
                float openangle=0;
                for(size_t clusterIdx = 0; clusterIdx < goodClusterVec.size(); clusterIdx++)
                {
                    art::Ptr<recob::Cluster>& cluster = goodClusterVec.at(clusterIdx);

                    // Container to hold the clusters to associate
                    std::vector<art::Ptr<recob::Cluster>> clusterPtrVec;
                    float deltaWire = fabs(cluster->StartWire() - cluster->EndWire());
                    
                    if(deltaWire > fMaximumMatchedLengthCut && cluster->StartWire() < cluster->EndWire())
                    {
                        clusterPtrVec.push_back(cluster);
                        // first cluster is greater than maximum length and  starts with StartWire?
                        for(size_t k2 = 0; k2 < goodClusterVec.size(); k2++)
                        {
                            if(k2 == clusterIdx) continue;
                            
                            art::Ptr<recob::Cluster>& cluster2 = goodClusterVec.at(k2);
                            
                            float deltaWire2 = fabs(cluster2->StartWire() - cluster2->EndWire());
                            
                            if(deltaWire2 >= deltaWire) continue;
                            
                            if((fabs(cluster->StartWire() - cluster2->StartWire()) <= fMaximumDistance) &&
                               (fabs(cluster->StartTick() - cluster2->StartTick()) <= fMaximumTime)       )
                            {   
				// Vertex found if two clusters are within fMaximumDistance and  fMaximumTime
                                if(cluster2->StartCharge() > cluster->StartCharge() || deltaWire > fMaximumLength)
                                {
				    // StartCharge of short cluster must be greater than StartChage of long cluster OR long cluster must be greater than fMaximumLength 
                                    if(cluster2->StartWire()<cluster2->EndWire()) openangle = fabs(cluster->StartAngle() - cluster2->StartAngle());
                                    else
                                    {
                                        if(cluster2->StartAngle()<0) openangle = fabs(cluster->StartAngle() - ((-1.)*TMath::Pi()+cluster2->StartAngle()));
                                        else openangle = fabs(cluster->StartAngle() - (TMath::Pi()+cluster2->StartAngle()));
                                    } 
                                    if(openangle > 0.1 && openangle < 1.57)
                                    {
                                        clusterPtrVec.push_back(cluster2);
                                    }
                                }
                            }
                            else if((fabs(cluster->StartWire() - cluster2->EndWire()) <= fMaximumDistance) && 
                                    (fabs(cluster->StartTick() - cluster2->EndTick()) <= fMaximumTime))
                            {
                                if(cluster2->EndCharge()>cluster->StartCharge() || deltaWire > fMaximumLength)
                                {
                                    if(cluster2->StartWire() < cluster2->EndWire())
                                    { 
                                        if(cluster2->EndAngle()<0) openangle = fabs(cluster->StartAngle() - ((-1.)*TMath::Pi()+cluster2->EndAngle()));
                                        else openangle = fabs(cluster->StartAngle() - (TMath::Pi()+cluster2->EndAngle()));
                                    }
                                    else openangle = fabs( cluster->StartAngle() - cluster2->EndAngle());
                                    if(openangle>0.1 && openangle<1.57)
                                    {
                                        clusterPtrVec.push_back(cluster2);
                                    }
                                }
                            }
                        } 
                    } //end startwire
                    else if(deltaWire > fMaximumMatchedLengthCut && cluster->StartWire() > cluster->EndWire())
                    {
                        clusterPtrVec.push_back(cluster);

                        // Starts with EndWire
                        for(size_t k2 = 0; k2 < goodClusterVec.size(); k2++)
                        {
                            if(k2 == clusterIdx) continue;
                            
                            art::Ptr<recob::Cluster>& cluster2 = goodClusterVec.at(k2);
                            
                            float deltaWire2 = fabs(cluster2->StartWire() - cluster2->EndWire());
                            
                            if(deltaWire2 >= deltaWire) continue;
                            if((fabs(cluster->EndWire() - cluster2->StartWire()) <= fMaximumDistance) && 
                                (fabs(cluster->EndTick() - cluster2->StartTick()) <= fMaximumTime))
                            {
                                if(cluster2->StartCharge()>cluster->EndCharge() || deltaWire > fMaximumLength)
                                {
                                    if(cluster2->StartWire() < cluster2->EndWire()) openangle=fabs(cluster->EndAngle() - cluster2->StartAngle());
                                    else
                                    {
                                        if(cluster2->StartAngle()<0) openangle=fabs(cluster->EndAngle() - ((-1.)*TMath::Pi()+cluster2->StartAngle()));
                                        else openangle=fabs(cluster->EndAngle() - (TMath::Pi()+cluster2->StartAngle()));
                                    }
                                    if(openangle>0.1 && openangle<1.57)
                                    {
                                        clusterPtrVec.push_back(cluster2);
                                    }
                                }
                            }
                            else if((fabs(cluster->EndWire() - cluster2->EndWire()) <= fMaximumDistance) && 
                                    (fabs(cluster->EndTick() - cluster2->EndTick()) <= fMaximumTime))
                            {
                                if(cluster2->EndCharge()>cluster->EndCharge() || deltaWire > fMaximumLength) 
                                {   
                                    if(cluster2->StartWire()<cluster2->EndWire())
                                    {
                                        if(cluster2->EndAngle()<0) openangle=fabs(cluster->EndAngle() -((-1.)*TMath::Pi()+cluster2->EndAngle()));
                                        else openangle=fabs(cluster->EndAngle() -(TMath::Pi()+cluster2->EndAngle()));
                                    }
                                    else openangle=fabs(cluster->EndAngle() - cluster2->EndAngle());
                                    if(openangle>0.1 && openangle<1.57)
                                    {
                                        clusterPtrVec.push_back(cluster2);
                                    }
                                }
                            }
                        }
                    } //end endwire

                    
                    // Handle the associations as everything related to the original cosmic tag for now
                    // I think we will need a better way going forward...
                    if(clusterPtrVec.size() > 1)
                    {
                	std::vector<art::Ptr<anab::CosmicTag> > cosmicVec = cosmicAssns.at(cluster.key());
                	if (!cosmicVec.empty())
                	{
                            art::Ptr<anab::CosmicTag>& cosmicTag(cosmicVec.front());
                            util::CreateAssn(*fMyProducerModule, event, cosmicTag, clusterPtrVec, *cosmicClusterAssociations);
                	                                    
				
			}
                    }
                }//end loop over good clusters
            }//end if num good clusters to proceed
        }//end make sure valid handles 
    }//end require valid handle
    
    // Add clusters and associations to event.
    event.put(std::move(cosmicClusterAssociations));
    //event.put(std::move(clusterPtrVec));
    
    return true;
}

} // namespace
