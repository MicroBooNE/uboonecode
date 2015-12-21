/**
 *  @file   Cluster2DNuAlg.h
 * 
 *  @brief  This is an algorithm for finding neutrino candidates using tracks and vertices
 * 
 */
#ifndef Cluster2DNuAlg_h
#define Cluster2DNuAlg_h

#include "uboone/TPCNeutrinoIDFilter/NeutrinoIDAlgBase.h"

// LArSoft includes
#include "larcore/Geometry/GeometryCore.h"
#include "lardata/DetectorInfo/DetectorProperties.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace neutrinoid
{

/**
 *  @brief  Cluster2DNuAlg class
 */
class Cluster2DNuAlg : virtual public NeutrinoIDAlgBase
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
    Cluster2DNuAlg(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~Cluster2DNuAlg();
    
    /**
     *  @brief a handler for the case where the algorithm control parameters are to be reset
     */
    virtual void reconfigure(fhicl::ParameterSet const&);
    
    /**
     *  @brief Set up for "beginJob" phase if requested
     */
    virtual void beginJob(art::ServiceHandle<art::TFileService>&);
    
    /**
     *  @brief Each algorithm may have different objects it wants "produced" so use this to
     *         let the top level producer module "know" what it is outputting
     */
    virtual void produces(art::EDProducer*);

    /**
     *  @brief Given the list of hits this will search for candidate Seed objects and return them
     */
    virtual bool findNeutrinoCandidates(art::Event&) const;

private:
    
    /**
     *  @ brief FHICL parameters.
     */
    std::string                fClusterModuleLabel;      ///< Producer of input tracks
    std::string                fCosmicModuleLabel;       ///< Producer of cosmic track tags
    // Below are cuts to apply, more or less in order
    size_t                     fPlaneToCheck;            ///< The plane to be analyzed
    size_t                     fMinimumHits;             ///< Minimum required hits in cluster
    float                      fMaximumTick;             ///< Maximum tick for cluster
    float                      fMinimumTick;             ///< Minimum tick for cluster
    float                      fMaximumWire;             ///< Maximum allowed wire in cluster
    float                      fMinimumWire;             ///< Minimum allowed wire in cluster
    float                      fMaximumAngle;            ///< Maximum allowed angle of long clusters
    float                      fMaximumLengthCut;        ///< Maximum length of cluster (in wires)
    float                      fMaximumMatchedLengthCut; ///< Maximum length of matched long cluster (in wires)
    float                      fMaximumLength;           ///< Minimum length of cluster (in wires)
    float                      fMinimumDeltaTicks;       ///< Minimum time for cluster (in ticks)
    float                      fMinimumDeltaWires;       ///< Minimum wire for cluster (in wires)
    float                      fMaxCosmicScore;          ///< Maximum allowed cosmic score
    size_t                     fMinCandidateClusters;    ///< Minimum candidate clusters to proceed
    float                      fMaximumDistance;         ///< Maximum distance between 2 clusters
    float                      fMaximumTime;             ///< Maximum time beteween 2 clusters
    
    art::EDProducer*           fMyProducerModule;        ///< The producer module driving us
    
    /// @{
    /**
     *  @brief Standard useful properties
     */
    geo::GeometryCore const*             m_geometry;           ///< pointer to the Geometry service
    detinfo::DetectorProperties const* m_detector;           ///< Pointer to the detector properties
    /// @}
};

} // namespace lar_cluster3d
#endif
