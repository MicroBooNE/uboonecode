/**
 *  @file   TrackPairPlusVertexAlg.h
 * 
 *  @brief  This is an algorithm for finding neutrino candidates using tracks and vertices
 * 
 */
#ifndef TrackPairPlusVertexAlg_h
#define TrackPairPlusVertexAlg_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"

#include "TPCNeutrinoIDFilter/NeutrinoIDAlgBase.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Utilities/DetectorProperties.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace neutrinoid
{

/**
 *  @brief  TrackPairPlusVertexAlg class
 */
class TrackPairPlusVertexAlg : virtual public NeutrinoIDAlgBase
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
    TrackPairPlusVertexAlg(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~TrackPairPlusVertexAlg();
    
    /**
     *  @brief a handler for the case where the algorithm control parameters are to be reset
     */
    virtual void reconfigure(fhicl::ParameterSet const&);
    
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
    std::string                fTrackModuleLabel;        ///< Producer of input tracks
    std::string                fVertexModuleLabel;       ///< Producer of input vertices
    std::string                fCosmicModuleLabel;       ///< Producer of cosmic track tags
    double                     fCosmicScoreCut;          ///< Cut value for possible cosmic tag scores
    double                     fNeutrinoVtxTrackDistCut; ///< Cut to select neutrino candidate
    
    art::EDProducer*           fMyProducerModule;        ///< The producer module driving us
    
    /**
     *  @brief Standard useful properties
     */
    geo::Geometry*             m_geometry;            //< pointer to the Geometry service
    util::DetectorProperties*  m_detector;            //< Pointer to the detector properties
};

} // namespace lar_cluster3d
#endif
