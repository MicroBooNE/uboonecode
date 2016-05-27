/**
 *  @file   NuMuCCInclusiveAlg.h
 * 
 *  @brief  This is an algorithm for finding neutrino candidates using tracks and vertices
 * 
 */
#ifndef ModNuMuCCInclusiveAlg_h
#define ModNuMuCCInclusiveAlg_h

#include "uboone/TPCNeutrinoIDFilter/Algorithms/NeutrinoIDAlgBase.h"

// LArSoft includes
#include "larcore/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// Root includes
#include "TH1D.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace neutrinoid
{

/**
 *  @brief  NuMuCCInclusiveAlg class
 */
class ModNuMuCCInclusiveAlg : virtual public NeutrinoIDAlgBase
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
    ModNuMuCCInclusiveAlg(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~ModNuMuCCInclusiveAlg();
    
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
    
    bool   inFV(double x, double y, double z) const;
    
    double FlashTrackDist(double flash, double start, double end) const;
    
    /**
     *  @ brief FHICL parameters.
     */
    std::string                fTrackModuleLabel;        ///< Producer of input tracks
    std::string                fVertexModuleLabel;       ///< Producer of input vertices
    std::string                fOpFlashModuleLabel;      ///< Producer of flashes
    
    double                     fDistToEdgeX;             ///< fiducial volume - x
    double                     fDistToEdgeY;             ///< fiducial volume - y
    double                     fDistToEdgeZ;             ///< fiducial volume - z
    
    double                     fFlashWidth;              ///< Cut on flash width
    double                     fBeamMin;                 ///< Cut on min beam time
    double                     fBeamMax;                 ///< Cut on max beam time
    double                     fPEThresh;                ///< Cut on PE threshold
    double                     fMinTrk2VtxDist;          ///< Minimum track to vertex distance
    double                     fMinTrackLen;             ///< Minimum track length
    bool                       fDoHists;                 ///< Fill histograms
    
    TH1D*                      fNFlashPerEvent;          ///< number of flashes per event
    TH1D*                      fFlashPE;                 ///< flash photoelectrons
    TH1D*                      fFlashTime;               ///< flash timing
    
    art::EDProducer*           fMyProducerModule;        ///< The producer module driving us
    
    /// @{
    /**
     *  @brief Standard useful properties
     */
    geo::GeometryCore const*            fGeometry;           ///< pointer to the Geometry service
    detinfo::DetectorProperties const*  fDetector;           ///< Pointer to the detector properties
    /// @}
};

} // namespace lar_cluster3d
#endif
