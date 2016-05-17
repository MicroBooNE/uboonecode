/**
 *  @file   AltNuMuCCInclusiveAlg.h
 * 
 *  @brief  This is an algorithm for finding neutrino candidates using tracks and vertices
 * 
 */
#ifndef AltNuMuCCInclusiveAlg_h
#define AltNuMuCCInclusiveAlg_h

#include "uboone/TPCNeutrinoIDFilter/Algorithms/NeutrinoIDAlgBase.h"

// Framework Includes
#include "art/Framework/Core/FindManyP.h"

// LArSoft includes
#include "larcore/Geometry/GeometryCore.h"
#include "lardata/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfo/DetectorClocks.h"

// Root includes
#include "TH1D.h"
#include "TH2D.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace recob
{
    class OpFlash;
    class PFParticle;
    class Track;
    class Vertex;
}

namespace neutrinoid
{

/**
 *  @brief  AltNuMuCCInclusiveAlg class
 */
class AltNuMuCCInclusiveAlg : virtual public NeutrinoIDAlgBase
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
    AltNuMuCCInclusiveAlg(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~AltNuMuCCInclusiveAlg();
    
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
    
    double projectedLength(const recob::Track* track) const;
    bool inFV(const TVector3&) const;
    bool endPointOK(const TVector3& pos) const;
    
    int traversePFParticleHierarchy(art::Handle<std::vector<recob::PFParticle>>& pfParticleHandle,
                                    size_t                                       pfParticleIdx,
                                    const art::FindManyP<recob::Track>&          trackAssns,
                                    const art::FindManyP<recob::Vertex>&         vertexAssns,
                                    std::vector<art::Ptr<recob::PFParticle>>&    pfParticleVec,
                                    std::vector<art::Ptr<recob::Track>>&         trackVec,
                                    std::vector<art::Ptr<recob::Vertex>>&        vertexVec) const;
    
    void   getBestFlashTrackDist(const std::vector<art::Ptr<recob::OpFlash>>&, double, double, size_t&, double&) const;
    double FlashTrackDistInZ(double flash, double start, double end) const;
    void   getTrackVertexDCA(const TVector3& vertex, const TVector3& trackStart, TVector3& trackDir, double& doca, double& arcLen) const;
    
    /**
     *  @ brief FHICL parameters.
     */
    std::string                fPFParticleModuleLabel;   ///< Producer of input PFParticles
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
    double                     fMinTrackLen;             ///< Minimum track length
    double                     fMaxTrackDoca;            ///< Maximum track to vertex doca
    double                     fMaxTrackArcLen;          ///< Maximum track arclen to trk/vtx doca
    bool                       fDoHists;                 ///< Fill histograms
    
    std::vector<TH1D*>         fTH1DVec;                 ///< histogram container
    std::vector<TH2D*>         fTH2DVec;                 ///< 2D histogram container
    
    art::EDProducer*           fMyProducerModule;        ///< The producer module driving us
    
    /// @{
    /**
     *  @brief Standard useful properties
     */
    geo::GeometryCore const*            fGeometry;           ///< pointer to the Geometry service
    detinfo::DetectorProperties const*  fDetector;           ///< Pointer to the detector properties
    detinfo::DetectorClocks const*      fClocks;             ///< Pointer to the clock services
    /// @}
};

} // namespace lar_cluster3d
#endif
