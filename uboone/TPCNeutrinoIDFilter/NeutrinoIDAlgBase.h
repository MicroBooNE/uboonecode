/**
 *  @file    NeutrinoIDAlgBase.h
 * 
 *  @brief   This is intended to define an interface to all Seed finder algorithms employed
 *           by the 3D clustering
 *
 *  @authors usher@slac.stanford.edu
 * 
 */
#ifndef NeutrinoIDAlgBase_h
#define NeutrinoIDAlgBase_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"

//------------------------------------------------------------------------------------------------------------------------------------------

// foward declaration for producers
namespace art
{
    class EDProducer;
    class Event;
}

namespace neutrinoid
{
/**
 *  @brief  NeutrinoIDAlgBase class
 */
class NeutrinoIDAlgBase
{
public:
    /**
     *  @brief Require that a handler is definied in case the algorithm control parameters are to be reset
     */
    virtual void reconfigure(fhicl::ParameterSet const &pset) = 0;

    /**
     *  @brief Each algorithm may have different objects it wants "produced" so use this to 
     *         let the top level producer module "know" what it is outputting
     */
    virtual void produces(art::EDProducer*) = 0;
    
    /**
     *  @brief Define the interface to take an input list of 3D hits and return seed candidates
     *         so hits are ordered along the axis
     */
    virtual bool findNeutrinoCandidates(art::Event & event) const = 0;

protected:

    /**
     *  @brief Define a comparator which will sort hits by arc length along a PCA axis
     */
//    struct Sort3DHitsByArcLen3D
//    {
//        bool operator()(const reco::ClusterHit3D* left, const reco::ClusterHit3D* right)
//        {
//            return left->getArclenToPoca() < right->getArclenToPoca();
//        }
//
//    };
    
private:
};

} // namespace lar_cluster3d
#endif
