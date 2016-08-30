/**
 *  @file    NeutrinoIDAlgFactory.cxx
 * 
 *  @brief   A small module to instantiate various algorithms inheriting from NeutrinoIDAlgBase
 * 
 *  @authors usher@slac.stanford.edu
 */

// Includes for algorithms
#include "uboone/TPCNeutrinoIDFilter/Algorithms/NeutrinoIDAlgFactory.h"
#include "uboone/TPCNeutrinoIDFilter/Algorithms/NuMuCCInclusiveAlg.h"
#include "uboone/TPCNeutrinoIDFilter/Algorithms/ModNuMuCCInclusiveAlg.h"
#include "uboone/TPCNeutrinoIDFilter/Algorithms/AltNuMuCCInclusiveAlg.h"
#include "uboone/TPCNeutrinoIDFilter/Algorithms/TrackPairPlusVertexAlg.h"
#include "uboone/TPCNeutrinoIDFilter/Algorithms/Cluster2DNuAlg.h"
#include "uboone/TPCNeutrinoIDFilter/Algorithms/NuMuCCSelectionIIAlg.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace neutrinoid
{
    
std::unique_ptr< NeutrinoIDAlgBase > NeutrinoIDAlgFactory::MakeNeutrinoIDAlg(fhicl::ParameterSet const& p)
{
    std::string algName = p.get<std::string>("NeutrinoIDAlgName");
    
    std::unique_ptr< NeutrinoIDAlgBase > ptr;
    
    if     (algName.compare("NuMuCCInclusiveAlg")==0)
    {
        std::unique_ptr< NeutrinoIDAlgBase > new_ptr(new NuMuCCInclusiveAlg(p));
        ptr.swap(new_ptr);
    }
    else if(algName.compare("ModNuMuCCInclusiveAlg")==0)
    {
        std::unique_ptr< NeutrinoIDAlgBase > new_ptr(new ModNuMuCCInclusiveAlg(p));
        ptr.swap(new_ptr);
    }
    else if(algName.compare("AltNuMuCCInclusiveAlg")==0)
    {
        std::unique_ptr< NeutrinoIDAlgBase > new_ptr(new AltNuMuCCInclusiveAlg(p));
        ptr.swap(new_ptr);
    }
    else if(algName.compare("TrackPairPlusVertexAlg")==0)
    {
        std::unique_ptr< NeutrinoIDAlgBase > new_ptr(new TrackPairPlusVertexAlg(p));
        ptr.swap(new_ptr);
    }
    else if(algName.compare("Cluster2DNuAlg")==0)
    {
        std::unique_ptr< NeutrinoIDAlgBase > new_ptr(new Cluster2DNuAlg(p));
        ptr.swap(new_ptr);
    }
    else if(algName.compare("NuMuCCSelectionIIAlg")==0){
        std::unique_ptr< NeutrinoIDAlgBase > new_ptr(new NuMuCCSelectionIIAlg(p));
        ptr.swap(new_ptr);
    }      
    else{
        std::cout << "Algname is ... " << algName << std::endl;
        throw std::runtime_error("ERROR in NeutrinoIDAlgFactory: No registered Neutrino ID with that name.");
    }
    
    return std::move(ptr);
}

} // namespace neutrinoid
