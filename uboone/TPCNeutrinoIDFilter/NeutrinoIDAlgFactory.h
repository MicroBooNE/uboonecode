/**
 *  @file   NeutrinoIDAlgBase.h
 * 
 *  @brief  This is intended to define an interface to all Seed finder algorithms employed
 *          by the 3D clustering
 * 
 */
#ifndef NeutrinoIDAlgFactory_h
#define NeutrinoIDAlgFActory_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"

// Abstract class include
#include "TPCNeutrinoIDFilter/NeutrinoIDAlgBase.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace neutrinoid
{
/**
 *  @brief  NeutrinoIDAlgFactory class
 */
class NeutrinoIDAlgFactory
{
public:
    NeutrinoIDAlgFactory() {}
   ~NeutrinoIDAlgFactory() {}
    
    std::unique_ptr< NeutrinoIDAlgBase > MakeNeutrinoIDAlg(fhicl::ParameterSet const& p);
private:
};

} // namespace lar_cluster3d
#endif
