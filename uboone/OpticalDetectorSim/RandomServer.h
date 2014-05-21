/**
 * \file RandomServer.h
 *
 * \ingroup OpticalDetector
 * 
 * \brief Class def header for a class RandomServer
 *
 * @author kazuhiro
 */

/** \addtogroup OpticalDetectorSim

    @{*/
#ifndef RANDOMSERVER_H
#define RANDOMSERVER_H

#include <iostream>
#include <TRandom3.h>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
//#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandGaussQ.h"

namespace opdet {
  /**
     \class RandomServer
     Simple class that interfaces with different random number generation scheme which may depend on the framework
  */
  class RandomServer : public TRandom3 {
    
  private:
    
    /// Default constructor
    RandomServer(){};
    
    /// Default destructor
    virtual ~RandomServer(){};
    
    static RandomServer* _me;
    
  public:
    
    static RandomServer& GetME()
    {
      
      if(!_me) _me = new RandomServer;
      
      return *_me;
    }
    
  public:
    
    /// Default constructor
    //RandomServer(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    
    /// Default destructor
    //virtual ~RandomServer();

    //void reconfigure(fhicl::ParameterSet const& pset);

    /*
    double Gaus(double mean, double sigma);

    double Uniform(double min=0, double max=1);

    int Poisson(double mean);
    */
  protected:

    /*
    ::CLHEP::RandFlat*    fFlatRandom;
    ::CLHEP::RandPoisson* fPoissonRandom;
    ::CLHEP::RandGaussQ*  fGausRandom;
    */
  };
}

//DECLARE_ART_SERVICE(opdet::RandomServer, LEGACY)

#endif
/** @} */ // end of doxygen group 

