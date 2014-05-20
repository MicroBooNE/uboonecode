#ifndef RANDOMSERVER_CXX
#define RANDOMSERVER_CXX

#include "RandomServer.h"

namespace opdet {

  RandomServer* RandomServer::_me = 0;

  //------------------------------------------------------
  RandomServer::RandomServer() : fFlatRandom    (nullptr),
				 fPoissonRandom (nullptr),
				 fGausRandom    (nullptr)
  //------------------------------------------------------
  {
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();

    fFlatRandom    = new CLHEP::RandFlat    (engine);
    fPoissonRandom = new CLHEP::RandPoisson (engine);
    fGausRandom    = new CLHEP::RandGaussQ  (engine);

  }

  //---------------------------
  RandomServer::~RandomServer()
  //---------------------------
  {
    delete fFlatRandom;
    delete fPoissonRandom;
    delete fGausRandom;

    fFlatRandom    = nullptr;
    fPoissonRandom = nullptr;
    fGausRandom    = nullptr;

  }


  //--------------------------------------------------
  double RandomServer::Gaus(double mean, double sigma)
  //--------------------------------------------------
  {
    return fGausRandom->fire(mean,sigma);
  }

  //--------------------------------------------------
  double RandomServer::Uniform(double min, double max)
  //--------------------------------------------------
  {
    return fFlatRandom->fire(min,max);
  }

  //------------------------------------
  int RandomServer::Poisson(double mean)
  //------------------------------------
  {
    return fPoissonRandom->fire(mean);
  }

}
#endif
