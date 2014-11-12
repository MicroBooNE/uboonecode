////////////////////////////////////////////////////////////////////////
// Class:       MCShowerFinder
// Module Type: producer
// File:        MCShowerFinder_module.cc
//
// Generated at Mon Aug 11 05:40:00 2014 by Kazuhiro Terao using artmod
// from cetpkgsupport v1_05_04.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "MCShowerRecoAlg.h"

#include <memory>

class MCShowerFinder;

class MCShowerFinder : public art::EDProducer {
public:
  explicit MCShowerFinder(fhicl::ParameterSet const & p);
  virtual ~MCShowerFinder();

  void produce(art::Event & e) override;


private:

  // Declare member data here.
  ::sim::MCShowerRecoAlg fMCShowerAlg;
};

MCShowerFinder::MCShowerFinder(fhicl::ParameterSet const & pset)
  : fMCShowerAlg(pset.get< fhicl::ParameterSet >("MCShowerRecoAlg"))
{
  produces< std::vector< sim::MCShower> >();
  // Call appropriate produces<>() functions here.
}

MCShowerFinder::~MCShowerFinder()
{
  // Clean up dynamic memory and other resources here.
}

void MCShowerFinder::produce(art::Event & evt)
{
  std::unique_ptr< std::vector<sim::MCShower> > outShowerArray(new std::vector<sim::MCShower>);

  fMCShowerAlg.RecoMCShower(evt);

  for(size_t shower_index=0; shower_index < fMCShowerAlg.NumShowers(); ++shower_index)

    outShowerArray->push_back(fMCShowerAlg.ShowerProfile(shower_index));
    
  evt.put(std::move(outShowerArray));

}

DEFINE_ART_MODULE(MCShowerFinder)
