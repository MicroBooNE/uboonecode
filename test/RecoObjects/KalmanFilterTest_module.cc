//
// Name:  KalmanFilterTest_module.cc
//
// Purpose: KalmanFilterTest module.
//
// Created:  12-Apr-2012  H. Greenlee

#include <iostream>
#include <cassert>

#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDAnalyzer.h"

#include "RecoObjects/KHitWireX.h"
#include "RecoObjects/SurfYZPlane.h"
#include "RecoObjects/PropYZPlane.h"
#include "RecoObjects/KHitMulti.h"

namespace trkf
{
  class KalmanFilterTest : public art::EDAnalyzer
  {
  public:

    // Constructor, destructor.

    explicit KalmanFilterTest(fhicl::ParameterSet const& pset);
    ~KalmanFilterTest();

    // Overrides.

    void beginJob();
    void analyze(const art::Event& evt);

  private:

    // Attributes.

  };

  DEFINE_ART_MODULE(KalmanFilterTest)

  KalmanFilterTest::KalmanFilterTest(const fhicl::ParameterSet& pset)
  : EDAnalyzer(pset)
  {}

  void KalmanFilterTest::beginJob()
  {
    // Make sure assert is enabled.

    bool assert_flag = false;
    assert((assert_flag = true, assert_flag));
    if ( ! assert_flag ) {
      std::cerr << "Assert is disabled" << std::endl;
      abort();
    }

    // Make a test track.

    std::shared_ptr<const trkf::Surface> psurf(new trkf::SurfYZPlane(0., 0., 1000., 0.));
    TrackVector vec(5);
    vec(0) = 10.;
    vec(1) = 0.;
    vec(2) = 0.5;
    vec(3) = 0.7;
    vec(4) = 1.;
    TrackError err(5);
    err.clear();
    err(0, 0) = 1000.;
    err(1, 1) = 1000.;
    err(2, 2) = 10.;
    err(3, 3) = 10.;
    err(4, 4) = 10.;
    KETrack tre(psurf, vec, err, trkf::Surface::FORWARD, 13);

    PropYZPlane prop(100., true);

    // Make some test measurements.

    int nsurf = 20;
    std::vector<std::shared_ptr<const trkf::KHitBase> > phits;
    for(int i=0; i<nsurf; ++i) {
      int channel = 0;
      if(i%3 == 0)
	channel = 1000 + 100*i;
      if(i%3 == 1)
	channel = 4000 + 100*i;
      else
	channel = 6000 + 100*i;
      trkf::KHitWireX hit(channel, 0., 0.);

      // Propagate track to measurement surface.

      KETrack treprop(tre);
      boost::optional<double> ok = prop.vec_prop(treprop, hit.getMeasSurface(), trkf::Propagator::UNKNOWN, false);
      assert(!!ok);
      double x = treprop.getVector()(0);
      phits.push_back(std::shared_ptr<const trkf::KHitBase>(new KHitWireX(channel, x, 0.1)));
    }

    // Make a new starting track.

    TrackVector vec2(5);
    vec2(0) = 0.;
    vec2(1) = 0.;
    vec2(2) = 0.4;
    vec2(3) = 0.6;
    vec2(4) = 1.;
    TrackError err2(5);
    err2.clear();
    err2(0, 0) = 1000.;
    err2(1, 1) = 1000.;
    err2(2, 2) = 10.;
    err2(3, 3) = 10.;
    err2(4, 4) = 10.;
    KETrack tre2(psurf, vec2, err2, trkf::Surface::FORWARD, 13);

    // Loop over measurements.

    for(unsigned int i=0; i<phits.size(); ++i) {
      const KHitBase& hit = *(phits[i]);
      std::cout << "\nMeasurement " << i << std::endl;
      std::cout << "Original track:" << std::endl;
      std::cout << tre2;
      hit.predict(tre2, &prop);
      std::cout << "Hit after prediction:" << std::endl;
      std::cout << hit << std::endl;
      hit.update(tre2);
      std::cout << "Track after update:" << std::endl;
      std::cout << tre2 << std::endl;
    }

    // Make another new starting track.

    KETrack tre3(psurf, vec2, err2, trkf::Surface::FORWARD, 13);

    // Add all measurements into a single composit measurement.

    std::cout << "\nCompound hit method." << std::endl;
    KHitMulti mhit;
    for(unsigned int i=0; i<phits.size(); ++i) {
      std::shared_ptr<const KHitBase>& phit = phits[i];
      mhit.addMeas(phit);
    }
    mhit.predict(tre3, &prop);
    mhit.update(tre3);
    std::cout << "Track after update:" << std::endl;
    std::cout << tre3 << std::endl;

    // Try again.

    tre3.setError(err2);
    mhit.predict(tre3, &prop);
    mhit.update(tre3);
    std::cout << "Track after second update:" << std::endl;
    std::cout << tre3 << std::endl;

    // Try again.

    tre3.setError(err2);
    mhit.predict(tre3, &prop);
    mhit.update(tre3);
    std::cout << "Track after third update:" << std::endl;
    std::cout << tre3 << std::endl;
    

    // Done (success).

    std::cout << "KalmanFilterTest: All tests passed." << std::endl;
  }

  KalmanFilterTest::~KalmanFilterTest()
  {}

  void KalmanFilterTest::analyze(const art::Event& /* evt */)
  {}
}
